#!/usr/bin/env python3
"""Serial pile-up over the whole cohort (one process). For large cohorts on a cluster, use the
SLURM path instead: scripts/pileup_one_sample.py (one job per sample) + scripts/aggregate_cohort.py.

Two read-quality filter sets are applied (see amplicon_lob.tables.FILTERS): raw (minbq1/minmq1,
all sequencer noise) and filt (minbq20/minmq30, Pisces callable). Produces under results/:
  pileup/<sample>.tsv       per sample, filt counts -> input to the LoB model
  pileup_raw/<sample>.tsv   per sample, raw counts  -> for the cohort aggregate
  cohort_alt.tsv            per (position, alt) with both regimes (_raw / _filt) -> plots
"""
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from amplicon_lob import panel, pileup, samples, tables  # noqa: E402

DATA = REPO / "data"
MANIFEST = DATA / "manifest.tsv"
RESULTS = REPO / "results"


def load_truth():
    """Ground truth exists only for the synthetic cohort; absent for production data."""
    gt = DATA / "synthetic" / "ground_truth.tsv"
    truth = {}
    if gt.exists():
        for line in gt.read_text().splitlines()[1:]:
            f = line.split("\t")
            truth[(f[2], int(f[3]), f[5])] = f[0]
    return truth


def main():
    amps = panel.load_amplicons()
    regions = panel.thick_regions(amps)
    refbases = panel.load_reference_bases()
    genes = panel.gene_by_chrom(amps)
    manifest = samples.load_manifest(MANIFEST)

    (RESULTS / "pileup").mkdir(parents=True, exist_ok=True)
    (RESULTS / "pileup_raw").mkdir(parents=True, exist_ok=True)
    frames = {k: [] for k in tables.FILTERS}
    for row in manifest:
        name = row["sample"]
        for key, (bq, mq) in tables.FILTERS.items():
            counts = pileup.pileup_bam(row["bam"], regions, min_bq=bq, min_mq=mq)
            df = tables.per_sample_table(name, counts, refbases, genes)
            frames[key].append(df)
            outdir = "pileup" if key == "filt" else "pileup_raw"
            df.to_csv(RESULTS / outdir / f"{name}.tsv", sep="\t", index=False)
    print(f"Wrote {len(manifest)} per-sample tables to results/pileup[_raw]/")

    masks = samples.load_masks(manifest)     # remove patient variants -> cohort table is blank
    agg = tables.merge_regimes(tables.aggregate(frames["raw"], masks),
                               tables.aggregate(frames["filt"], masks))
    truth = load_truth()
    agg["truth"] = [truth.get((c, p, a), "") for c, p, a in
                    zip(agg["chrom"], agg["pos"], agg["alt"])]
    agg.to_csv(RESULTS / "cohort_alt.tsv", sep="\t", index=False)
    print(f"Wrote cohort_alt.tsv  ({len(agg)} position x alt rows, raw + filt)")

    print("\nTop non-reference signals (raw vs filtered mean VAF):")
    cols = ["gene", "chrom", "pos", "ref", "alt", "mean_vaf_raw", "mean_vaf_filt",
            "strand_frac_fwd_filt", "truth"]
    top = agg.sort_values("mean_vaf_raw", ascending=False).head(10)[cols]
    with pd.option_context("display.width", 200, "display.max_columns", None):
        print(top.to_string(index=False, formatters={
            "mean_vaf_raw": "{:.4f}".format, "mean_vaf_filt": "{:.4f}".format,
            "strand_frac_fwd_filt": "{:.2f}".format}))


if __name__ == "__main__":
    main()
