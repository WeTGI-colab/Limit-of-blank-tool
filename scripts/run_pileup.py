#!/usr/bin/env python3
"""Build per-base pile-up tables for the synthetic cohort under two QC regimes.

Two read-quality filter sets are applied, mirroring the laboratory:
  raw   -- minbq=1,  minmq=1   : essentially unfiltered, i.e. all sequencer noise.
  filt  -- minbq=20, minmq=30  : the Pisces calling thresholds, i.e. the noise that survives
                                  the laboratory's quality filters.

Produces under results/:
  pileup/<sample>.tsv   per sample, the *filt* (callable) counts -> input to the LoB model.
  cohort_alt.tsv        one row per (position, alternate base) with across-sample statistics
                        for BOTH regimes (columns suffixed _raw / _filt) -> input to the plots.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from amplicon_lob import panel, pileup, samples  # noqa: E402

DATA = REPO / "data"
MANIFEST = DATA / "manifest.tsv"
RESULTS = REPO / "results"

FILTERS = {"raw": (1, 1), "filt": (20, 30)}     # name -> (min_base_quality, min_mapping_quality)


def load_truth():
    """Ground truth is only available for the synthetic cohort; absent for production data."""
    gt = DATA / "synthetic" / "ground_truth.tsv"
    truth = {}
    if gt.exists():
        for line in gt.read_text().splitlines()[1:]:
            f = line.split("\t")
            truth[(f[2], int(f[3]), f[5])] = f[0]
    return truth


def per_sample_table(name, counts, refbases, genes):
    rows = []
    for (chrom, pos0), c in sorted(counts.items()):
        pos1 = pos0 + 1
        ref = refbases.get((chrom, pos1))
        if ref is None:
            continue
        totals = {b: c[b][0] + c[b][1] for b in "ACGT"}
        depth = sum(totals.values())
        nonref = depth - totals[ref]
        row = {"sample": name, "gene": genes[chrom], "chrom": chrom, "pos": pos1,
               "ref": ref, "depth": depth,
               "A": totals["A"], "C": totals["C"], "G": totals["G"], "T": totals["T"],
               "nonref": nonref, "nonref_vaf": nonref / depth if depth else 0.0}
        for b in "ACGT":
            row[f"{b}_fwd"], row[f"{b}_rev"] = c[b][0], c[b][1]
        rows.append(row)
    return pd.DataFrame(rows)


def aggregate(frames):
    """Across-sample aggregate at (position, alternate base) resolution."""
    long = []
    for df in frames:
        for _, r in df.iterrows():
            for alt in "ACGT":
                if alt == r["ref"] or r["depth"] == 0:
                    continue
                ac = r[f"{alt}_fwd"] + r[f"{alt}_rev"]
                long.append({"chrom": r["chrom"], "pos": r["pos"], "gene": r["gene"],
                             "ref": r["ref"], "alt": alt, "depth": r["depth"],
                             "alt_fwd": r[f"{alt}_fwd"], "alt_rev": r[f"{alt}_rev"],
                             "vaf": ac / r["depth"]})
    L = pd.DataFrame(long)

    def summarise(g):
        fwd, rev = g["alt_fwd"].sum(), g["alt_rev"].sum()
        return pd.Series({
            "n_samples": len(g),
            "mean_vaf": g["vaf"].mean(),
            "alt_reads_mean": (fwd + rev) / len(g),
            "mean_depth": g["depth"].mean(),
            "strand_frac_fwd": fwd / (fwd + rev) if (fwd + rev) else np.nan,
        })
    return L.groupby(["chrom", "pos", "gene", "ref", "alt"]).apply(
        summarise, include_groups=False).reset_index()


def main():
    amps = panel.load_amplicons()
    regions = panel.thick_regions(amps)
    refbases = panel.load_reference_bases()
    genes = panel.gene_by_chrom(amps)
    truth = load_truth()
    manifest = samples.load_manifest(MANIFEST)

    (RESULTS / "pileup").mkdir(parents=True, exist_ok=True)
    frames = {k: [] for k in FILTERS}
    for row in manifest:
        name = row["sample"]
        for key, (bq, mq) in FILTERS.items():
            counts = pileup.pileup_bam(row["bam"], regions, min_bq=bq, min_mq=mq)
            frames[key].append(per_sample_table(name, counts, refbases, genes))
        frames["filt"][-1].to_csv(RESULTS / "pileup" / f"{name}.tsv", sep="\t", index=False)
    print(f"Wrote {len(manifest)} per-sample (filt) tables to results/pileup/")

    keys = ["chrom", "pos", "gene", "ref", "alt"]
    agg = aggregate(frames["raw"]).merge(aggregate(frames["filt"]), on=keys,
                                         how="outer", suffixes=("_raw", "_filt"))
    numeric = [c for c in agg.columns if c not in keys]
    agg[numeric] = agg[numeric].fillna(0)
    agg["truth"] = [truth.get((c, p, a), "") for c, p, a in
                    zip(agg["chrom"], agg["pos"], agg["alt"])]
    agg = agg.sort_values(["chrom", "pos", "alt"])
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
