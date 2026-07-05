#!/usr/bin/env python3
"""Pile-up for a SINGLE sample -- the unit of work for one SLURM array task.

Selects one manifest row (by 1-based --index, e.g. $SLURM_ARRAY_TASK_ID, or by --sample name),
runs the pile-up under both quality regimes, and writes:
  results/pileup/<sample>.tsv       filtered (minbq20/minmq30) -> LoB model
  results/pileup_raw/<sample>.tsv   raw (minbq1/minmq1)        -> cohort aggregate

Existing outputs are skipped unless --force, so a re-submitted array only redoes failed tasks.
After all tasks finish, combine with scripts/aggregate_cohort.py.
"""
import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from amplicon_lob import panel, pileup, samples, tables  # noqa: E402

MANIFEST = REPO / "data" / "manifest.tsv"
RESULTS = REPO / "results"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--index", type=int, help="1-based manifest row (e.g. $SLURM_ARRAY_TASK_ID)")
    g.add_argument("--sample", help="sample id from the manifest")
    ap.add_argument("--force", action="store_true", help="recompute even if outputs exist")
    args = ap.parse_args()

    manifest = samples.load_manifest(MANIFEST)
    if args.index is not None:
        if not (1 <= args.index <= len(manifest)):
            sys.exit(f"--index {args.index} out of range 1..{len(manifest)}")
        row = manifest[args.index - 1]
    else:
        hits = [r for r in manifest if r["sample"] == args.sample]
        if not hits:
            sys.exit(f"sample {args.sample} not in manifest")
        row = hits[0]

    amps = panel.load_amplicons()
    regions = panel.thick_regions(amps)
    refbases = panel.load_reference_bases()
    genes = panel.gene_by_chrom(amps)
    name = row["sample"]

    for key, (bq, mq) in tables.FILTERS.items():
        outdir = RESULTS / ("pileup" if key == "filt" else "pileup_raw")
        outdir.mkdir(parents=True, exist_ok=True)
        out = outdir / f"{name}.tsv"
        if out.exists() and not args.force:
            print(f"skip (exists): {out.relative_to(REPO)}")
            continue
        counts = pileup.pileup_bam(row["bam"], regions, min_bq=bq, min_mq=mq)
        tables.per_sample_table(name, counts, refbases, genes).to_csv(out, sep="\t", index=False)
        print(f"wrote {out.relative_to(REPO)}")


if __name__ == "__main__":
    main()
