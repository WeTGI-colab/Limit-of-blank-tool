#!/usr/bin/env python3
"""Combine the per-sample pile-up tables (from the SLURM array) into the cohort table.

Reads every results/pileup/*.tsv (filtered) and results/pileup_raw/*.tsv (raw) and writes
results/cohort_alt.tsv, identical in format to the serial run_pileup.py output. Run after the
per-sample array finishes; then run scripts/run_lob.py and scripts/plot_regions.py.
"""
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from amplicon_lob import tables  # noqa: E402

DATA = REPO / "data"
RESULTS = REPO / "results"


def load_frames(directory):
    return [pd.read_csv(f, sep="\t") for f in sorted(directory.glob("*.tsv"))]


def load_truth():
    gt = DATA / "synthetic" / "ground_truth.tsv"
    truth = {}
    if gt.exists():
        for line in gt.read_text().splitlines()[1:]:
            f = line.split("\t")
            truth[(f[2], int(f[3]), f[5])] = f[0]
    return truth


def main():
    raw = load_frames(RESULTS / "pileup_raw")
    filt = load_frames(RESULTS / "pileup")
    if not raw or not filt:
        sys.exit("No per-sample tables found in results/pileup[_raw]/ -- run the array first.")
    if len(raw) != len(filt):
        print(f"WARNING: {len(raw)} raw vs {len(filt)} filt per-sample tables "
              f"(some tasks may have failed).")

    agg = tables.merge_regimes(tables.aggregate(raw), tables.aggregate(filt))
    truth = load_truth()
    agg["truth"] = [truth.get((c, p, a), "") for c, p, a in
                    zip(agg["chrom"], agg["pos"], agg["alt"])]
    agg.to_csv(RESULTS / "cohort_alt.tsv", sep="\t", index=False)
    print(f"Aggregated {len(filt)} samples -> results/cohort_alt.tsv "
          f"({len(agg)} position x alt rows)")


if __name__ == "__main__":
    main()
