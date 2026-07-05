#!/usr/bin/env python3
"""Per-run breakdown of the flagged artefact sites, to expose run-level batch effects.

For each systematic-artefact site flagged in results/lob_table.tsv (or the top --top by blank
VAF), report the non-reference VAF across the samples of each sequencing run (grouped by the
manifest 'run' column). A uniform signal across runs points to an instrument-level artefact; a
spike in one or a few runs points to a run/batch effect. Run IDs are flowcell IDs, not
patient-identifiable, so the output is shareable.

Outputs:
  results/artefact_by_run.tsv   long: (site, run) -> n_samples, mean_vaf, max_vaf, mean_depth
  console pivot                 rows = site, cols = run, values = mean non-reference VAF (%)
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from amplicon_lob import samples  # noqa: E402

MANIFEST = REPO / "data" / "manifest.tsv"
RESULTS = REPO / "results"
PILEUP = RESULTS / "pileup"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top", type=int, default=0,
                    help="use the top-N flagged sites by blank VAF (0 = all flagged)")
    ap.add_argument("--min-depth", type=int, default=500,
                    help="ignore sites below this mean depth")
    args = ap.parse_args()

    run_of = {r["sample"]: r["run"] for r in samples.load_manifest(MANIFEST)}
    lob = pd.read_csv(RESULTS / "lob_table.tsv", sep="\t")
    sites = lob[lob["flagged"] == True]                                     # noqa: E712
    if "mean_depth" in sites.columns:
        sites = sites[sites["mean_depth"] >= args.min_depth]
    sites = sites.sort_values("blank_mean_vaf", ascending=False)
    if args.top:
        sites = sites.head(args.top)
    targets = [(r.chrom, int(r.pos), r.alt, r.sub) for r in sites.itertuples()]
    if not targets:
        sys.exit("No flagged sites to break down.")
    want_pos = {p for _, p, _, _ in targets}

    recs = []
    for tsv in sorted(PILEUP.glob("*.tsv")):
        sample = tsv.stem
        run = run_of.get(sample, "?")
        df = pd.read_csv(tsv, sep="\t")
        df = df[df["pos"].isin(want_pos)].set_index(["chrom", "pos"])
        for chrom, pos, alt, sub in targets:
            if (chrom, pos) not in df.index:
                continue
            row = df.loc[(chrom, pos)]
            depth = row["depth"]
            ac = row[f"{alt}_fwd"] + row[f"{alt}_rev"]
            recs.append({"site": f"{chrom}:{pos} {sub}", "run": run, "sample": sample,
                         "depth": depth, "vaf": ac / depth if depth else 0.0})
    R = pd.DataFrame(recs)

    g = R.groupby(["site", "run"]).agg(
        n_samples=("vaf", "size"), mean_vaf=("vaf", "mean"),
        max_vaf=("vaf", "max"), mean_depth=("depth", "mean")).reset_index()
    g.to_csv(RESULTS / "artefact_by_run.tsv", sep="\t", index=False)
    print(f"Wrote results/artefact_by_run.tsv  ({len(targets)} sites x {R['run'].nunique()} runs)")

    piv = (R.groupby(["site", "run"])["vaf"].mean() * 100).unstack().round(3)
    print("\nMean non-reference VAF (%) per run  [rows = artefact site, cols = run]:")
    with pd.option_context("display.width", 250, "display.max_columns", None):
        print(piv.to_string())


if __name__ == "__main__":
    main()
