"""Shared per-sample and cohort table construction (used by the serial and the SLURM paths)."""
import numpy as np
import pandas as pd

FILTERS = {"raw": (1, 1), "filt": (20, 30)}     # name -> (min_base_quality, min_mapping_quality)


def per_sample_table(name, counts, refbases, genes):
    """One row per genomic position: depth, per-strand depth, and per-base strand-split counts."""
    rows = []
    for (chrom, pos0), c in sorted(counts.items()):
        pos1 = pos0 + 1
        ref = refbases.get((chrom, pos1))
        if ref is None:
            continue
        totals = {b: c[b][0] + c[b][1] for b in "ACGT"}
        depth = sum(totals.values())
        depth_fwd = sum(c[b][0] for b in "ACGT")
        depth_rev = sum(c[b][1] for b in "ACGT")
        nonref = depth - totals[ref]
        row = {"sample": name, "gene": genes[chrom], "chrom": chrom, "pos": pos1,
               "ref": ref, "depth": depth, "depth_fwd": depth_fwd, "depth_rev": depth_rev,
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
                             "depth_fwd": r["depth_fwd"], "depth_rev": r["depth_rev"],
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
            "mean_depth_fwd": g["depth_fwd"].mean(),
            "mean_depth_rev": g["depth_rev"].mean(),
            "strand_frac_fwd": fwd / (fwd + rev) if (fwd + rev) else np.nan,
        })
    return L.groupby(["chrom", "pos", "gene", "ref", "alt"]).apply(
        summarise, include_groups=False).reset_index()


def merge_regimes(agg_raw, agg_filt):
    """Merge the raw and filt aggregates into one cohort table (columns suffixed _raw/_filt)."""
    keys = ["chrom", "pos", "gene", "ref", "alt"]
    m = agg_raw.merge(agg_filt, on=keys, how="outer", suffixes=("_raw", "_filt"))
    numeric = [c for c in m.columns if c not in keys]
    m[numeric] = m[numeric].fillna(0)
    return m.sort_values(["chrom", "pos", "alt"])
