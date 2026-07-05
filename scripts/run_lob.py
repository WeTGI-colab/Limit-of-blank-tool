#!/usr/bin/env python3
"""Positional Limit-of-Blank (LoB) model for amplicon background artefacts.

Method
------
1. Masking. For every sample, non-reference observations that coincide with a call in that
   sample's VCF (genuine variants / known germline) are removed. What remains is the *blank*:
   background signal that is not biology.

2. Background model. For each (position, substitution) the blank is described across the cohort
   by the per-sample non-reference fractions. A substitution-class error floor (the random
   instrument error) is taken as the median blank rate of all positions carrying that
   substitution, robust to the few systematically noisy sites.

3. Systematic-site detection. A position/substitution is flagged as a systematic artefact when
   its blank rate rises materially above the class floor (fold-change threshold) and is
   statistically distinguishable from it (one-sided t-test on the per-sample fractions,
   Bonferroni-corrected). Strand asymmetry is reported alongside as corroborating evidence.

4. Positional Limit of Blank (CLSI EP17). LoB is estimated per position/substitution as
   mean_blank + 1.645 x SD_blank, i.e. the 95th percentile of the blank distribution. A
   candidate variant in a new sample must exceed the LoB at its own locus rather than a flat
   panel-wide VAF threshold.

Output: results/lob_table.tsv, plus a validation summary against the synthetic ground truth.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from amplicon_lob import samples  # noqa: E402

SYN = REPO / "data" / "synthetic"
MANIFEST = REPO / "data" / "manifest.tsv"
RESULTS = REPO / "results"

FOLD_THRESHOLD = 3.0      # blank rate must exceed the class floor by at least this factor
MIN_FLAG_VAF = 1e-3       # and reach at least this absolute mean blank fraction
MIN_DEPTH = 500           # positions below this mean depth are too low-coverage to trust
Z_95 = 1.645              # one-sided 95th percentile (EP17 LoB)


def load_sample_masks(manifest):
    """{sample: set((chrom, pos1, alt))} -- every VCF record is masked out of the blank.

    All records are used (not only FILTER=PASS); multi-allelic ALTs are split; chromosome
    names are normalised to the panel convention so they match the pile-up tables.
    """
    masks = {}
    for row in manifest:
        keys = set()
        vcf = row["vcf"]
        if vcf and Path(vcf).exists():
            for line in Path(vcf).read_text().splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                f = line.split("\t")
                chrom, pos = samples.normalize_chrom(f[0]), int(f[1])
                for alt in f[4].split(","):
                    keys.add((chrom, pos, alt))
        masks[row["sample"]] = keys
    return masks


def load_blank_observations(masks, manifest):
    """Long blank table: one row per (sample, position, alt) after masking genuine calls."""
    rows = []
    for row_m in manifest:
        sample = row_m["sample"]
        tsv = RESULTS / "pileup" / f"{sample}.tsv"
        if not tsv.exists():
            continue
        df = pd.read_csv(tsv, sep="\t")
        mask = masks.get(sample, set())
        for _, r in df.iterrows():
            for alt in "ACGT":
                if alt == r["ref"] or r["depth"] == 0:
                    continue
                if (r["chrom"], int(r["pos"]), alt) in mask:
                    continue                              # genuine call -> not blank
                acount = r[f"{alt}_fwd"] + r[f"{alt}_rev"]
                rows.append({
                    "sample": sample, "chrom": r["chrom"], "pos": int(r["pos"]),
                    "ref": r["ref"], "alt": alt, "sub": f"{r['ref']}>{alt}",
                    "depth": int(r["depth"]), "alt_count": int(acount),
                    "alt_fwd": int(r[f"{alt}_fwd"]), "alt_rev": int(r[f"{alt}_rev"]),
                    "vaf": acount / r["depth"],
                })
    return pd.DataFrame(rows)


def load_truth():
    """Ground truth exists only for the synthetic cohort; returns {} for production data."""
    gt = SYN / "ground_truth.tsv"
    truth = {}
    if gt.exists():
        for line in gt.read_text().splitlines()[1:]:
            f = line.split("\t")
            truth[(f[2], int(f[3]), f[5])] = f[0]
    return truth


def summarise(blank):
    """Per (chrom, pos, ref, alt, sub): blank statistics, strand fraction and LoB."""
    def agg(g):
        vafs = g["vaf"].to_numpy()
        fwd, rev = g["alt_fwd"].sum(), g["alt_rev"].sum()
        return pd.Series({
            "n_blank": len(g),
            "pooled_k": int(g["alt_count"].sum()),
            "pooled_n": int(g["depth"].sum()),
            "mean_depth": g["depth"].mean(),
            "blank_mean_vaf": vafs.mean(),
            "blank_sd_vaf": vafs.std(ddof=1) if len(vafs) > 1 else 0.0,
            "strand_frac_fwd": fwd / (fwd + rev) if (fwd + rev) else np.nan,
        })
    s = blank.groupby(["chrom", "pos", "ref", "alt", "sub"]).apply(
        agg, include_groups=False).reset_index()
    s["pooled_rate"] = s["pooled_k"] / s["pooled_n"]
    s["lob_vaf"] = s["blank_mean_vaf"] + Z_95 * s["blank_sd_vaf"]
    return s


def detect_systematic(summary, blank):
    """Flag systematic artefacts relative to a substitution-class error floor."""
    floor = summary.groupby("sub")["pooled_rate"].median().to_dict()
    summary["class_floor"] = summary["sub"].map(floor)
    summary["fold_over_floor"] = summary["pooled_rate"] / summary["class_floor"]

    per_site = {k: g["vaf"].to_numpy() for k, g in
                blank.groupby(["chrom", "pos", "alt"])}
    n_tests = len(summary)
    alpha = 0.05 / max(1, n_tests)

    pvals, flags = [], []
    for _, r in summary.iterrows():
        vafs = per_site.get((r["chrom"], r["pos"], r["alt"]), np.array([]))
        f0 = r["class_floor"]
        if len(vafs) > 1 and vafs.std(ddof=1) > 0:
            t, p_two = stats.ttest_1samp(vafs, f0)
            p = p_two / 2 if t > 0 else 1.0            # one-sided (greater)
        else:
            p = 0.0 if r["pooled_rate"] > f0 else 1.0
        flagged = bool(r["fold_over_floor"] >= FOLD_THRESHOLD
                       and r["blank_mean_vaf"] >= MIN_FLAG_VAF
                       and r["mean_depth"] >= MIN_DEPTH and p < alpha)
        pvals.append(p)
        flags.append(flagged)
    summary["pval"] = pvals
    summary["flagged"] = flags
    return summary, alpha


def validate(summary, truth):
    summary["truth"] = [truth.get((c, p, a), "") for c, p, a in
                        zip(summary["chrom"], summary["pos"], summary["alt"])]
    art = summary[summary["truth"] == "artefact"]
    real = summary[summary["truth"] == "real"]
    flagged = summary[summary["flagged"]]

    tp = int((art["flagged"]).sum())
    fn = int((~art["flagged"]).sum())
    fp = int(((summary["truth"] != "artefact") & summary["flagged"]).sum())

    print("\n=== Validation against synthetic ground truth ===")
    print(f"Ground-truth artefacts detected : {tp}/{len(art)}  (false negatives={fn})")
    print(f"False positives (flagged, not a true artefact): {fp}")
    print(f"Genuine variants flagged as artefact (should be 0, they are VCF-masked): "
          f"{int(real['flagged'].sum())}/{len(real)}")

    cols = ["gene_pos", "sub", "blank_mean_vaf", "fold_over_floor",
            "strand_frac_fwd", "lob_vaf", "flagged", "truth"]
    show = summary.copy()
    show["gene_pos"] = show["chrom"] + ":" + show["pos"].astype(str)
    view = show[(show["flagged"]) | (show["truth"] != "")].sort_values(
        "blank_mean_vaf", ascending=False)
    print("\nFlagged sites and ground-truth loci:")
    with pd.option_context("display.width", 200, "display.max_columns", None):
        print(view[cols].to_string(index=False, formatters={
            "blank_mean_vaf": "{:.5f}".format, "fold_over_floor": "{:.1f}".format,
            "strand_frac_fwd": "{:.2f}".format, "lob_vaf": "{:.5f}".format}))


def main():
    manifest = samples.load_manifest(MANIFEST)
    masks = load_sample_masks(manifest)
    blank = load_blank_observations(masks, manifest)
    truth = load_truth()
    summary = summarise(blank)
    summary, alpha = detect_systematic(summary, blank)

    summary["truth"] = [truth.get((c, p, a), "") for c, p, a in
                        zip(summary["chrom"], summary["pos"], summary["alt"])]
    out = summary.sort_values(["chrom", "pos", "alt"])
    out.to_csv(RESULTS / "lob_table.tsv", sep="\t", index=False)
    n_flagged = int(out["flagged"].sum())
    print(f"Wrote lob_table.tsv  ({len(out)} position x alt rows, "
          f"{n_flagged} flagged systematic; Bonferroni alpha={alpha:.2e})")
    if truth:
        validate(out, truth)
    else:
        show = out[out["flagged"]].copy()
        show["locus"] = show["chrom"] + ":" + show["pos"].astype(str)
        print(f"\n{len(show)} systematic artefact site(s) flagged:")
        cols = ["locus", "sub", "n_blank", "blank_mean_vaf", "fold_over_floor",
                "strand_frac_fwd", "lob_vaf"]
        with pd.option_context("display.width", 200, "display.max_columns", None):
            print(show.sort_values("blank_mean_vaf", ascending=False)[cols].to_string(
                index=False, formatters={"blank_mean_vaf": "{:.5f}".format,
                                         "fold_over_floor": "{:.1f}".format,
                                         "strand_frac_fwd": "{:.2f}".format,
                                         "lob_vaf": "{:.5f}".format}))


if __name__ == "__main__":
    main()
