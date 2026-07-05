#!/usr/bin/env python3
"""Per-amplicon visualisation of the cohort BLANK (patient variants already masked upstream).

cohort_alt.tsv is built with masking, so these plots show only background / artefact signal --
no genuine patient mutations. One figure per amplicon, stacked top to bottom (shared x-axis):

  1. plot  -- passing QC (minbq20/minmq30): per alternate base, the mean non-reference VAF (%)
              on each strand -- forward = filled circle, reverse = triangle.
  2. table -- per position, grouped by strand: forward total depth then, per base, the reads
              supporting it and the % they represent (VAF); reverse the same, below.

Off-scale values (> YMAX%) are flagged red in the table. Genomic x-axis (GRCh38).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from amplicon_lob import panel  # noqa: E402

RESULTS = REPO / "results"
PLOTS = RESULTS / "plots"
BASE_COLOR = {"A": "#2ca02c", "C": "#1f77b4", "G": "#ff7f0e", "T": "#d62728"}  # IGV-like
BASE_OFFSET = {"A": -0.24, "C": -0.08, "G": 0.08, "T": 0.24}
YMAX = 1.0
RED = "#d62728"


def add_strand_vafs(sub):
    for tag in ("raw", "filt"):
        for s in ("fwd", "rev"):
            d = sub[f"mean_depth_{s}_{tag}"].replace(0, np.nan)
            sub[f"vaf_{s}_{tag}"] = (sub[f"mean_alt_{s}_{tag}"] / d).fillna(0.0)
    return sub


def draw_plot(ax, sub, xi, tag, panel_label, legend=False):
    for base in "ACGT":
        s = sub[sub["alt"] == base]
        xs = [xi[p] + BASE_OFFSET[base] for p in s["pos"]]
        ax.scatter(xs, s[f"vaf_fwd_{tag}"] * 100, marker="o", s=15,
                   color=BASE_COLOR[base], alpha=0.85, edgecolors="none")
        ax.scatter(xs, s[f"vaf_rev_{tag}"] * 100, marker="^", s=22,
                   color=BASE_COLOR[base], alpha=0.85, edgecolors="none")
    ax.set_ylim(0, YMAX)
    ax.set_yticks([round(0.25 * k, 2) for k in range(int(YMAX / 0.25) + 1)])
    ax.set_ylabel("% non-reference")
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.5)
    ax.text(0.004, 0.93, panel_label, transform=ax.transAxes, fontsize=7, va="top",
            fontweight="bold")
    if legend:
        handles = [Line2D([0], [0], marker="o", ls="", color=c, label=f"alt {b}")
                   for b, c in BASE_COLOR.items()]
        handles += [Line2D([0], [0], marker="o", ls="", color="grey", label="forward"),
                    Line2D([0], [0], marker="^", ls="", color="grey", label="reverse")]
        ax.legend(handles=handles, fontsize=6, ncol=6, loc="upper right", framealpha=0.9)
    plt.setp(ax.get_xticklabels(), visible=False)


def build_by_pos(sub):
    by_pos = {}
    for p, g in sub.groupby("pos"):
        rec = {"ref": g["ref"].iloc[0],
               "depth_fwd": g["mean_depth_fwd_filt"].iloc[0],
               "depth_rev": g["mean_depth_rev_filt"].iloc[0], "alts": {}}
        for _, r in g.iterrows():
            rec["alts"][r["alt"]] = {"fwd": r["mean_alt_fwd_filt"], "rev": r["mean_alt_rev_filt"],
                                     "vaf_fwd": r["vaf_fwd_filt"], "vaf_rev": r["vaf_rev_filt"]}
        by_pos[p] = rec
    return by_pos


def draw_table(tb, positions, xi, by_pos, fs):
    labels = ["fwd depth", "fwd A", "fwd C", "fwd G", "fwd T",
              "rev depth", "rev A", "rev C", "rev G", "rev T"]
    ys = list(np.linspace(0.94, 0.06, len(labels)))
    tb.set_ylim(0, 1)
    tb.set_yticks(ys)
    tb.set_yticklabels(labels, fontsize=6.5)
    for lbl in tb.get_yticklabels():
        b = lbl.get_text()[-1]
        if b in BASE_COLOR:
            lbl.set_color(BASE_COLOR[b])
    tb.set_xlim(-0.6, len(positions) - 0.4)
    yof = dict(zip(labels, ys))
    for p in positions:
        i = xi[p]
        rec = by_pos[p]
        tb.text(i, yof["fwd depth"], f"{rec['depth_fwd']:.0f}", ha="center", va="center",
                fontsize=fs, rotation=0)
        tb.text(i, yof["rev depth"], f"{rec['depth_rev']:.0f}", ha="center", va="center",
                fontsize=fs, rotation=0)
        for base in "ACGT":
            for strand, dep in (("fwd", rec["depth_fwd"]), ("rev", rec["depth_rev"])):
                key = f"{strand} {base}"
                a = rec["alts"].get(base, {})
                reads = a.get(strand, 0.0)
                vaf = a.get(f"vaf_{strand}", 0.0) * 100
                if base == rec["ref"]:
                    txt, color = "·", "#bbbbbb"
                elif reads < 0.5:
                    txt, color = "·", "#dddddd"       # no supporting reads on this strand
                else:
                    txt = f"{reads:.0f}/{vaf:.2f}"
                    color = RED if vaf > YMAX else "black"
                tb.text(i, yof[key], txt, ha="center", va="center", fontsize=fs,
                        rotation=0, color=color)
    for i in range(len(positions) + 1):
        tb.axvline(i - 0.5, color="#eeeeee", lw=0.4)
    tb.tick_params(axis="y", length=0)


def plot_amplicon(amp, agg):
    lo, hi = amp["thick_start"], amp["thick_end"]
    sub = agg[(agg["chrom"] == amp["chrom"]) & (agg["pos"] > lo) & (agg["pos"] <= hi)].copy()
    if sub.empty:
        return None
    sub = add_strand_vafs(sub)
    positions = sorted(sub["pos"].unique())
    xi = {p: i for i, p in enumerate(positions)}
    by_pos = build_by_pos(sub)
    n = len(positions)
    fs = 5.5 if n <= 60 else 4.5 if n <= 130 else 3.5

    fig = plt.figure(figsize=(min(40, max(7, n * 0.34)), 6.4))
    gs = GridSpec(2, 1, height_ratios=[2.8, 3.0], hspace=0.08)
    ax_top = fig.add_subplot(gs[0])
    ax_tbl = fig.add_subplot(gs[1], sharex=ax_top)

    ax_top.set_title(f"{amp['name']} — {n} bases length amplicon  "
                     f"(blank: patient variants removed)")
    draw_plot(ax_top, sub, xi, "filt", "passing QC  (minbq20 / minmq30)", legend=True)
    draw_table(ax_tbl, positions, xi, by_pos, fs)
    ax_tbl.set_xlim(-0.6, n - 0.4)
    ax_tbl.set_xticks(range(n))
    ax_tbl.set_xticklabels([f"{p:,}" for p in positions], rotation=90, fontsize=fs)
    ax_tbl.set_xlabel(f"{amp['chrom']} genomic coordinate (GRCh38)  "
                      f"[table cells: reads / VAF%]")

    (PLOTS / "amplicons").mkdir(parents=True, exist_ok=True)
    out = PLOTS / "amplicons" / f"{amp['name']}.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_strand_diagnostic(agg):
    fig, ax = plt.subplots(figsize=(7, 5))
    s = agg[agg["mean_vaf_raw"] > 0]
    ax.scatter(s["mean_vaf_raw"] * 100, s["strand_frac_fwd_raw"], s=10, color="#888888",
               alpha=0.5, edgecolors="none")
    ax.axhspan(0.35, 0.65, color="#1f77b4", alpha=0.06)
    ax.axhline(0.5, color="grey", ls=":", lw=1)
    ax.set_xscale("log")
    ax.set_xlim(1e-2, 10)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("% non-reference (raw blank, cohort mean, log)")
    ax.set_ylabel("forward-strand fraction of alt reads")
    ax.set_title("Strand bias vs VAF (blank): balanced ~0.5, artefacts at the extremes")
    fig.tight_layout()
    out = PLOTS / "strand_vs_vaf.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def main():
    agg = pd.read_csv(RESULTS / "cohort_alt.tsv", sep="\t")
    amps = panel.load_amplicons()
    PLOTS.mkdir(parents=True, exist_ok=True)
    n = 0
    for amp in amps:
        if plot_amplicon(amp, agg) is not None:
            n += 1
    plot_strand_diagnostic(agg)
    print(f"Wrote {n} per-amplicon plots to results/plots/amplicons/ + strand_vs_vaf.png")


if __name__ == "__main__":
    main()
