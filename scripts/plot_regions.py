#!/usr/bin/env python3
"""Per-amplicon visualisation of the cohort BLANK (patient variants already masked upstream).

One figure per amplicon, stacked top to bottom (shared genomic x-axis):

  1. table -- quality 0  (minbq1 / minmq1, all reads): per position, by strand -- forward
              depth then per-base reads/VAF, reverse the same below.
  2. plot  -- quality 0: per alternate base, mean non-reference VAF (%) per strand
              (forward = filled circle, reverse = triangle).
  3. plot  -- quality 20 (minbq20 / minmq30, caller-callable): same.
  4. table -- quality 20: same layout as the top table.

A strand's VAF is only plotted where that strand has >= MIN_STRAND_DEPTH coverage (low-coverage
strands give spurious VAFs). Off-scale values (> YMAX%) are flagged red in the tables. The
"strand" is the amplicon direction (_F_/_R_ from the read CO tag), not the BAM reverse flag.
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
MIN_STRAND_DEPTH = 500   # only plot a strand's VAF where that strand is deep enough to trust
TAG_LABEL = {"raw": "quality 0  (minbq1 / minmq1)", "filt": "quality 20  (minbq20 / minmq30)"}


def add_strand_vafs(sub):
    for tag in ("raw", "filt"):
        for s in ("fwd", "rev"):
            d = sub[f"mean_depth_{s}_{tag}"].replace(0, np.nan)
            sub[f"vaf_{s}_{tag}"] = (sub[f"mean_alt_{s}_{tag}"] / d).fillna(0.0)
    return sub


def draw_plot(ax, sub, xi, tag, legend=False):
    for base in "ACGT":
        s = sub[sub["alt"] == base]
        sf = s[s[f"mean_depth_fwd_{tag}"] >= MIN_STRAND_DEPTH]
        ax.scatter([xi[p] + BASE_OFFSET[base] for p in sf["pos"]], sf[f"vaf_fwd_{tag}"] * 100,
                   marker="o", s=15, color=BASE_COLOR[base], alpha=0.85, edgecolors="none")
        sr = s[s[f"mean_depth_rev_{tag}"] >= MIN_STRAND_DEPTH]
        ax.scatter([xi[p] + BASE_OFFSET[base] for p in sr["pos"]], sr[f"vaf_rev_{tag}"] * 100,
                   marker="^", s=22, color=BASE_COLOR[base], alpha=0.85, edgecolors="none")
    ax.set_ylim(0, YMAX)
    ax.set_yticks([round(0.25 * k, 2) for k in range(int(YMAX / 0.25) + 1)])
    ax.set_ylabel("% non-ref")
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.5)
    ax.text(0.004, 0.9, TAG_LABEL[tag], transform=ax.transAxes, fontsize=7, va="top",
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
        rec = {"ref": g["ref"].iloc[0]}
        for tag in ("raw", "filt"):
            alts = {}
            for _, r in g.iterrows():
                alts[r["alt"]] = {"fwd": r[f"mean_alt_fwd_{tag}"], "rev": r[f"mean_alt_rev_{tag}"],
                                  "vaf_fwd": r[f"vaf_fwd_{tag}"], "vaf_rev": r[f"vaf_rev_{tag}"]}
            rec[tag] = {"depth_fwd": g[f"mean_depth_fwd_{tag}"].iloc[0],
                        "depth_rev": g[f"mean_depth_rev_{tag}"].iloc[0], "alts": alts}
        by_pos[p] = rec
    return by_pos


def draw_table(tb, positions, xi, by_pos, tag, fs, label):
    labels = ["fwd depth", "fwd A", "fwd C", "fwd G", "fwd T",
              "rev depth", "rev A", "rev C", "rev G", "rev T"]
    ys = list(np.linspace(0.93, 0.07, len(labels)))
    yof = dict(zip(labels, ys))
    tb.set_ylim(0, 1)
    tb.set_yticks(ys)
    tb.set_yticklabels(labels, fontsize=6.5)
    for lbl in tb.get_yticklabels():
        b = lbl.get_text()[-1]
        if b in BASE_COLOR:
            lbl.set_color(BASE_COLOR[b])
    tb.set_xlim(-0.6, len(positions) - 0.4)
    tb.text(0.004, 0.99, label, transform=tb.transAxes, fontsize=7, va="top", fontweight="bold")
    for p in positions:
        i = xi[p]
        rec = by_pos[p]
        d = rec[tag]
        tb.text(i, yof["fwd depth"], f"{d['depth_fwd']:.0f}", ha="center", va="center",
                fontsize=fs, rotation=0)
        tb.text(i, yof["rev depth"], f"{d['depth_rev']:.0f}", ha="center", va="center",
                fontsize=fs, rotation=0)
        for base in "ACGT":
            for strand in ("fwd", "rev"):
                a = d["alts"].get(base, {})
                reads = a.get(strand, 0.0)
                vaf = a.get(f"vaf_{strand}", 0.0) * 100
                if base == rec["ref"]:
                    txt, color = "·", "#bbbbbb"
                elif reads < 0.5:
                    txt, color = "·", "#dddddd"
                else:
                    txt = f"{reads:.0f}/{vaf:.2f}"
                    color = RED if vaf > YMAX else "black"
                tb.text(i, yof[f"{strand} {base}"], txt, ha="center", va="center",
                        fontsize=fs, rotation=0, color=color)
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

    fig = plt.figure(figsize=(min(40, max(7, n * 0.34)), 10.5))
    gs = GridSpec(4, 1, height_ratios=[2.2, 2.3, 2.3, 2.2], hspace=0.10)
    ax_t0 = fig.add_subplot(gs[0])
    ax_p0 = fig.add_subplot(gs[1], sharex=ax_t0)
    ax_p1 = fig.add_subplot(gs[2], sharex=ax_t0)
    ax_t1 = fig.add_subplot(gs[3], sharex=ax_t0)

    ax_t0.set_title(f"{amp['name']} — {n} bases length amplicon  (blank: patient variants removed)")
    draw_table(ax_t0, positions, xi, by_pos, "raw", fs, "quality 0  (minbq1)")
    plt.setp(ax_t0.get_xticklabels(), visible=False)
    draw_plot(ax_p0, sub, xi, "raw", legend=True)
    draw_plot(ax_p1, sub, xi, "filt")
    draw_table(ax_t1, positions, xi, by_pos, "filt", fs, "quality 20  (minbq20)")

    ax_t1.set_xlim(-0.6, n - 0.4)
    ax_t1.set_xticks(range(n))
    ax_t1.set_xticklabels([f"{p:,}" for p in positions], rotation=90, fontsize=fs)
    ax_t1.set_xlabel(f"{amp['chrom']} genomic coordinate (GRCh38)  [table cells: reads / VAF%]")

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
    ax.set_ylabel("forward-amplicon fraction of alt reads")
    ax.set_title("Strand (amplicon F/R) bias vs VAF (blank)")
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
