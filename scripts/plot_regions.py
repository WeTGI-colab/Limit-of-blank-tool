#!/usr/bin/env python3
"""Per-amplicon visualisation of the cohort pile-up under two QC regimes.

One figure per amplicon, stacked top to bottom (shared genomic x-axis):

  1. table  -- % non-reference per alternate base counting ONLY reads passing the laboratory
               filters (minbq=20, minmq=30): the noise that survives QC.
  2. plot   -- points of that filtered % (the table above), y-axis capped at YMAX%.
  3. plot   -- points of the raw % (minbq=1, minmq=1): all sequencer noise.
  4. table  -- raw read depth and the reads supporting each alternate base (all reads).

In both tables the reference base is shown as '.'; a value is coloured red when its VAF
exceeds the y-axis maximum (off-scale above the plot).

Also writes results/plots/strand_vs_vaf.png (global strand-bias diagnostic, raw counts).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.gridspec import GridSpec  # noqa: E402
import pandas as pd  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
from amplicon_lob import panel  # noqa: E402

RESULTS = REPO / "results"
PLOTS = RESULTS / "plots"
BASE_COLOR = {"A": "#2ca02c", "C": "#1f77b4", "G": "#ff7f0e", "T": "#d62728"}  # IGV-like
BASE_OFFSET = {"A": -0.24, "C": -0.08, "G": 0.08, "T": 0.24}
YMAX = 1.0            # y-axis maximum (%) for both plots
RED = "#d62728"


def draw_plot(ax, sub, xi, valcol, panel_label, legend=False):
    for base in "ACGT":
        s = sub[sub["alt"] == base]
        ax.scatter([xi[p] + BASE_OFFSET[base] for p in s["pos"]], s[valcol] * 100,
                   s=14, color=BASE_COLOR[base], alpha=0.85, edgecolors="none",
                   label=f"alt {base}")
    for _, r in sub[sub["truth"].isin(["artefact", "real"])].iterrows():
        v = r[valcol] * 100
        if v > YMAX:
            continue
        mark = RED if r["truth"] == "artefact" else "#1f77b4"
        ax.scatter([xi[r["pos"]] + BASE_OFFSET[r["alt"]]], [v], s=40, facecolors="none",
                   edgecolors=mark, linewidths=1.2)
        ax.annotate(r["truth"], (xi[r["pos"]] + BASE_OFFSET[r["alt"]], v),
                    textcoords="offset points", xytext=(0, 4), ha="center",
                    fontsize=6, color=mark, rotation=90)
    ax.set_ylim(0, YMAX)
    ax.set_yticks([0.5 * k for k in range(int(YMAX / 0.5) + 1)])
    ax.set_ylabel("% non-reference")
    ax.grid(axis="y", ls=":", lw=0.5, alpha=0.5)
    ax.text(0.004, 0.94, panel_label, transform=ax.transAxes, fontsize=7,
            va="top", fontweight="bold")
    if legend:
        ax.legend(fontsize=6, ncol=4, loc="upper right", framealpha=0.9)
    plt.setp(ax.get_xticklabels(), visible=False)


def draw_table(tb, positions, xi, by_pos, mode, fs):
    """mode 'filt_pct' -> depth + per-base %; mode 'raw_reads' -> depth + per-base read counts."""
    row_y = {"depth": 0.88, "A": 0.70, "C": 0.53, "G": 0.36, "T": 0.19}
    depth_key = "depth_filt" if mode == "filt_pct" else "depth_raw"
    tb.set_ylim(0, 1)
    tb.set_yticks(list(row_y.values()))
    tb.set_yticklabels(["depth", "alt A", "alt C", "alt G", "alt T"], fontsize=7)
    for lbl, key in zip(tb.get_yticklabels(), ["depth", "A", "C", "G", "T"]):
        if key in BASE_COLOR:
            lbl.set_color(BASE_COLOR[key])
    tb.set_xlim(-0.6, len(positions) - 0.4)
    for p in positions:
        i = xi[p]
        rec = by_pos[p]
        tb.text(i, row_y["depth"], f"{rec[depth_key]:.0f}", ha="center", va="center",
                fontsize=fs, rotation=0)
        for base in "ACGT":
            if base == rec["ref"]:
                txt, color = "·", "#bbbbbb"
            else:
                a = rec["alts"].get(base, {})
                over = a.get("filt_pct" if mode == "filt_pct" else "raw_pct", 0) > YMAX
                if mode == "filt_pct":
                    txt = f"{a.get('filt_pct', 0):.2f}"
                else:
                    txt = f"{a.get('raw_reads', 0):.1f}"
                color = RED if over else "black"
            tb.text(i, row_y[base], txt, ha="center", va="center", fontsize=fs,
                    rotation=0, color=color)
    for i in range(len(positions) + 1):
        tb.axvline(i - 0.5, color="#eeeeee", lw=0.4)
    tb.tick_params(axis="y", length=0)


def build_by_pos(sub):
    by_pos = {}
    for p, g in sub.groupby("pos"):
        ref = g["ref"].iloc[0]
        alts = {}
        for _, r in g.iterrows():
            alts[r["alt"]] = {
                "raw_pct": r["mean_vaf_raw"] * 100, "filt_pct": r["mean_vaf_filt"] * 100,
                "raw_reads": r["alt_reads_mean_raw"], "truth": r["truth"],
            }
        by_pos[p] = {"ref": ref, "depth_raw": g["mean_depth_raw"].iloc[0],
                     "depth_filt": g["mean_depth_filt"].iloc[0], "alts": alts}
    return by_pos


def plot_amplicon(amp, agg):
    lo, hi = amp["thick_start"], amp["thick_end"]
    sub = agg[(agg["chrom"] == amp["chrom"]) & (agg["pos"] > lo) & (agg["pos"] <= hi)].copy()
    if sub.empty:
        return None
    positions = sorted(sub["pos"].unique())
    xi = {p: i for i, p in enumerate(positions)}
    by_pos = build_by_pos(sub)
    n = len(positions)
    fs = 6 if n <= 60 else 5 if n <= 130 else 4

    fig = plt.figure(figsize=(min(34, max(6, n * 0.30)), 8.4))
    gs = GridSpec(4, 1, height_ratios=[1.5, 2.6, 2.6, 1.5], hspace=0.08)
    ax_ttbl = fig.add_subplot(gs[0])
    ax_top = fig.add_subplot(gs[1], sharex=ax_ttbl)
    ax_bot = fig.add_subplot(gs[2], sharex=ax_ttbl)
    ax_btbl = fig.add_subplot(gs[3], sharex=ax_ttbl)

    ax_ttbl.set_title(f"{amp['name']} — {n} bases length amplicon")
    draw_table(ax_ttbl, positions, xi, by_pos, "filt_pct", fs)
    plt.setp(ax_ttbl.get_xticklabels(), visible=False)
    draw_plot(ax_top, sub, xi, "mean_vaf_filt", "passing QC  (minbq20 / minmq30)", legend=True)
    draw_plot(ax_bot, sub, xi, "mean_vaf_raw", "raw  (all reads, minbq1 / minmq1)")

    draw_table(ax_btbl, positions, xi, by_pos, "raw_reads", fs)
    ax_btbl.set_xlim(-0.6, n - 0.4)
    ax_btbl.set_xticks(range(n))
    ax_btbl.set_xticklabels([f"{p:,}" for p in positions], rotation=90, fontsize=fs)
    arte = set(sub[sub["truth"] == "artefact"]["pos"])
    real = set(sub[sub["truth"] == "real"]["pos"])
    for lbl, p in zip(ax_btbl.get_xticklabels(), positions):
        if p in arte:
            lbl.set_color(RED)
        elif p in real:
            lbl.set_color("#1f77b4")
    ax_btbl.set_xlabel(f"{amp['chrom']} genomic coordinate (GRCh38)")

    (PLOTS / "amplicons").mkdir(parents=True, exist_ok=True)
    out = PLOTS / "amplicons" / f"{amp['name']}.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_strand_diagnostic(agg):
    fig, ax = plt.subplots(figsize=(7, 5))
    for key, label, color, size, alpha in [
            ("", "baseline / error floor", "#bbbbbb", 8, 0.4),
            ("real", "genuine variant", "#1f77b4", 70, 0.95),
            ("artefact", "systematic artefact", RED, 70, 0.95)]:
        s = agg[(agg["truth"] == key) & (agg["mean_vaf_raw"] > 0)]
        ax.scatter(s["mean_vaf_raw"] * 100, s["strand_frac_fwd_raw"], s=size, color=color,
                   label=label, alpha=alpha, edgecolors="none")
    ax.axhspan(0.35, 0.65, color="#1f77b4", alpha=0.06)
    ax.axhline(0.5, color="grey", ls=":", lw=1)
    ax.set_xscale("log")
    ax.set_xlim(1e-2, 10)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("% non-reference (raw, cohort mean, log)")
    ax.set_ylabel("forward-strand fraction of alt reads")
    ax.set_title("Artefact vs genuine variant: strand bias against VAF")
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout()
    out = PLOTS / "strand_vs_vaf.png"
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def main():
    agg = pd.read_csv(RESULTS / "cohort_alt.tsv", sep="\t")
    agg["truth"] = agg["truth"].fillna("").astype(str).replace("nan", "")
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
