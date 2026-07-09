#!/usr/bin/env python3
"""Heatmap of the per-run artefact VAF (blank): rows = flagged sites, columns = runs.

Reads results/artefact_by_run.tsv (from artefact_by_run.py) and draws a heatmap of the mean
non-reference VAF (%) per site per run on a log colour scale, so run-level batch effects (a
whole column lighting up, e.g. an oxidation-damaged flow-cell) stand out at a glance. Sites and
runs are ordered worst-first. Cells at or above 1% are annotated.

Output: results/plots/artefact_by_run.png
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
RESULTS = REPO / "results"


def main():
    tsv = RESULTS / "artefact_by_run.tsv"
    if not tsv.exists():
        sys.exit(f"{tsv} not found -- run scripts/artefact_by_run.py first.")
    df = pd.read_csv(tsv, sep="\t")
    piv = df.pivot_table(index="site", columns="run", values="mean_vaf", aggfunc="mean") * 100
    # order worst-first (sites by their peak run, runs by their peak site)
    piv = piv.loc[piv.max(axis=1).sort_values(ascending=False).index]
    piv = piv[piv.max(axis=0).sort_values(ascending=False).index]

    vmax = float(np.nanmax(piv.values))
    data = np.where(piv.values > 0, piv.values, np.nan)
    nrows, ncols = piv.shape

    fig, ax = plt.subplots(figsize=(max(8, ncols * 0.42), max(5, nrows * 0.42)))
    cmap = plt.get_cmap("YlOrRd").copy()
    cmap.set_bad("#f2f2f2")
    im = ax.imshow(data, aspect="auto", cmap=cmap, norm=LogNorm(vmin=0.01, vmax=vmax))

    ax.set_xticks(range(ncols))
    ax.set_xticklabels(piv.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(nrows))
    ax.set_yticklabels(piv.index, fontsize=7)
    ax.set_xticks(np.arange(-0.5, ncols), minor=True)
    ax.set_yticks(np.arange(-0.5, nrows), minor=True)
    ax.grid(which="minor", color="white", lw=0.5)
    ax.tick_params(which="minor", length=0)

    for i in range(nrows):
        for j in range(ncols):
            v = piv.values[i, j]
            if np.isfinite(v) and v >= 1.0:
                ax.text(j, i, f"{v:.1f}", ha="center", va="center", fontsize=6,
                        color="black" if v < vmax * 0.5 else "white")

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    ticks = [t for t in (0.01, 0.1, 1, 10) if t <= vmax]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels([f"{t:g}%" for t in ticks])
    cbar.set_label("non-reference VAF (blank, log)")
    ax.set_title("Per-run artefact VAF — a bright column is a run/batch effect")
    fig.tight_layout()
    (RESULTS / "plots").mkdir(parents=True, exist_ok=True)
    out = RESULTS / "plots" / "artefact_by_run.png"
    fig.savefig(out, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out.relative_to(REPO)}  ({nrows} sites x {ncols} runs)")


if __name__ == "__main__":
    main()
