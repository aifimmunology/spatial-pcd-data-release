#!/usr/bin/env python3
"""
Stage 05 (diagnostic) — how much does the QC cutoff change the UMAP?

Sweeps a per-sample drop budget. At each level it derives the largest per-sample
cutoffs that fit the budget, filters, and runs the SAME pipeline as
03_process.py (_common.process, not a copy), so the differences come from the
cutoff and nothing else. Seeds are fixed, so what you see is not run-to-run noise.

Metrics, to argue from rather than eyeball:
    n_spots, median_counts   what survived
    n_clusters               Leiden clusters at fixed resolution
    depth_rho                |Spearman| between UMI counts and the UMAP axes —
                             how much of the layout is just sequencing depth
    low_depth_frac           fraction of spots in clusters whose median count is
                             below LOW_DEPTH: depth-driven clusters
    batch_entropy            mean per-cluster sample mixing, 0-1; low means
                             clusters are sample-specific despite Harmony

Reads:  h5ad/qc.h5ad
Writes: figures/05_umap_sweep_*.{png,pdf}, tables/05_cutoff_sweep_*.csv

Roughly one 03_process.py run per budget — about 12 min for 11 budgets.

    python visium/scripts/spatial/05_cutoff_sweep.py
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from _common import (
    QC_H5AD,
    SAMPLE_KEY,
    cutoffs_for_budget,
    dropped_mask,
    process,
    read_h5ad,
    save_fig,
    save_table,
)

BUDGETS = list(range(0, 21, 2))  # per-sample drop budgets, %
LOW_DEPTH = 50                   # a cluster below this median UMI is depth-driven
NCOLS = 4

qc = read_h5ad(QC_H5AD)
n_samples = qc.obs[SAMPLE_KEY].nunique()

results, metrics, cutoff_rows = {}, [], []
for budget in BUDGETS:
    cutoffs = cutoffs_for_budget(qc.obs, budget / 100.0)
    drop = dropped_mask(qc.obs, cutoffs)
    sub = qc[~drop].copy()

    for sample, t in cutoffs.items():
        sel = qc.obs[SAMPLE_KEY] == sample
        cutoff_rows.append({"budget_pct": budget, "sample": sample, **t,
                            "n_spots": int(sel.sum()), "dropped": int(drop[sel].sum())})

    process(sub)

    results[budget] = pd.DataFrame({
        "umap1": sub.obsm["X_umap"][:, 0],
        "umap2": sub.obsm["X_umap"][:, 1],
        "sample": sub.obs[SAMPLE_KEY].to_numpy(),
        "leiden": sub.obs["leiden"].to_numpy(),
        "total_counts": sub.obs["total_counts"].to_numpy(),
    })

    counts = sub.obs["total_counts"].to_numpy()
    umap = sub.obsm["X_umap"]
    clusters = sub.obs["leiden"]
    median_per_cluster = sub.obs.groupby(clusters, observed=True)["total_counts"].median()
    # Sample mixing per cluster, normalised to [0, 1]; Harmony should push it up
    entropies = []
    for _, rows in sub.obs.groupby(clusters, observed=True):
        p = rows[SAMPLE_KEY].value_counts(normalize=True)
        p = p[p > 0]
        entropies.append(-(p * np.log(p)).sum() / np.log(n_samples))

    metrics.append({
        "budget_pct": budget,
        "dropped_pct": round(100 * drop.mean(), 2),
        "n_spots": int(sub.n_obs),
        "median_counts": float(np.median(counts)),
        "n_clusters": int(clusters.nunique()),
        "depth_rho": round(max(abs(spearmanr(counts, umap[:, 0]).statistic),
                               abs(spearmanr(counts, umap[:, 1]).statistic)), 3),
        "low_depth_frac": round(float(clusters.isin(
            median_per_cluster[median_per_cluster < LOW_DEPTH].index).mean()), 3),
        "batch_entropy": round(float(np.mean(entropies)), 3),
    })
    print(metrics[-1])
    del sub

metrics = pd.DataFrame(metrics).set_index("budget_pct")
print("\n=== sweep summary ===")
print(metrics.to_string())
save_table(metrics, "05_cutoff_sweep_metrics.csv")
save_table(pd.DataFrame(cutoff_rows).set_index(["budget_pct", "sample"]),
           "05_cutoff_sweep_cutoffs.csv")

# --- One UMAP grid per colouring -------------------------------------------
TITLES = {
    "counts": "coloured by sequencing depth — is the layout just depth?",
    "sample": "coloured by sample — are samples mixing, or splitting apart?",
    "leiden": "coloured by Leiden cluster",
}
nrows = int(np.ceil(len(BUDGETS) / NCOLS))

for color, subtitle in TITLES.items():
    fig, axes = plt.subplots(nrows, NCOLS, figsize=(4.2 * NCOLS, 4.0 * nrows),
                             squeeze=False)
    axes = axes.flatten()
    for ax, budget in zip(axes, BUDGETS):
        df = results[budget]
        m = metrics.loc[budget]
        values = {
            "counts": np.log10(df.total_counts + 1),
            "sample": pd.Categorical(df["sample"]).codes,
            "leiden": pd.Categorical(df["leiden"]).codes,
        }[color]
        ax.scatter(df.umap1, df.umap2, s=1.5, c=values,
                   cmap="viridis" if color == "counts" else "tab20",
                   edgecolors="none")
        ax.set_title(f"drop <={budget}%  ({m.dropped_pct:.1f}% dropped)\n"
                     f"{int(m.n_spots):,} spots - median {int(m.median_counts)} UMI - "
                     f"{int(m.n_clusters)} clusters\n"
                     f"depth rho={m.depth_rho:.2f} - low-depth {m.low_depth_frac:.0%}",
                     fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])
    for ax in axes[len(BUDGETS):]:
        ax.axis("off")

    fig.suptitle(f"UMAP vs QC cutoff — {subtitle}\n"
                 "same pipeline and seed at every level; only the cutoff changes",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save_fig(f"05_umap_sweep_{color}")
