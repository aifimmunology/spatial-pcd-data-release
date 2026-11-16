"""
fh1_pretreatment AZU1-removed run — figure pass on the Pearson-residual clustering.

Mirrors the 03 figure set (diagnostic UMAP, per-subject highlights, per-resolution
cluster UMAPs + highlight grids, sweep grid) but drawn from
combined_pearson_clustered.h5ad, so the Pearson version of this run has the same
figures as the total-count version for side-by-side comparison. Plotting only —
the embedding and clusters already exist in the object.

Outputs go to a SEPARATE folder (figures/archive/azu1_removed_pearson/); the
total-count figures are left untouched. Read-only w.r.t. data.

    python xenium/fh1_pretreatment/scripts/archive/azu1_removed/07_pearson_figures.py
"""

import math
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

DATA_DIR = "/home/workspace/PCD_data/xenium/fh1_pretreatment/data/archive"
RUN_DIR = os.path.join(DATA_DIR, "azu1_removed")
FIG_DIR = "/home/workspace/PCD_data/xenium/fh1_pretreatment/figures/archive/azu1_removed_pearson"
CLUSTERS_FIG_DIR = os.path.join(FIG_DIR, "resolution_sweep", "clusters")
HIGHLIGHTS_FIG_DIR = os.path.join(FIG_DIR, "resolution_sweep", "highlights")
for d in (FIG_DIR, CLUSTERS_FIG_DIR, HIGHLIGHTS_FIG_DIR):
    os.makedirs(d, exist_ok=True)

RESOLUTIONS = [round(0.1 * i, 2) for i in range(1, 12)]
SWEEP_KEYS = [f"leiden_res{r:.2f}" for r in RESOLUTIONS]
HIGHLIGHT_NCOLS = 4


def save(fig, name):
    """Rasterize scatter points (light PDF), then save png + pdf. The post-hoc
    set_rasterized on each collection is the version-robust approach."""
    for ax in fig.axes:
        for coll in ax.collections:
            coll.set_rasterized(True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def save_to(fig, subdir, name):
    for ax in fig.axes:
        for coll in ax.collections:
            coll.set_rasterized(True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(subdir, f"{name}.{ext}"), dpi=150, bbox_inches="tight")
    plt.close(fig)


combined = sc.read_h5ad(os.path.join(RUN_DIR, "combined_pearson_clustered.h5ad"))
print(f"Loaded Pearson object: {combined.n_obs:,} cells x {combined.n_vars} genes")
combined.obs["cohort_subject"] = combined.obs["cohort_subject"].astype("category")

# ── Diagnostic UMAP ─────────────────────────────────────────────────────────────
fig = sc.pl.umap(combined,
                 color=["total_counts", "n_genes_by_counts", "cohort_subject", "treatment_timepoint"],
                 wspace=0.3, ncols=2, show=False, return_fig=True, size=6)
fig.suptitle("Pearson — Harmony diagnostic")
save(fig, "umap_harmony_diagnostic")
print("saved umap_harmony_diagnostic")

# ── Per-subject highlight grid ──────────────────────────────────────────────────
subjects = list(combined.obs["cohort_subject"].cat.categories)
nrows = math.ceil((len(subjects) + 1) / HIGHLIGHT_NCOLS)
fig, axes = plt.subplots(nrows, HIGHLIGHT_NCOLS, figsize=(HIGHLIGHT_NCOLS * 3, nrows * 3))
axes = axes.ravel()
sc.pl.umap(combined, color="cohort_subject", ax=axes[0], size=6,
           legend_loc="right margin", legend_fontsize=6, show=False)
axes[0].set_title("all subjects")
for ax, subject in zip(axes[1:], subjects):
    sc.pl.umap(combined, color="cohort_subject", groups=[subject], ax=ax, size=6,
               legend_loc=None, na_in_legend=False, show=False)
    ax.set_title(subject)
for ax in axes[len(subjects) + 1:]:
    ax.axis("off")
fig.suptitle("Pearson — cohort_subject")
fig.tight_layout()
save(fig, "umap_subject_highlights")
print("saved umap_subject_highlights")

# ── Per-resolution cluster UMAPs ────────────────────────────────────────────────
for key in SWEEP_KEYS:
    fig = sc.pl.umap(combined, color=key, legend_loc="on data", legend_fontsize=12,
                     size=6, show=False, return_fig=True)
    fig.axes[0].set_title(f"Pearson — {key}")
    save_to(fig, CLUSTERS_FIG_DIR, f"umap_{key}")
print(f"saved {len(SWEEP_KEYS)} per-resolution UMAPs")

# ── Per-resolution cluster-highlight grids ──────────────────────────────────────
for key in SWEEP_KEYS:
    clusters = list(combined.obs[key].cat.categories)
    nrows = math.ceil((len(clusters) + 1) / HIGHLIGHT_NCOLS)
    fig, axes = plt.subplots(nrows, HIGHLIGHT_NCOLS, figsize=(HIGHLIGHT_NCOLS * 3, nrows * 3))
    axes = axes.ravel()
    sc.pl.umap(combined, color=key, ax=axes[0], size=6, legend_loc="on data",
               legend_fontsize=6, show=False)
    axes[0].set_title("all clusters")
    for ax, cluster in zip(axes[1:], clusters):
        sc.pl.umap(combined, color=key, groups=[cluster], ax=ax, size=6,
                   legend_loc=None, na_in_legend=False, show=False)
        ax.set_title(f"cluster {cluster}")
    for ax in axes[len(clusters) + 1:]:
        ax.axis("off")
    fig.suptitle(f"Pearson — {key}")
    fig.tight_layout()
    save_to(fig, HIGHLIGHTS_FIG_DIR, f"umap_{key}_highlights")
print(f"saved {len(SWEEP_KEYS)} cluster-highlight grids")

# ── Resolution sweep grid ───────────────────────────────────────────────────────
fig = sc.pl.umap(combined, color=SWEEP_KEYS, wspace=0.05, ncols=3,
                 legend_loc="on data", show=False, return_fig=True, size=6)
fig.suptitle("Pearson — resolution sweep")
save(fig, "umap_resolution_sweep")
print("saved umap_resolution_sweep")

print(f"\nOutputs → {FIG_DIR}")
