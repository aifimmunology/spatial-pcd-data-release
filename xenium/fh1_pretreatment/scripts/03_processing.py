import math
import os
import warnings

import matplotlib.pyplot as plt
import scanpy as sc
import harmonypy as hm

warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

# Rasterize the scatter points in vector output (points become an embedded image,
# text/axes stay vector): shrinks each UMAP PDF from ~7 MB to ~150 KB. scanpy
# hardcodes rasterized=settings._vector_friendly in its scatter call (default
# False), so this flag is the switch — passing rasterized= to sc.pl.umap directly
# raises a TypeError. Setting the flag rather than set_figure_params() avoids that
# function's font-size/dpi side effects.
sc.settings._vector_friendly = True

DATA_DIR      = '/home/workspace/PCD_data/xenium/fh1_pretreatment/data'
OBJECTS_DIR   = os.path.join(DATA_DIR, 'objects')
SWEEP_DIR     = os.path.join(DATA_DIR, 'resolution_sweep')
FIG_DIR       = '/home/workspace/PCD_data/xenium/fh1_pretreatment/figures'
CLUSTERS_DIR  = os.path.join(FIG_DIR, 'clusters')
HIGHLIGHTS_DIR = os.path.join(FIG_DIR, 'highlights')
os.makedirs(OBJECTS_DIR, exist_ok=True)
os.makedirs(SWEEP_DIR, exist_ok=True)
os.makedirs(CLUSTERS_DIR, exist_ok=True)
os.makedirs(HIGHLIGHTS_DIR, exist_ok=True)

BATCH_KEY       = 'cohort_subject'
RESOLUTIONS     = [round(0.1 * i, 2) for i in range(1, 12)]  # 0.1, 0.2, …, 1.1
HIGHLIGHT_NCOLS = 4

combined = sc.read_h5ad(os.path.join(OBJECTS_DIR, 'combined_filtered.h5ad'))
print(f'Loaded: {combined.n_obs:,} cells x {combined.n_vars} genes')

# ── Pearson-residual embedding ────────────────────────────────────────────────
# PCA runs on analytic Pearson residuals (the scanpy-native form of SCTransform):
# each count is replaced by how far it sits from what a cell with that total depth
# would be expected to show, which replaces normalize_total + log1p + scale in one
# depth- and variance-stabilized step. Only the embedding is built from the
# residuals; combined.X stays log-normalized so marker tests and expression plots
# read normally.
# The residuals live in a throwaway object `resid`; `combined` itself is not altered here.
resid = sc.AnnData(combined.layers['counts'].copy())
sc.experimental.pp.normalize_pearson_residuals(resid)
sc.pp.pca(resid, n_comps=30, random_state=0)
# Adds obsm['X_pca']: the 30 components from the Pearson residuals, not from .X.
combined.obsm['X_pca'] = resid.obsm['X_pca']

# ── Harmony batch correction on patient ───────────────────────────────────────
# Harmony shifts each patient's cells in PCA space so that patient identity alone
# does not drive the clustering. All 13 samples are the same MM Pre-Treatment
# timepoint, so patient (cohort_subject) is the batch here.
ho = hm.run_harmony(combined.obsm['X_pca'].astype('float64'), combined.obs,
                    BATCH_KEY, random_state=0)
z = ho.Z_corr  # shape may be (n_pcs x n_cells) or (n_cells x n_pcs) depending on version
# Adds obsm['X_pca_harmony']: the batch-corrected coordinates. obsm['X_pca'] stays uncorrected.
combined.obsm['X_pca_harmony'] = z.T if z.shape[0] != combined.n_obs else z

# Adds obsp['connectivities'] and obsp['distances'] (each cell's 15 nearest neighbors, measured
# in the Harmony space) and uns['neighbors'] (the settings used).
sc.pp.neighbors(combined, use_rep='X_pca_harmony', n_neighbors=15, n_pcs=30, random_state=0)
# Adds obsm['X_umap']: a 2D position per cell computed from the neighbor graph, plus uns['umap'].
sc.tl.umap(combined, min_dist=0.5, random_state=0)

# ── Diagnostic UMAP: verify batch correction ──────────────────────────────────
# Cells should intermix by patient after Harmony; strong per-patient separation
# here would mean the correction did not converge.
diag = ['total_counts', 'n_genes_by_counts', BATCH_KEY]
if 'treatment_timepoint' in combined.obs:
    diag.append('treatment_timepoint')
fig = sc.pl.umap(combined, color=diag, ncols=2, wspace=0.3, size=6,
                 show=False, return_fig=True)
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'umap_diagnostic.{ext}'), dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved umap_diagnostic to {FIG_DIR}')

# ── Per-patient highlight grid ────────────────────────────────────────────────
# Contact sheet: first panel colors every patient; each following panel isolates
# one patient (its cells colored, all others greyed via groups=) to show where
# each patient sits in the shared embedding after Harmony.
# Changes the dtype of obs['cohort_subject'] to categorical; the values are unchanged.
# Plotting a categorical column also adds uns['cohort_subject_colors'].
combined.obs[BATCH_KEY] = combined.obs[BATCH_KEY].astype('category')
groups = list(combined.obs[BATCH_KEY].cat.categories)
nrows = math.ceil((len(groups) + 1) / HIGHLIGHT_NCOLS)
fig, axes = plt.subplots(nrows, HIGHLIGHT_NCOLS, figsize=(HIGHLIGHT_NCOLS * 3, nrows * 3))
axes = axes.ravel()
sc.pl.umap(combined, color=BATCH_KEY, ax=axes[0], size=6, legend_loc='on data',
           legend_fontsize=6, show=False)
axes[0].set_title('all')
for ax, g in zip(axes[1:], groups):
    sc.pl.umap(combined, color=BATCH_KEY, groups=[g], ax=ax, size=6, legend_loc=None,
               na_in_legend=False, show=False)
    ax.set_title(str(g))
for ax in axes[len(groups) + 1:]:
    ax.axis('off')
fig.suptitle(BATCH_KEY)
fig.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'umap_batch_highlights.{ext}'), dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved umap_batch_highlights to {FIG_DIR}')

# ── Leiden resolution sweep + Wilcoxon marker tables ──────────────────────────
sweep_keys = []
for res in RESOLUTIONS:
    key = f'leiden_res{res:.2f}'
    sweep_keys.append(key)
    # Adds obs[key]: a cluster number per cell at this resolution. Cluster numbers are only
    # meaningful within this run and resolution.
    sc.tl.leiden(combined, resolution=res, flavor='igraph', n_iterations=2,
                 directed=False, random_state=0, key_added=key)
    print(f'resolution {res:.2f}: {combined.obs[key].nunique()} clusters')
    # Adds uns['rank_{key}']: marker-gene statistics for each cluster, also written to a CSV below.
    sc.tl.rank_genes_groups(combined, key, method='wilcoxon', key_added=f'rank_{key}')
    deg = sc.get.rank_genes_groups_df(combined, group=None, key=f'rank_{key}')
    deg.to_csv(os.path.join(SWEEP_DIR, f'deg_{key}.csv'), index=False)
print(f'\nWrote {len(RESOLUTIONS)} DEG tables to {SWEEP_DIR}')

# ── Resolution-sweep UMAP grid ────────────────────────────────────────────────
fig = sc.pl.umap(combined, color=sweep_keys, ncols=3, wspace=0.05, size=6,
                 legend_loc='on data', show=False, return_fig=True)
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'umap_resolution_sweep.{ext}'), dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved umap_resolution_sweep to {FIG_DIR}')

# ── Per-resolution cluster UMAPs + cluster-highlight grids ─────────────────────
for key in sweep_keys:
    fig = sc.pl.umap(combined, color=key, legend_loc='on data', legend_fontsize=12,
                     size=6, show=False, return_fig=True)
    fig.axes[0].set_title(key)
    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(CLUSTERS_DIR, f'umap_{key}.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)

    clusters = list(combined.obs[key].cat.categories)
    nrows = math.ceil((len(clusters) + 1) / HIGHLIGHT_NCOLS)
    fig, axes = plt.subplots(nrows, HIGHLIGHT_NCOLS, figsize=(HIGHLIGHT_NCOLS * 3, nrows * 3))
    axes = axes.ravel()
    sc.pl.umap(combined, color=key, ax=axes[0], size=6, legend_loc='on data',
               legend_fontsize=6, show=False)
    axes[0].set_title('all')
    for ax, cluster in zip(axes[1:], clusters):
        sc.pl.umap(combined, color=key, groups=[cluster], ax=ax, size=6, legend_loc=None,
                   na_in_legend=False, show=False)
        ax.set_title(str(cluster))
    for ax in axes[len(clusters) + 1:]:
        ax.axis('off')
    fig.suptitle(key)
    fig.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(HIGHLIGHTS_DIR, f'umap_{key}_highlights.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)
print(f'Saved per-resolution UMAPs and highlight grids to {FIG_DIR}')

# Saves everything added above (X_pca, X_pca_harmony, X_umap, neighbor graph, leiden columns,
# rank_* results) as a new file; .X is log-normalized and layers['counts'] holds raw counts.
out_path = os.path.join(OBJECTS_DIR, 'combined_clustered.h5ad')
combined.write(out_path)
print(f'\nWrote {out_path}')
