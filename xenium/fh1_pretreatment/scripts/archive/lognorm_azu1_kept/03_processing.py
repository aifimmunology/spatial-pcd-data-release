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
# False), so this flag is the switch — and passing rasterized= to sc.pl.umap
# directly raises a TypeError. Setting the flag rather than set_figure_params()
# avoids that function's font-size/dpi side effects.
sc.settings._vector_friendly = True

DATA_DIR      = '/home/workspace/PCD_data/xenium/fh1_pretreatment/data/archive/lognorm_azu1_kept'
OBJECTS_DIR   = os.path.join(DATA_DIR, 'objects')
FIGURES_DIR   = '/home/workspace/PCD_data/xenium/fh1_pretreatment/figures/archive/lognorm_azu1_kept'
SWEEP_DIR       = os.path.join(DATA_DIR, 'resolution_sweep')
SWEEP_FIG_DIR   = os.path.join(FIGURES_DIR, 'resolution_sweep')
CLUSTERS_FIG_DIR   = os.path.join(SWEEP_FIG_DIR, 'clusters')
HIGHLIGHTS_FIG_DIR = os.path.join(SWEEP_FIG_DIR, 'highlights')
os.makedirs(OBJECTS_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(SWEEP_DIR, exist_ok=True)
os.makedirs(CLUSTERS_FIG_DIR, exist_ok=True)
os.makedirs(HIGHLIGHTS_FIG_DIR, exist_ok=True)

RESOLUTIONS = [round(0.1 * i, 2) for i in range(1, 12)]  # 0.1, 0.2, …, 1.1
HIGHLIGHT_NCOLS = 4

combined = sc.read_h5ad(os.path.join(OBJECTS_DIR, 'combined_filtered.h5ad'))
print(f'Loaded: {combined.n_obs:,} cells x {combined.n_vars} genes')

# ── PCA ───────────────────────────────────────────────────────────────────────
sc.pp.pca(combined, n_comps=30, use_highly_variable=False, random_state=0)

# ── Harmony batch correction ──────────────────────────────────────────────────
# Harmony adjusts the 30 PCA components to remove per-patient technical variation
# while preserving biological signal. We correct by 'cohort_subject' (patient):
# each patient's tissue shares biopsy, processing, and Xenium-run effects. All 13
# samples are the same MM Pre-Treatment timepoint, so patient is the batch here.
ho = hm.run_harmony(
    combined.obsm['X_pca'].astype('float64'),
    combined.obs,
    'cohort_subject',
    random_state=0,
)

# Z_corr shape may be (n_pcs x n_cells) or (n_cells x n_pcs) depending on version
z = ho.Z_corr
if z.shape[0] != combined.n_obs:
    z = z.T
combined.obsm['X_pca_harmony'] = z
combined.obsm['X_pca'] = z  # downstream steps use X_pca by default

sc.pp.neighbors(combined, use_rep='X_pca', n_neighbors=15, n_pcs=30, random_state=0)
sc.tl.umap(combined, min_dist=0.3, random_state=0)

# ── Diagnostic UMAP: verify batch correction ──────────────────────────────────
# Cells should intermix by patient after Harmony; strong per-patient separation
# here would mean the correction did not converge.
fig = sc.pl.umap(combined,
                 color=['total_counts', 'n_genes_by_counts', 'cohort_subject', 'treatment_timepoint'],
                 wspace=0.3, ncols=2, show=False, return_fig=True,
                 size=6)
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIGURES_DIR, f'umap_harmony_diagnostic.{ext}'),
                dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved umap_harmony_diagnostic to {FIGURES_DIR}')

# ── Per-subject highlight grid ────────────────────────────────────────────────
# Same contact-sheet idea as the per-cluster grids below, but split by patient
# instead of by cluster. The first panel colours every patient; each following
# panel isolates one patient — that patient's cells coloured, every other cell
# greyed out — to show where each patient sits in the shared embedding after
# Harmony. groups= does the greying (non-selected cells render as light-grey).
combined.obs['cohort_subject'] = combined.obs['cohort_subject'].astype('category')
subjects = list(combined.obs['cohort_subject'].cat.categories)
nrows = math.ceil((len(subjects) + 1) / HIGHLIGHT_NCOLS)
fig, axes = plt.subplots(nrows, HIGHLIGHT_NCOLS,
                         figsize=(HIGHLIGHT_NCOLS * 3, nrows * 3))
axes = axes.ravel()

sc.pl.umap(combined, color='cohort_subject', ax=axes[0], size=6,
           legend_loc='right margin', legend_fontsize=6, show=False)
axes[0].set_title('all subjects')

for ax, subject in zip(axes[1:], subjects):
    sc.pl.umap(combined, color='cohort_subject', groups=[subject], ax=ax, size=6,
               legend_loc=None, na_in_legend=False, show=False)
    ax.set_title(subject)

for ax in axes[len(subjects) + 1:]:
    ax.axis('off')

fig.suptitle('cohort_subject')
fig.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIGURES_DIR, f'umap_subject_highlights.{ext}'),
                dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved umap_subject_highlights to {FIGURES_DIR}')

# ── Resolution sweep ──────────────────────────────────────────────────────────
sweep_keys = []
for res in RESOLUTIONS:
    key = f'leiden_res{res:.2f}'
    sweep_keys.append(key)

    sc.tl.leiden(combined, flavor='igraph', n_iterations=2, directed=False,
                 resolution=res, key_added=key, random_state=0)

    n_clusters = combined.obs[key].nunique()
    print(f'resolution {res:.2f}: {n_clusters} clusters')

    sc.tl.rank_genes_groups(combined, key, method='wilcoxon',
                            key_added=f'rank_{key}')
    deg = sc.get.rank_genes_groups_df(combined, group=None, key=f'rank_{key}')
    deg.to_csv(os.path.join(SWEEP_DIR, f'deg_{key}_wilcoxon.csv'), index=False)

print(f'\nWrote {len(RESOLUTIONS)} DEG tables to {SWEEP_DIR}')

# ── Per-resolution UMAPs ──────────────────────────────────────────────────────
# One standalone UMAP per resolution, colored by that resolution's clusters, to
# accompany each DEG table. The cluster columns already exist in .obs from the
# sweep above, so each plot just recolors the precomputed UMAP.
for key in sweep_keys:
    fig = sc.pl.umap(combined, color=key, legend_loc='on data',
                     legend_fontsize=12, size=6, show=False, return_fig=True)
    fig.axes[0].set_title(key)
    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(CLUSTERS_FIG_DIR, f'umap_{key}.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)
print(f'Saved {len(sweep_keys)} per-resolution UMAPs to {CLUSTERS_FIG_DIR}')

# ── Per-resolution cluster-highlight grids ────────────────────────────────────
# For each resolution, a contact sheet of small UMAPs. The first panel is the
# full clustering (every cluster colored) as a reference; each following panel
# isolates one cluster — that cluster's cells colored, every other cell greyed
# out. scanpy's groups= argument does the greying (non-selected cells render as
# light-grey 'NA'); we build the axis grid ourselves because a single sc.pl.umap
# call grids different variables per panel, not different highlights of one
# variable.
for key in sweep_keys:
    clusters = list(combined.obs[key].cat.categories)
    nrows = math.ceil((len(clusters) + 1) / HIGHLIGHT_NCOLS)
    fig, axes = plt.subplots(nrows, HIGHLIGHT_NCOLS,
                             figsize=(HIGHLIGHT_NCOLS * 3, nrows * 3))
    axes = axes.ravel()

    sc.pl.umap(combined, color=key, ax=axes[0], size=6,
               legend_loc='on data', legend_fontsize=6, show=False)
    axes[0].set_title('all clusters')

    for ax, cluster in zip(axes[1:], clusters):
        sc.pl.umap(combined, color=key, groups=[cluster], ax=ax, size=6,
                   legend_loc=None, na_in_legend=False, show=False)
        ax.set_title(f'cluster {cluster}')

    for ax in axes[len(clusters) + 1:]:
        ax.axis('off')

    fig.suptitle(key)
    fig.tight_layout()
    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(HIGHLIGHTS_FIG_DIR, f'umap_{key}_highlights.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)
print(f'Saved {len(sweep_keys)} cluster-highlight grids to {HIGHLIGHTS_FIG_DIR}')

# ── Resolution sweep UMAP grid ────────────────────────────────────────────────
fig = sc.pl.umap(combined, color=sweep_keys, wspace=0.05, ncols=3,
                 legend_loc='on data', show=False, return_fig=True,
                 size=6)
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIGURES_DIR, f'umap_resolution_sweep.{ext}'),
                dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved umap_resolution_sweep to {FIGURES_DIR}')

# ── Save ──────────────────────────────────────────────────────────────────────
out_path = os.path.join(OBJECTS_DIR, 'combined_clustered.h5ad')
combined.write_h5ad(out_path)
print(f'\nWrote {out_path}')
