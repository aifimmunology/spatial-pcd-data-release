import logging
import os

import matplotlib.pyplot as plt
import scipy.sparse as sp
import scanpy as sc
import harmonypy as hm

DATA_DIR    = '/home/workspace/PCD_data/xenium/commercial/data'
OBJECTS_DIR = os.path.join(DATA_DIR, 'objects')
SWEEP_DIR   = os.path.join(DATA_DIR, 'resolution_sweep')
LOGS_DIR    = os.path.join(DATA_DIR, 'logs')
FIGURES_DIR = '/home/workspace/PCD_data/xenium/commercial/figures/processing'

for directory in (OBJECTS_DIR, SWEEP_DIR, LOGS_DIR, FIGURES_DIR):
    os.makedirs(directory, exist_ok=True)

logging.basicConfig(
    filename=os.path.join(LOGS_DIR, '03_processing.log'),
    level=logging.INFO,
    format='%(asctime)s %(message)s',
)

RESOLUTIONS = [round(0.1 * i, 2) for i in range(1, 12)]  # 0.1, 0.2, …, 1.1


def save_figure(fig, name):
    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(FIGURES_DIR, f'{name}.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    logging.info('saved %s to %s', name, FIGURES_DIR)


combined = sc.read_h5ad(os.path.join(OBJECTS_DIR, 'combined_filtered.h5ad'))

# Overwrites .X: rescales every cell to the dataset's median total counts (default target_sum).
sc.pp.normalize_total(combined, inplace=True)
# Changes .X's storage format only (dense to sparse), not its values; most of the 377 genes
# have zero counts in any given cell.
combined.X = sp.csr_matrix(combined.X)
# log1p-transforms .X in place (values become log(count+1)); .X comes from normalize_total above.
sc.pp.log1p(combined)
# Adds layers['lognorm']: snapshot of .X after log1p, before scaling (log-normalized, unscaled).
combined.layers["lognorm"] = combined.X.copy()
# Clips .X values to max 10 and scales each gene to unit variance; zero_center=False keeps .X sparse.
sc.pp.scale(combined, zero_center=False, max_value=10)
# Adds obsm['X_pca']: 30 principal components computed from all 377 genes.
sc.pp.pca(combined, n_comps=30, use_highly_variable=False, random_state=0)

# Harmony shifts the PCA coordinates so the six samples (batches) overlap.
harmony = hm.run_harmony(
    combined.obsm['X_pca'].astype('float64'),
    combined.obs,
    'sample',
    random_state=0,
)
# harmonypy versions differ on whether Z_corr is cells x PCs or PCs x cells.
harmony_pca = harmony.Z_corr
if harmony_pca.shape[0] != combined.n_obs:
    harmony_pca = harmony_pca.T
# Adds obsm['X_pca_harmony']: the batch-corrected coordinates. obsm['X_pca'] stays uncorrected.
combined.obsm['X_pca_harmony'] = harmony_pca

# Adds obsp['connectivities'] and obsp['distances'] (each cell's 15 nearest neighbors, measured
# in the Harmony space) and uns['neighbors'] (the settings used).
sc.pp.neighbors(combined, use_rep='X_pca_harmony', n_neighbors=15, n_pcs=30, random_state=0)
# Adds obsm['X_umap']: a 2D position per cell computed from the neighbor graph, plus uns['umap'].
sc.tl.umap(combined, random_state=0)

# Cells should intermix by sample; strong per-sample separation means Harmony did not converge.
fig = sc.pl.umap(combined, color=['total_counts', 'n_genes_by_counts', 'sample'],
                 wspace=0.2, show=False, return_fig=True,
                 size=6, rasterized=True)
save_figure(fig, 'umap_harmony_diagnostic')

sweep_keys = [f'leiden_res{res:.2f}' for res in RESOLUTIONS]
for res, key in zip(RESOLUTIONS, sweep_keys):
    # Adds an obs column (leiden_res{res}) assigning each cell to a cluster at this resolution.
    sc.tl.leiden(combined, flavor='igraph', n_iterations=2, directed=False,
                 resolution=res, key_added=key, random_state=0)
    logging.info('resolution %.2f: %d clusters', res, combined.obs[key].nunique())

    # Adds uns['rank_{key}']: marker-gene statistics for each cluster, then writes them to a CSV.
    sc.tl.rank_genes_groups(combined, key, method='wilcoxon', key_added=f'rank_{key}')
    deg = sc.get.rank_genes_groups_df(combined, group=None, key=f'rank_{key}')
    deg_path = os.path.join(SWEEP_DIR, f'deg_{key}_wilcoxon.csv')
    deg.to_csv(deg_path, index=False)
    logging.info('wrote %s', deg_path)

fig = sc.pl.umap(combined, color=sweep_keys, wspace=0.05, ncols=3,
                 legend_loc='on data', show=False, return_fig=True,
                 size=6, rasterized=True)
save_figure(fig, 'umap_resolution_sweep')

out_path = os.path.join(OBJECTS_DIR, 'combined_clustered.h5ad')
# Saves everything added above (X_pca, X_pca_harmony, X_umap, neighbor graph, leiden columns,
# rank_* results) as a new file; .X holds normalized values and layers['counts'] the raw counts.
combined.write_h5ad(out_path)
logging.info('wrote %s', out_path)
