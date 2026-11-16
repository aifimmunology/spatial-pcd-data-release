"""
fh1_pretreatment AZU1-removed run — Wilcoxon DEGs for the Pearson-residual clusters.

Pearson residuals (X) only drove the clustering. DE runs on log-normalized counts
rebuilt from layers['counts'], identical to combined_clustered.X, so these tables
compare directly with resolution_sweep/deg_*_wilcoxon.csv; only the clusters differ.

    python xenium/fh1_pretreatment/scripts/archive/azu1_removed/08_pearson_degs.py
"""

import os
import warnings

import scanpy as sc

warnings.filterwarnings('ignore')

RUN_DIR = '/home/workspace/PCD_data/xenium/fh1_pretreatment/data/archive/azu1_removed'
OUT_DIR = os.path.join(RUN_DIR, 'resolution_sweep_pearson')
os.makedirs(OUT_DIR)  # no exist_ok: refuse to write into an existing run's folder

RESOLUTIONS = [round(0.1 * i, 2) for i in range(1, 12)]  # 0.1, 0.2, …, 1.1

# backed='r' leaves the 1.6 GB dense residual matrix on disk; only counts/obs/var are loaded.
pear = sc.read_h5ad(os.path.join(RUN_DIR, 'combined_pearson_clustered.h5ad'), backed='r')
lognorm = sc.AnnData(pear.layers['counts'].copy(), obs=pear.obs.copy(), var=pear.var.copy())
sc.pp.normalize_total(lognorm)
sc.pp.log1p(lognorm)

for res in RESOLUTIONS:
    key = f'leiden_res{res:.2f}'
    sc.tl.rank_genes_groups(lognorm, key, method='wilcoxon', use_raw=False,
                            key_added=f'rank_{key}')
    deg = sc.get.rank_genes_groups_df(lognorm, group=None, key=f'rank_{key}')
    deg.to_csv(os.path.join(OUT_DIR, f'deg_{key}_wilcoxon.csv'), index=False)

print(f'Wrote {len(RESOLUTIONS)} DEG tables to {OUT_DIR}')
