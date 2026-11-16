import os
import warnings

import scanpy as sc

warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

DATA_DIR = '/home/workspace/PCD_data/xenium/commercial/data'
OBJECTS_DIR = os.path.join(DATA_DIR, 'objects')

KEY = 'leiden_res0.60'

# WARNING: this mapping is valid only for the specific clustering in
# combined_clustered.h5ad. Leiden cluster numbers are not stable — if
# 03_processing.py is re-run, re-read the new marker-gene plots and
# rebuild this dict from scratch before running this script.
cluster_to_cell_type_res06 = {
    "0":  "Erythroid Precursor 1",
    "1":  "Erythroid Precursor 2",
    "2":  "Neoplastic Plasma Cell 1 (MYC/CXCR4-high)",
    "3":  "Neoplastic Plasma Cell 2 (LILRB4-high)",
    "4":  "Endothelial",
    "5":  "Stromal Fibroblast",
    "6":  "Neoplastic Plasma Cell 3",
    "7":  "Cytotoxic T Cell",
    "8":  "Macrophage",
    "9":  "Monocyte",
    "10": "Proliferating Myeloid Precursor",
}

combined = sc.read_h5ad(os.path.join(OBJECTS_DIR, 'combined_clustered.h5ad'))
print(f'Loaded: {combined.n_obs:,} cells')

missing = set(combined.obs[KEY].unique()) - set(cluster_to_cell_type_res06.keys())
if missing:
    raise ValueError(
        f'Clusters in {KEY} not covered by mapping: {sorted(missing)}. '
        'Re-run 03_processing.py and rebuild cluster_to_cell_type_res06.'
    )

# Adds obs['cell_type_res06']: each cell gets the name its leiden_res0.60 cluster maps to.
# No other data changes.
combined.obs['cell_type_res06'] = (
    combined.obs[KEY].map(cluster_to_cell_type_res06).astype('category')
)
print('Cell-type counts (res 0.6):')
print(combined.obs['cell_type_res06'].value_counts().to_string())

out_path = os.path.join(OBJECTS_DIR, 'combined_clustered.h5ad')
# Overwrites combined_clustered.h5ad in place with the new cell_type_res06 column included.
combined.write_h5ad(out_path)
print(f'Updated {out_path}')
