import os
import warnings

import scanpy as sc

warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

DATA_DIR = '/home/workspace/PCD_data/xenium/commercial/data'
OBJECTS_DIR = os.path.join(DATA_DIR, 'objects')

KEY = 'leiden_res1.10'

# WARNING: this mapping is valid only for the specific clustering in
# combined_clustered.h5ad. Leiden cluster numbers are not stable — if
# 03_processing.py is re-run, re-read the new marker-gene plots and
# rebuild this dict from scratch before running this script.
cluster_to_cell_type_res1_1 = {
    "0":  "Erythroid Precursor 1",
    "1":  "Erythroid Precursor 2",
    "2":  "Erythroid Precursor 2",
    "3":  "Neoplastic Plasma Cell 2 (BCL2L11-high, apoptosis-primed)",
    "4":  "Neoplastic Plasma Cell 3",
    "5":  "Neoplastic Plasma Cell 6",
    "6":  "Stromal Fibroblast",
    "7":  "Neoplastic Plasma Cell 2 (BCL2L11-high, apoptosis-primed)",
    "8":  "Proliferating Myeloid Precursor",
    "9":  "Macrophage",
    "10": "Monocyte",
    "11": "Macrophage",
    "12": "Cytotoxic T Cell",
    "13": "Neoplastic Plasma Cell 5 (MYC-high, possibly proliferative)",
    "14": "Neoplastic Plasma Cell 4",
    "15": "Neoplastic Plasma Cell 1",
    "16": "Endothelial",
    "17": "TREM2+ Macrophage (Tumor-Supportive, Adjacent)",
}

combined = sc.read_h5ad(os.path.join(OBJECTS_DIR, 'combined_clustered.h5ad'))
print(f'Loaded: {combined.n_obs:,} cells')

missing = set(combined.obs[KEY].unique()) - set(cluster_to_cell_type_res1_1.keys())
if missing:
    raise ValueError(
        f'Clusters in {KEY} not covered by mapping: {sorted(missing)}. '
        'Re-run 03_processing.py and rebuild cluster_to_cell_type_res1_1.'
    )

# Adds obs['cell_type_res1_1']: each cell gets the name its leiden_res1.10 cluster maps to.
# Several clusters share a name, so the 18 clusters become 15 labels. No other data changes.
combined.obs['cell_type_res1_1'] = (
    combined.obs[KEY].map(cluster_to_cell_type_res1_1).astype('category')
)
print('Cell-type counts (res 1.1):')
print(combined.obs['cell_type_res1_1'].value_counts().to_string())

out_path = os.path.join(OBJECTS_DIR, 'combined_clustered.h5ad')
# Overwrites combined_clustered.h5ad in place with the new cell_type_res1_1 column included.
combined.write_h5ad(out_path)
print(f'Updated {out_path}')
