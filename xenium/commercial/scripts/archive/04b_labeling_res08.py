import os
import warnings

import scanpy as sc

warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

DATA_DIR = '/home/workspace/PCD_data/xenium/commercial/data'
OBJECTS_DIR = os.path.join(DATA_DIR, 'objects')

KEY = 'leiden_res0.80'

# WARNING: this mapping is valid only for the specific clustering in
# combined_clustered.h5ad. Leiden cluster numbers are not stable — if
# 03_processing.py is re-run, re-read the new marker-gene plots and
# rebuild this dict from scratch before running this script.
cluster_to_cell_type_res08 = {
    "0":  "Erythroid Precursor 1",
    "1":  "Erythroid Precursor 2",
    "2":  "Endothelial",
    "3":  "Neoplastic Plasma Cell 1",
    "4":  "Neoplastic Plasma Cell 2 (BCL2L11-high, apoptosis-primed)",
    "5":  "Neoplastic Plasma Cell 3",
    "6":  "TREM2+ Macrophage (Tumor-Supportive, Adjacent)",
    "7":  "Monocyte",
    "8":  "Neoplastic Plasma Cell 4",
    "9":  "Macrophage",
    "10": "Cytotoxic T Cell",
    "11": "Neoplastic Plasma Cell 5 (MYC-high, possibly proliferative)",
    "12": "Plasmablast (B-to-Plasma Cell Transitional)",
    "13": "Stromal Fibroblast",
    "14": "Neoplastic Plasma Cell 6",
    "15": "Lymphocyte/Dendritic Cell/Vessel (Intravascular)",
    "16": "Proliferating Myeloid Precursor",
    "17": "Lipid-Laden Macrophage (Tumor-Adjacent)",
}

combined = sc.read_h5ad(os.path.join(OBJECTS_DIR, 'combined_clustered.h5ad'))
print(f'Loaded: {combined.n_obs:,} cells')

missing = set(combined.obs[KEY].unique()) - set(cluster_to_cell_type_res08.keys())
if missing:
    raise ValueError(
        f'Clusters in {KEY} not covered by mapping: {sorted(missing)}. '
        'Re-run 03_processing.py and rebuild cluster_to_cell_type_res08.'
    )

combined.obs['cell_type_res08'] = (
    combined.obs[KEY].map(cluster_to_cell_type_res08).astype('category')
)
print('Cell-type counts (res 0.8):')
print(combined.obs['cell_type_res08'].value_counts().to_string())

out_path = os.path.join(OBJECTS_DIR, 'combined_clustered.h5ad')
combined.write_h5ad(out_path)
print(f'Updated {out_path}')
