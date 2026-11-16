import json
import os
import warnings

import matplotlib.pyplot as plt
import scanpy as sc

warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

DATA_DIR = '/home/workspace/PCD_data/xenium/commercial/data'
OBJECTS_DIR = os.path.join(DATA_DIR, 'objects')
FIG_DIR  = '/home/workspace/PCD_data/xenium/commercial/figures/spatial_celltypes'
os.makedirs(FIG_DIR, exist_ok=True)

combined = sc.read_h5ad(os.path.join(OBJECTS_DIR, 'combined_clustered.h5ad'))
print(f'Loaded: {combined.n_obs:,} cells')

palette_path = os.path.join(OBJECTS_DIR, 'celltypes_palette.json')
with open(palette_path) as f:
    palettes = json.load(f)
# Overwrites uns['{col}_colors'] in memory with the palette saved by 06, so the colors match
# the UMAP. The h5ad on disk is not rewritten.
for col, colors in palettes.items():
    combined.uns[f'{col}_colors'] = colors

ANNOTATION_COLS = [
    ('cell_type_res1_1', 'res1_1'),
]

for sample in sorted(combined.obs['sample'].unique()):
    sub = combined[combined.obs['sample'] == sample].copy()
    # Adds obsm['X_spatial'] to this per-sample copy only: the same coordinates as obsm['spatial']
    # under the name sc.pl.embedding expects (it looks for X_{basis}). Values are unchanged.
    sub.obsm['X_spatial'] = sub.obsm['spatial']

    for col, label in ANNOTATION_COLS:
        fig = sc.pl.embedding(
            sub,
            basis='spatial',
            color=col,
            legend_loc='right margin',
            show=False,
            return_fig=True,
            size=6,
        )
        for ax in fig.axes:
            for coll in ax.collections:
                coll.set_rasterized(True)
        fname = f'spatial_celltypes_{label}_{sample}'
        for ext in ['png', 'pdf']:
            fig.savefig(os.path.join(FIG_DIR, f'{fname}.{ext}'),
                        dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f'  Saved {fname}')

print(f'Done. Figures saved to {FIG_DIR}')
