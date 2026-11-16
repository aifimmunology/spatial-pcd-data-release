import json
import os
import warnings

import matplotlib.pyplot as plt
import scanpy as sc

warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

DATA_DIR = '/home/workspace/PCD_data/xenium/commercial/data'
OBJECTS_DIR = os.path.join(DATA_DIR, 'objects')
FIG_DIR  = '/home/workspace/PCD_data/xenium/commercial/figures/umap_celltypes'
os.makedirs(FIG_DIR, exist_ok=True)

combined = sc.read_h5ad(os.path.join(OBJECTS_DIR, 'combined_clustered.h5ad'))
print(f'Loaded: {combined.n_obs:,} cells')

ANNOTATION_COLS = [
    ('cell_type_res1_1', 'res1_1'),
]

for col, label in ANNOTATION_COLS:
    fig = sc.pl.umap(
        combined,
        color=col,
        legend_loc='on data',
        legend_fontsize=8,
        show=False,
        return_fig=True,
        size=6,
    )
    for ax in fig.axes:
        for coll in ax.collections:
            coll.set_rasterized(True)
    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(FIG_DIR, f'umap_celltypes_{label}_ondata.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)

    fig = sc.pl.umap(
        combined,
        color=col,
        legend_loc='right margin',
        show=False,
        return_fig=True,
        size=6,
    )
    for ax in fig.axes:
        for coll in ax.collections:
            coll.set_rasterized(True)
    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(FIG_DIR, f'umap_celltypes_{label}_legend.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f'Saved {label} UMAP figures to {FIG_DIR}')

# Plotting a categorical column made scanpy add uns['{col}_colors'] in memory. That is the only
# change to the object; it is saved to a separate JSON for 07, and the h5ad is not rewritten.
palette = {col: list(combined.uns[f'{col}_colors']) for col, _ in ANNOTATION_COLS if f'{col}_colors' in combined.uns}
palette_path = os.path.join(OBJECTS_DIR, 'celltypes_palette.json')
with open(palette_path, 'w') as f:
    json.dump(palette, f, indent=2)
print(f'Saved color palette to {palette_path}')
