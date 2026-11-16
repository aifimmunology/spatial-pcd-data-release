import os
import warnings

import matplotlib.pyplot as plt
import scanpy as sc

warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

DATA_DIR = '/home/workspace/PCD_data/xenium/commercial/data'
OBJECTS_DIR = os.path.join(DATA_DIR, 'objects')
FIG_DIR  = '/home/workspace/PCD_data/xenium/commercial/figures/umap_by_sample'
os.makedirs(FIG_DIR, exist_ok=True)

combined = sc.read_h5ad(os.path.join(OBJECTS_DIR, 'combined_clustered.h5ad'))
print(f'Loaded: {combined.n_obs:,} cells')

fig = sc.pl.umap(
    combined,
    color='sample',
    legend_loc='right margin',
    show=False,
    return_fig=True,
    size=6,
)
for ax in fig.axes:
    for coll in ax.collections:
        coll.set_rasterized(True)
for ext in ['png', 'pdf']:
    fig.savefig(os.path.join(FIG_DIR, f'umap_by_sample.{ext}'),
                dpi=150, bbox_inches='tight')
plt.close(fig)
print(f'Saved umap_by_sample to {FIG_DIR}')
