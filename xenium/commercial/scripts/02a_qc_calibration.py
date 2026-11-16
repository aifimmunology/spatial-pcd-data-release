import logging
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

DATA_DIR    = '/home/workspace/PCD_data/xenium/commercial/data'
SAMPLES_DIR = os.path.join(DATA_DIR, 'samples')
LOGS_DIR    = os.path.join(DATA_DIR, 'logs')
FIGURES_DIR = '/home/workspace/PCD_data/xenium/commercial/figures/qc_calibration'

for directory in (LOGS_DIR, FIGURES_DIR):
    os.makedirs(directory, exist_ok=True)

logging.basicConfig(
    filename=os.path.join(LOGS_DIR, '02a_qc_calibration.log'),
    level=logging.INFO,
    format='%(asctime)s %(message)s',
)

SAMPLE_NAMES = [
    'TIS05393',
    'TIS04788',
    'TIS05392',
    'TIS05394',
    'TIS04789',
    'TIS05397_008',
]

QC_COLUMNS = ['total_counts', 'n_genes_by_counts']


def save_figure(fig, name):
    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(FIGURES_DIR, f'{name}.{ext}'), dpi=150, bbox_inches='tight')
    plt.close(fig)
    logging.info('saved %s to %s', name, FIGURES_DIR)


samples = {}
for name in SAMPLE_NAMES:
    adata = sc.read_h5ad(os.path.join(SAMPLES_DIR, f'{name}.h5ad'))

    # Adds per-cell obs columns total_counts and n_genes_by_counts (in memory only; nothing is saved).
    sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
    samples[name] = adata

# One row per sample: counts histogram, genes histogram, and knee plot (rank vs counts).
# A sharp bend in the knee marks a natural cutoff between real cells and debris.
fig, axes = plt.subplots(len(samples), 3, figsize=(15, 3.2 * len(samples)), squeeze=False)
for row_axes, (name, adata) in zip(axes, samples.items()):
    counts = adata.obs['total_counts'].to_numpy()
    genes  = adata.obs['n_genes_by_counts'].to_numpy()

    row_axes[0].hist(np.log10(counts + 1), bins=60, color='steelblue')
    row_axes[0].set_xlabel('log10(counts + 1)')
    row_axes[0].set_title(f'{name} — counts  (median {np.median(counts):.0f})')

    row_axes[1].hist(np.log10(genes + 1), bins=60, color='steelblue')
    row_axes[1].set_xlabel('log10(genes + 1)')
    row_axes[1].set_title(f'{name} — genes  (median {np.median(genes):.0f})')

    row_axes[2].plot(np.arange(1, len(counts) + 1), np.sort(counts)[::-1], lw=1.5)
    row_axes[2].set_xscale('log')
    row_axes[2].set_yscale('log')
    row_axes[2].set_xlabel('cell rank')
    row_axes[2].set_ylabel('counts')
    row_axes[2].set_title(f'{name} — knee')
fig.tight_layout()
save_figure(fig, 'qc_distributions')

# One row per sample, one column per QC metric, plotted at tissue coordinates (color capped at p99).
# Low counts at tissue edges are normal; a diffuse low-count haze points to a run problem.
fig, axes = plt.subplots(
    len(samples), len(QC_COLUMNS),
    figsize=(5.5 * len(QC_COLUMNS), 5 * len(samples)),
    squeeze=False,
)
for i, (name, adata) in enumerate(samples.items()):
    for j, column in enumerate(QC_COLUMNS):
        sc.pl.spatial(
            adata, color=column, vmax='p99', spot_size=50,
            ax=axes[i, j], show=False,
            title=f'{name} — {column}  (median {adata.obs[column].median():.0f})',
        )
        for coll in axes[i, j].collections:
            coll.set_rasterized(True)
fig.tight_layout()
save_figure(fig, 'spatial_qc_maps')

summary = pd.DataFrame({
    name: {
        'n_cells':       adata.n_obs,
        'median_counts': int(np.median(adata.obs['total_counts'])),
        'median_genes':  int(np.median(adata.obs['n_genes_by_counts'])),
    }
    for name, adata in samples.items()
}).T
logging.info('per-sample QC summary:\n%s', summary.to_string())
