import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc

warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

DATA_DIR    = '/home/workspace/PCD_data/xenium/fh1_pretreatment/data/archive/lognorm_azu1_kept'
OBJECTS_DIR = os.path.join(DATA_DIR, 'objects')
FIGURES_DIR = '/home/workspace/PCD_data/xenium/fh1_pretreatment/figures/archive/lognorm_azu1_kept'
os.makedirs(FIGURES_DIR, exist_ok=True)

combined = sc.read_h5ad(os.path.join(OBJECTS_DIR, 'combined_raw.h5ad'))
sc.pp.calculate_qc_metrics(combined, percent_top=None, inplace=True)
# n_genes_by_counts: number of genes with at least one count in the cell
combined.obs['n_genes_by_counts'] = np.asarray((combined.X > 0).sum(axis=1)).flatten()
print(f'Loaded: {combined.n_obs:,} cells x {combined.n_vars} genes')

# ── Per-sample distribution plots, one figure per timepoint ───────────────────
# Three panels per sample: log10 count histogram, log10 gene histogram, knee plot.
# The knee plot (rank vs counts on log-log axes) shows where sequencing depth
# falls off — a sharp bend indicates a natural threshold between real cells and
# empty/debris captures.

for timepoint in sorted(combined.obs['treatment_timepoint'].unique()):
    tp_cells = combined[combined.obs['treatment_timepoint'] == timepoint]
    sample_names = sorted(tp_cells.obs['sample'].unique())

    fig, axes = plt.subplots(len(sample_names), 3,
                             figsize=(15, 3.2 * len(sample_names)), squeeze=False)

    for row_axes, name in zip(axes, sample_names):
        adata = tp_cells[tp_cells.obs['sample'] == name]
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

    fig.suptitle(timepoint, fontsize=14, y=1.001)
    fig.tight_layout()
    slug = timepoint.replace(' ', '_').replace('/', '_')
    for ext in ['png', 'pdf']:
        fig.savefig(os.path.join(FIGURES_DIR, f'qc_distributions_{slug}.{ext}'),
                    dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved qc_distributions_{slug}')

# ── Summary table ─────────────────────────────────────────────────────────────
summary = combined.obs.groupby('sample', observed=True).agg(
    n_cells=('total_counts', 'size'),
    median_counts=('total_counts', 'median'),
    median_genes=('n_genes_by_counts', 'median'),
)
summary = summary.join(
    combined.obs.groupby('sample', observed=True)[['cohort_subject', 'treatment_timepoint']].first()
)
summary = summary.sort_values(['treatment_timepoint', 'cohort_subject'])
print('\nPer-sample QC summary:')
print(summary.to_string())
