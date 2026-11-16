"""
fh1_pretreatment AZU1-removed run — QC filtering.

Same as xenium/fh1_pretreatment/scripts/archive/lognorm_azu1_kept/02b_qc_filtering.py, with two differences:
  1. AZU1 is dropped (QC flagged it as background / off-target signal).
  2. Outputs go to the separate objects_azu1_removed/ folder; the original run's
     objects/ is read-only here (only combined_raw.h5ad is read) and never written.
"""

import os
import warnings

import pandas as pd
import scipy.sparse as sp
import scanpy as sc

warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

DATA_DIR = '/home/workspace/PCD_data/xenium/fh1_pretreatment/data/archive'
SRC_OBJECTS_DIR = os.path.join(DATA_DIR, 'lognorm_azu1_kept', 'objects')            # read-only: source of combined_raw
RUN_DIR = os.path.join(DATA_DIR, 'azu1_removed')       # new run outputs
os.makedirs(RUN_DIR, exist_ok=True)

MIN_COUNTS = 15
MIN_GENES  = 10

combined = sc.read_h5ad(os.path.join(SRC_OBJECTS_DIR, 'combined_raw.h5ad'))
print(f'Loaded: {combined.n_obs:,} cells x {combined.n_vars} genes')

sc.pp.calculate_qc_metrics(combined, percent_top=None, inplace=True)

# Preserve raw integer counts before normalization so they remain available for
# methods that require non-normalized input (e.g. negative-binomial models).
combined.layers['counts'] = combined.X.copy()

n_before = combined.n_obs
before_per_sample = combined.obs['sample'].value_counts()

sc.pp.filter_cells(combined, min_counts=MIN_COUNTS)
sc.pp.filter_cells(combined, min_genes=MIN_GENES)

# QC (04_qc_diagnostics.py) flagged AZU1 as background / off-target signal: a
# promyelocyte granule gene sitting at a flat ~1-count floor in ~89% of all cells
# (6% of total counts) — far above the negative-control-probe background, and
# implausible in plasma/erythroid/T cells. Drop it after cell filtering (so the
# per-cell QC totals still include it) but before normalization and PCA, so it
# neither inflates the size factor nor steers the components. Slicing also removes
# it from layers['counts'].
combined = combined[:, combined.var_names != 'AZU1'].copy()

n_after = combined.n_obs
n_removed = n_before - n_after
print(f'\nFiltered {n_removed:,} cells ({100 * n_removed / n_before:.1f}% removed)')
print(f'Remaining: {n_after:,} cells x {combined.n_vars} genes (AZU1 removed)')

after_per_sample = combined.obs['sample'].value_counts()
removed_per_sample = (before_per_sample - after_per_sample).rename('removed')
pct_per_sample = (removed_per_sample / before_per_sample * 100).rename('pct_removed').round(1)
print('\nRemoved per sample:')
print(pd.concat([before_per_sample.rename('before'), after_per_sample.rename('after'),
                 removed_per_sample, pct_per_sample], axis=1).to_string())

# Normalize each cell's counts to a common library size (the dataset median).
sc.pp.normalize_total(combined, inplace=True)

# Log-transform (log(1+x)) so PCA is not dominated by a few high-count genes.
# Raw integer counts remain in layers['counts'] for count-based methods.
sc.pp.log1p(combined)

if not sp.issparse(combined.X):
    combined.X = sp.csr_matrix(combined.X)

out_path = os.path.join(RUN_DIR, 'combined_filtered.h5ad')
combined.write_h5ad(out_path)
print(f'\nWrote {out_path}')
