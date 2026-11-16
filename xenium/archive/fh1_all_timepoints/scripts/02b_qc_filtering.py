import os
import warnings

import pandas as pd
import scipy.sparse as sp
import scanpy as sc

warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

DATA_DIR = '/home/workspace/PCD_data/xenium/archive/fh1_all_timepoints/data'
OBJECTS_DIR = os.path.join(DATA_DIR, 'objects')
os.makedirs(OBJECTS_DIR, exist_ok=True)

# Adjust these based on 02a_qc_calibration.py plots.
# Look for the inflection point in the knee curves.
MIN_COUNTS = 15
MIN_GENES  = 10

combined = sc.read_h5ad(os.path.join(OBJECTS_DIR, 'combined_raw.h5ad'))
print(f'Loaded: {combined.n_obs:,} cells x {combined.n_vars} genes')

sc.pp.calculate_qc_metrics(combined, percent_top=None, inplace=True)

# Preserve raw integer counts before normalization so they remain available
# for methods that require non-normalized input (e.g. negative-binomial models).
combined.layers['counts'] = combined.X.copy()

n_before = combined.n_obs
before_per_sample = combined.obs['sample'].value_counts()

sc.pp.filter_cells(combined, min_counts=MIN_COUNTS)
sc.pp.filter_cells(combined, min_genes=MIN_GENES)

n_after = combined.n_obs
n_removed = n_before - n_after
print(f'\nFiltered {n_removed:,} cells ({100 * n_removed / n_before:.1f}% removed)')
print(f'Remaining: {n_after:,} cells')

after_per_sample = combined.obs['sample'].value_counts()
removed_per_sample = (before_per_sample - after_per_sample).rename('removed')
pct_per_sample = (removed_per_sample / before_per_sample * 100).rename('pct_removed').round(1)
print('\nRemoved per sample:')
print(pd.concat([before_per_sample.rename('before'), after_per_sample.rename('after'),
                 removed_per_sample, pct_per_sample], axis=1).to_string())

# Normalize each cell's counts to a common library size (the dataset median).
# This makes transcript levels comparable across samples captured at different depths.
sc.pp.normalize_total(combined, inplace=True)

if not sp.issparse(combined.X):
    combined.X = sp.csr_matrix(combined.X)

out_path = os.path.join(OBJECTS_DIR, 'combined_filtered.h5ad')
combined.write_h5ad(out_path)
print(f'\nWrote {out_path}')
