import os
import warnings

import pandas as pd
import scipy.sparse as sp
import scanpy as sc

warnings.filterwarnings('ignore')
sc.settings.verbosity = 1

DATA_DIR    = '/home/workspace/PCD_data/xenium/fh1_pretreatment/data'
OBJECTS_DIR = os.path.join(DATA_DIR, 'objects')
os.makedirs(OBJECTS_DIR, exist_ok=True)

# Adjust these based on the 02a_qc_calibration.py knee plots if the inflection
# point for this sample set suggests different cutoffs.
MIN_COUNTS = 15
MIN_GENES  = 10

# AZU1 behaved like background / off-target signal rather than real expression
# (the evidence is in figures/archive/lognorm_azu1_kept/qc_diagnostics/), so it is
# dropped before normalization — its counts must not inflate any cell's library size.
DROP_GENES = ['AZU1']

combined = sc.read_h5ad(os.path.join(OBJECTS_DIR, 'combined_raw.h5ad'))
print(f'Loaded: {combined.n_obs:,} cells x {combined.n_vars} genes')

# Adds obs total_counts and n_genes_by_counts (plus per-gene var columns), computed on the
# unfiltered raw counts. .X is not changed.
sc.pp.calculate_qc_metrics(combined, percent_top=None, inplace=True)

# Adds layers['counts'] (copy of the raw counts). Preserve raw integer counts before normalization so they remain available for
# methods that require non-normalized input (e.g. negative-binomial models, and
# the Pearson residuals in 03_processing.py).
combined.layers['counts'] = combined.X.copy()

n_before = combined.n_obs
before_per_sample = combined.obs['sample'].value_counts()

# Both calls drop whole rows (cells) in place, on raw counts, before normalization.
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

# Drops the DROP_GENES columns from .X, layers['counts'] and var (the only gene removal in the
# pipeline); total_counts and n_genes_by_counts in obs were computed earlier and still include them.
combined = combined[:, [g for g in combined.var_names if g not in set(DROP_GENES)]].copy()
print(f'\nDropped {DROP_GENES}: {combined.n_vars} genes remain')

# Overwrites .X: rescales each cell's counts to a common library size (the dataset median) so
# transcript levels are comparable across samples captured at different depths.
sc.pp.normalize_total(combined, inplace=True)

# Overwrites .X again with log(1+x) of the normalized values, so PCA is not dominated by a few high-count genes and
# the finer distinctions between related cell types survive into the embedding.
# Raw integer counts remain in layers['counts'] for count-based methods.
sc.pp.log1p(combined)

# Changes .X's storage format only (dense to sparse), not its values.
if not sp.issparse(combined.X):
    combined.X = sp.csr_matrix(combined.X)

# Saves the filtered object as a new file: .X is log-normalized, layers['counts'] is raw.
out_path = os.path.join(OBJECTS_DIR, 'combined_filtered.h5ad')
combined.write(out_path)
print(f'\nWrote {out_path}')
