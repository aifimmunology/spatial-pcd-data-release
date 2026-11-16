import logging
import os

import pandas as pd
import scanpy as sc

DATA_DIR    = '/home/workspace/PCD_data/xenium/commercial/data'
OBJECTS_DIR = os.path.join(DATA_DIR, 'objects')
LOGS_DIR    = os.path.join(DATA_DIR, 'logs')

for directory in (OBJECTS_DIR, LOGS_DIR):
    os.makedirs(directory, exist_ok=True)

logging.basicConfig(
    filename=os.path.join(LOGS_DIR, '02b_qc_filtering.log'),
    level=logging.INFO,
    format='%(asctime)s %(message)s',
)

# Set from the 02a knee plots.
MIN_COUNTS = 15
MIN_GENES  = 10

combined = sc.read_h5ad(os.path.join(OBJECTS_DIR, 'combined_raw.h5ad'))
# Adds per-cell obs columns total_counts and n_genes_by_counts (plus per-gene var columns),
# computed on the unfiltered raw counts. .X is not changed.
sc.pp.calculate_qc_metrics(combined, percent_top=None, inplace=True)

# Adds layers['counts']: a copy of the raw integer counts, since normalize_total overwrites .X.
combined.layers['counts'] = combined.X.copy()

before_per_sample = combined.obs['sample'].value_counts()

# Removes cells below MIN_COUNTS total counts, then cells below MIN_GENES detected genes.
# This runs before normalizing; after normalization every cell's total is the median,
# so the count threshold would no longer remove anything.
# Both calls drop whole rows (cells) in place; genes and all other columns are untouched.
sc.pp.filter_cells(combined, min_counts=MIN_COUNTS)
sc.pp.filter_cells(combined, min_genes=MIN_GENES)

after_per_sample = combined.obs['sample'].value_counts()
removed_per_sample = before_per_sample - after_per_sample
filtering_summary = pd.DataFrame({
    'before':      before_per_sample,
    'after':       after_per_sample,
    'removed':     removed_per_sample,
    'pct_removed': (100 * removed_per_sample / before_per_sample).round(1),
})
logging.info(
    'filtered with MIN_COUNTS=%d MIN_GENES=%d: %d -> %d cells\n%s',
    MIN_COUNTS, MIN_GENES, before_per_sample.sum(), after_per_sample.sum(),
    filtering_summary.to_string(),
)


out_path = os.path.join(OBJECTS_DIR, 'combined_filtered.h5ad')
# Saves the filtered object (still raw counts in .X, copy in layers['counts']) as a new file.
combined.write_h5ad(out_path)
logging.info('wrote %s', out_path)
