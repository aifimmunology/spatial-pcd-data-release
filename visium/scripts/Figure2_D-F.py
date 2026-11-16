import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq
import hisepy as hp
import glob, os, zipfile
import anndata as ad

INPUT_ROOT = os.path.expanduser("~/input")
OUTS_BASE  = "/home/workspace/PCD_data/visium/data/outs"
os.makedirs(OUTS_BASE, exist_ok=True)

all_zips = glob.glob(os.path.join(INPUT_ROOT, "**", "*.zip"), recursive=True)

adatas = {}
for sample_id, uuid in SAMPLES.items():
    sample_outs = os.path.join(OUTS_BASE, sample_id)
    adata = sq.read.visium(sample_outs)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    adatas[sample_id] = adata

combined_adatas = sc.concat(adatas, label='dataset_name', join='outer')
combined_adatas.obs_names_make_unique()

combined_adatas.uns['spatial'] = {}
for sample_id, adata in adatas.items():
    combined_adatas.uns['spatial'].update(adata.uns['spatial'])


combined_adatas.var['mt'] = combined_adatas.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(
    combined_adatas, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True
)


# --- Filter spots ---
# Adjust thresholds after reviewing the violin/scatter plots above
MIN_GENES  = 20
MIN_COUNTS = 50
MAX_MT_PCT = 20

combined_adatas = combined_adatas[combined_adatas.obs['n_genes_by_counts']   > MIN_GENES ].copy()
combined_adatas = combined_adatas[combined_adatas.obs['total_counts']        > MIN_COUNTS].copy()
combined_adatas = combined_adatas[combined_adatas.obs['pct_counts_mt']       < MAX_MT_PCT].copy()

combined_adatas.layers["counts"] = combined_adatas.X.copy()
sc.pp.normalize_total(combined_adatas, inplace=True)
sc.pp.log1p(combined_adatas)

# Ensure X stays sparse after normalization
if not sp.issparse(combined_adatas.X):
    combined_adatas.X = sp.csr_matrix(combined_adatas.X)
    
sc.pp.highly_variable_genes(combined_adatas)
sc.pp.pca(combined_adatas)
sc.pp.neighbors(combined_adatas)
sc.tl.umap(combined_adatas)
sc.tl.leiden(combined_adatas)