#!/usr/bin/env python3
"""
Stage 03 — normalize -> log1p -> HVG -> PCA -> Harmony -> kNN -> UMAP -> Leiden,
then marker genes per cluster.

The pipeline itself is _common.process(), shared with 05_cutoff_sweep.py so the
sweep and the real analysis cannot drift apart. Raw counts are kept in
layers['counts'].

Reads:  h5ad/filtered.h5ad
Writes: h5ad/processed.h5ad, tables/03_marker_genes.csv

    python visium/scripts/spatial/03_process.py
"""

import scanpy as sc

from _common import (
    FILTERED_H5AD,
    PROCESSED_H5AD,
    process,
    read_h5ad,
    save_table,
    write_h5ad,
)

RESOLUTION = 0.8

adata = read_h5ad(FILTERED_H5AD)
process(adata, resolution=RESOLUTION)

print(f"{adata.obs['leiden'].nunique()} Leiden clusters at resolution {RESOLUTION}")
print(adata.obs["leiden"].value_counts().sort_index().to_string())

sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
save_table(sc.get.rank_genes_groups_df(adata, group=None), "03_marker_genes.csv")

write_h5ad(adata, PROCESSED_H5AD)
