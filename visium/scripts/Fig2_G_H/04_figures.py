#!/usr/bin/env python3
"""
Stage 04 — figures from the processed object. No analysis happens here.

    04_pca_variance_ratio   PC elbow
    04_umap_qc              UMAP by QC metrics + sample
    04_umap_leiden          UMAP by counts / genes / leiden / sample
    04_marker_genes         rank_genes_groups panel
    04_spatial_leiden_*     one per sample, coloured by cluster
    04_spatial_genes_*      one per sample, GENES side by side

Cluster ids are only valid for the 03_process.py run that produced them, so no
cell-type mapping happens here — annotate from tables/03_marker_genes.csv.

Reads:  h5ad/processed.h5ad
Writes: figures/04_*.{png,pdf}

    python visium/scripts/spatial/04_figures.py
"""

import scanpy as sc

from _common import PROCESSED_H5AD, SAMPLE_KEY, read_h5ad, save_fig

GENES = ["MZB1", "PECAM1", "GNLY"]

adata = read_h5ad(PROCESSED_H5AD)

sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
save_fig("04_pca_variance_ratio")

sc.pl.umap(adata, color=["n_genes_by_counts", "pct_counts_mt", SAMPLE_KEY],
           wspace=0.4, show=False)
save_fig("04_umap_qc")

sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "leiden", SAMPLE_KEY],
           wspace=0.4, show=False)
save_fig("04_umap_leiden")

sc.pl.rank_genes_groups(adata, n_genes=10, sharey=False, show=False)
save_fig("04_marker_genes")

for sample in adata.obs[SAMPLE_KEY].cat.categories:
    spots = adata[adata.obs[SAMPLE_KEY] == sample]
    library = str(spots.obs["library_id"].iloc[0])

    sc.pl.spatial(spots, color="leiden", library_id=library, title=sample, show=False)
    save_fig(f"04_spatial_leiden_{sample}")

    sc.pl.spatial(spots, color=GENES, library_id=library, show=False)
    save_fig(f"04_spatial_genes_{sample}")
