#!/usr/bin/env python3
"""
The one shared file for the Visium spatial pipeline: samples, paths, the
filtering rule, the processing steps, and IO helpers.

Everything else is a flat script that runs top to bottom and hands the next
stage an .h5ad:

    00_load.py      HISE zips -> outs/ -> combined            -> raw.h5ad
    01_qc.py        QC metrics + overview figure              -> qc.h5ad
    02_qc_apply.py  spot filtering on the cutoffs             -> filtered.h5ad
    03_process.py   normalize / HVG / PCA / kNN / UMAP/leiden -> processed.h5ad
    04_figures.py   figures                                   -> figures/*.png
    05_cutoff_sweep.py  diagnostic: how much does the cutoff move the UMAP?

Run them directly:  python visium/scripts/spatial/01_qc.py
Everything lands under visium/data/spatial-analysis/ (git-ignored).
"""

import os

import matplotlib
import pandas as pd

# Scripts run headless; set the backend before pyplot is imported anywhere.
matplotlib.use("Agg")

# pandas >= 3.0 infers an Arrow-backed `str` dtype for string columns/indices,
# which anndata 0.12 has no h5ad writer for. Restore object dtype for this
# process; drop this once anndata gains str-dtype support.
pd.set_option("future.infer_string", False)

import anndata as ad  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import scanpy as sc  # noqa: E402
import scipy.sparse as sp  # noqa: E402

# --- Samples ---------------------------------------------------------------
# Space Ranger output-zip UUIDs from the Visium tracking sheet.
SAMPLES = {
    "TIS05393": "f7f5aa28-4469-42cf-8db6-c75f8af37aaa",  # TIS05393-001-011, EXP-00971
    "TIS05394": "fa2b2e91-ff8d-4a48-8089-f575ec043e1e",  # TIS05394-001-010, EXP-00971
    "TIS05397": "099ef102-c537-4dd5-bef4-7a7bb0bfe316",  # TIS05397-001,     EXP-00963
    "TIS04788": "3de083f7-7e41-43b0-a3a3-a9a7af96a5db",  # TIS04788-001-004, EXP-00947
    "TIS04789": "53bfa8c0-e5c1-4227-bb32-63c3f1244d41",  # TIS04789-001-013, EXP-00948
    "TIS05392": "cc9ab020-6cab-444a-b28e-05b7643fc4af",  # TIS05392-001,     EXP-00963
}

SAMPLE_KEY = "dataset_name"
QC_METRICS = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
QC_LABELS = ["Genes per spot", "UMI counts per spot", "% mitochondrial"]

# --- Paths -----------------------------------------------------------------
REPO = "/home/workspace/PCD_data"
DATA = os.path.join(REPO, "visium", "data")
INPUT_ROOT = os.path.expanduser("~/input")  # HISE read-only FUSE cache

ANALYSIS = os.path.join(DATA, "spatial-analysis")
OUTS_BASE = os.path.join(DATA, "outs")
H5AD_DIR = os.path.join(ANALYSIS, "h5ad")
FIG_DIR = os.path.join(ANALYSIS, "figures")
TABLE_DIR = os.path.join(ANALYSIS, "tables")

RAW_H5AD = os.path.join(H5AD_DIR, "raw.h5ad")
QC_H5AD = os.path.join(H5AD_DIR, "qc.h5ad")
FILTERED_H5AD = os.path.join(H5AD_DIR, "filtered.h5ad")
PROCESSED_H5AD = os.path.join(H5AD_DIR, "processed.h5ad")

for _d in (OUTS_BASE, H5AD_DIR, FIG_DIR, TABLE_DIR):
    os.makedirs(_d, exist_ok=True)


# --- Cutoffs ---------------------------------------------------------------
# Per-sample QC cutoffs. A spot is DROPPED when it fails any of its sample's
# cuts (see dropped_mask). Edit the numbers here, re-run 01_qc.py to see the
# consequences, then 02_qc_apply.py to apply them.
#
# Set 2026-08-06 under a hard constraint: DROP AT MOST 10% OF SPOTS PER SAMPLE.
# Each sample's cut is the largest threshold whose union drop (counts OR genes
# OR mito) still fits under 10% of that sample, derived from the sample's own
# quantiles, walking observed values because counts are discrete with heavy ties
# at 1-5 UMIs. Realised drop: TIS05393 9.98%, TIS05394 9.93%, TIS04789 8.95%,
# TIS04788 8.61%, TIS05397 7.24%, TIS05392 6.12%. Overall 91.9% kept.
#
# CONSEQUENCE, on the record: the 10% cap binds well below the distribution
# structure on the four weaker captures. TIS05397 / TIS04788 / TIS04789 /
# TIS05392 are bimodal (a near-empty mode at 4-19 UMIs and a tissue mode at
# 113-992) and the valley sits at 40-131 counts — far above what a 10% budget
# allows. Cutting at 0-4 keeps the near-empty mode almost entirely: median
# counts among KEPT spots is 7 (TIS05392), 16 (TIS04788), 25 (TIS05397) against
# 302 and 170 for the two good captures. Those spots carry little signal and are
# still clustered and embedded downstream — at resolution 0.8 they fragment the
# graph into ~75 Leiden clusters. Run 05_cutoff_sweep.py to see the trade.
CUTOFFS = {
    # unimodal captures; the 10% cap is the binding constraint
    "TIS05393": {"min_genes": 14, "min_counts": 19, "max_pct_mt": 20},
    "TIS05394": {"min_genes": 22, "min_counts": 26, "max_pct_mt": 20},
    # bimodal; the cap allows far less than the valley would want
    "TIS05397": {"min_genes": 4, "min_counts": 4, "max_pct_mt": 20},   # valley at 131
    "TIS04788": {"min_genes": 2, "min_counts": 2, "max_pct_mt": 20},   # valley at 88
    "TIS04789": {"min_genes": 4, "min_counts": 4, "max_pct_mt": 20},   # valley at 58
    # no threshold above 0 fits the cap (>10% of spots sit at 1 UMI)
    "TIS05392": {"min_genes": 0, "min_counts": 0, "max_pct_mt": 20},
}


# --- Filtering -------------------------------------------------------------
def dropped_mask(obs, cutoffs):
    """
    Per-spot drop mask — the single definition of the filtering rule, used by
    the 01 preview, the 02 apply step and the 05 sweep.

    A spot is dropped when it fails ANY of its sample's cuts:
        n_genes_by_counts <= min_genes
        total_counts      <= min_counts
        pct_counts_mt     >= max_pct_mt
    """
    mask = pd.Series(False, index=obs.index)
    for sample, t in cutoffs.items():
        sel = obs[SAMPLE_KEY] == sample
        mask |= sel & (
            (obs["n_genes_by_counts"] <= t["min_genes"])
            | (obs["total_counts"] <= t["min_counts"])
            | (obs["pct_counts_mt"] >= t["max_pct_mt"])
        )
    return mask


def cutoffs_frame(cutoffs):
    """Per-sample cutoffs as a DataFrame, for printing and CSV output."""
    return pd.DataFrame(cutoffs).T[["min_genes", "min_counts", "max_pct_mt"]]


def cutoffs_for_budget(obs, budget, max_pct_mt=20.0):
    """
    Per-sample cutoffs that drop at most `budget` (a fraction) of EACH sample.
    Used only by 05_cutoff_sweep.py.

    Counts are discrete with heavy ties at 1-5 UMIs, so a plain quantile call
    overshoots the cap; this walks the observed values and takes the largest
    (min_counts, min_genes) pair whose union drop still fits. budget=0 yields
    cutoffs that drop nothing at all.
    """
    import numpy as np

    out = {}
    for sample, d in obs.groupby(SAMPLE_KEY, observed=True):
        if budget <= 0:
            # -1 is below every observed value, 101 above every percentage
            out[sample] = {"min_genes": -1, "min_counts": -1, "max_pct_mt": 101.0}
            continue

        counts = d["total_counts"].to_numpy()
        genes = d["n_genes_by_counts"].to_numpy()
        mt_fail = d["pct_counts_mt"].to_numpy() >= max_pct_mt

        best = (-1, -1)
        cand_c = np.unique(counts[counts <= np.quantile(counts, min(0.5, budget * 3))])
        cand_g = np.unique(genes[genes <= np.quantile(genes, min(0.5, budget * 3))])
        for tc in cand_c:
            for tg in cand_g[::-1]:
                if ((counts <= tc) | (genes <= tg) | mt_fail).mean() <= budget:
                    if tc > best[0] or (tc == best[0] and tg > best[1]):
                        best = (tc, tg)
                    break
        # If even the mito cut alone exceeds the budget, drop it too
        mt_cut = max_pct_mt if mt_fail.mean() <= budget else 101.0
        out[sample] = {
            "min_genes": int(best[1]),
            "min_counts": int(best[0]),
            "max_pct_mt": float(mt_cut),
        }
    return out


# --- Processing ------------------------------------------------------------
def process(adata, resolution=0.8, n_comps=50, n_pcs=20, n_neighbors=15, seed=0):
    """
    normalize -> log1p -> HVG -> PCA -> Harmony -> kNN -> UMAP -> Leiden.

    Raw counts are preserved in layers['counts'] first. Mutates and returns
    adata. Marker genes are not computed here — 03_process.py adds them.

    Kept in this shared file so the 05 sweep runs the SAME pipeline the real
    analysis runs; copies would drift and the sweep would stop being evidence.

    Harmony corrects X_pca over obs['dataset_name'] via harmonypy directly (the
    scanpy wrapper is skipped — version mismatch). X_pca is left intact and the
    corrected embedding kept in obsm['X_pca_harmony']; the graph uses the
    latter. To cluster on raw PCs instead, point use_rep at 'X_pca' below.
    """
    import harmonypy as hm

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    adata.X = sp.csr_matrix(adata.X)  # normalize_total can densify X

    sc.pp.highly_variable_genes(adata)
    sc.pp.pca(adata, n_comps=n_comps)

    ho = hm.run_harmony(adata.obsm["X_pca"].astype("float64"), adata.obs, SAMPLE_KEY)
    # Z_corr orientation varies by harmonypy version (n_pcs x n_cells or reverse)
    z = ho.Z_corr
    adata.obsm["X_pca_harmony"] = z.T if z.shape[0] != adata.n_obs else z

    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_pcs=n_pcs,
                    n_neighbors=n_neighbors, random_state=seed)
    sc.tl.umap(adata, random_state=seed)
    sc.tl.leiden(adata, resolution=resolution, flavor="igraph", n_iterations=2,
                 directed=False, random_state=seed)
    return adata


# --- IO --------------------------------------------------------------------
def read_h5ad(path):
    adata = ad.read_h5ad(path)
    adata.obs[SAMPLE_KEY] = adata.obs[SAMPLE_KEY].astype("category")
    return adata


def write_h5ad(adata, path):
    adata.write_h5ad(path)


# Every figure is written once per format: PNG to look at, PDF for slides.
FIG_FORMATS = ("png", "pdf")


def save_fig(name, dpi=300):
    """Save the current figure into FIG_DIR in each format, then close it."""
    for fmt in FIG_FORMATS:
        plt.savefig(os.path.join(FIG_DIR, f"{name}.{fmt}"), dpi=dpi, bbox_inches="tight")
    plt.close("all")


def save_table(df, name):
    df.to_csv(os.path.join(TABLE_DIR, name))
