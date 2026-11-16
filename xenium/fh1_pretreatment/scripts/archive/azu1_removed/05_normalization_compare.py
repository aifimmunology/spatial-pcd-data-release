"""
fh1_pretreatment AZU1-removed run — normalization comparison: the total-count pipeline
(normalize_total -> log1p -> PCA) vs a composition-robust path (analytic Pearson
residuals -> PCA), the scanpy-native form of SCTransform. Everything downstream
(PCA dims, Harmony, neighbours, Leiden sweep) is identical, so any difference is
the normalization.

Pearson residuals are computed from RAW counts and fed straight to PCA: they
replace normalize_total, log1p AND sc.pp.scale (a residual is already
depth-corrected and variance-stabilized per gene).

All inputs and outputs live under the objects_azu1_removed/ run; the original
objects/ folder is never touched.

    python xenium/fh1_pretreatment/scripts/archive/azu1_removed/05_normalization_compare.py
"""

import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy as hm
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

DATA_DIR = "/home/workspace/PCD_data/xenium/fh1_pretreatment/data/archive"
RUN_DIR = os.path.join(DATA_DIR, "azu1_removed")
FILTERED = os.path.join(RUN_DIR, "combined_filtered.h5ad")
PEARSON_OUT = os.path.join(RUN_DIR, "combined_pearson_clustered.h5ad")
FIG_DIR = "/home/workspace/PCD_data/xenium/fh1_pretreatment/figures/archive/azu1_removed/normalization_compare"
TAB_DIR = os.path.join(RUN_DIR, "normalization_compare")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TAB_DIR, exist_ok=True)

BATCH_KEY = "cohort_subject"
N_COMPS = 30
N_NEIGHBORS = 15
RESOLUTIONS = [round(0.1 * i, 2) for i in range(1, 12)]
REF_RES = 0.6
PLASMA_IDENTITY = [
    "MZB1", "POU2AF1", "DERL3", "SDC1", "FCRL5", "TNFRSF17", "IRF4", "XBP1",
    "CD38", "PRDM1", "SLAMF7", "ITM2C", "CD79B", "CD79A", "TENT5C", "CD63",
]


def build(counts_ad, method, seed):
    ad = counts_ad.copy()
    if method == "lognorm":
        sc.pp.normalize_total(ad)
        sc.pp.log1p(ad)
    else:
        sc.experimental.pp.normalize_pearson_residuals(ad)
    sc.pp.pca(ad, n_comps=N_COMPS, use_highly_variable=False, random_state=seed)
    ho = hm.run_harmony(ad.obsm["X_pca"].astype("float64"), ad.obs, BATCH_KEY,
                        random_state=seed)
    z = ho.Z_corr
    ad.obsm["X_pca_harmony"] = z.T if z.shape[0] != ad.n_obs else z
    sc.pp.neighbors(ad, use_rep="X_pca_harmony", n_neighbors=N_NEIGHBORS,
                    n_pcs=N_COMPS, random_state=seed)
    sc.tl.umap(ad, random_state=seed)
    labels = {}
    for res in RESOLUTIONS:
        sc.tl.leiden(ad, resolution=res, flavor="igraph", n_iterations=2,
                     directed=False, random_state=seed, key_added=f"leiden_res{res:.2f}")
        labels[res] = ad.obs[f"leiden_res{res:.2f}"].to_numpy().copy()
    return ad, labels


def main():
    print("=== fh1_pretreatment (AZU1-removed) normalization comparison: total-count vs Pearson ===")
    filtered = sc.read_h5ad(FILTERED)
    counts_ad = sc.AnnData(filtered.layers["counts"].copy(),
                           obs=filtered.obs.copy(), var=filtered.var.copy())
    print(f"Loaded: {counts_ad.n_obs:,} cells x {counts_ad.n_vars} genes")

    # Plasma score, computed once from a log-norm view so both embeddings are
    # coloured by the same identity signal.
    ln = counts_ad.copy()
    sc.pp.normalize_total(ln)
    sc.pp.log1p(ln)
    ident = [g for g in PLASMA_IDENTITY if g in ln.var_names]
    sc.tl.score_genes(ln, ident, score_name="plasma_score")
    plasma_score = ln.obs["plasma_score"].to_numpy()

    print("Building total-count path (seed 0)…")
    ad_log, lab_log = build(counts_ad, "lognorm", 0)
    print("Building total-count path (seed 1, noise floor)…")
    _, lab_log2 = build(counts_ad, "lognorm", 1)
    print("Building Pearson-residual path (seed 0)…")
    ad_pr, lab_pr = build(counts_ad, "pearson", 0)

    rows = []
    for res in RESOLUTIONS:
        rows.append({
            "resolution": res,
            "ari_pearson_vs_lognorm": adjusted_rand_score(lab_log[res], lab_pr[res]),
            "nmi_pearson_vs_lognorm": normalized_mutual_info_score(lab_log[res], lab_pr[res]),
            "ari_noise_floor": adjusted_rand_score(lab_log[res], lab_log2[res]),
            "nmi_noise_floor": normalized_mutual_info_score(lab_log[res], lab_log2[res]),
            "n_clusters_lognorm": len(np.unique(lab_log[res])),
            "n_clusters_pearson": len(np.unique(lab_pr[res])),
        })
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(TAB_DIR, "normalization_similarity.csv"), index=False)
    print(df.round(3).to_string(index=False))

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, metric in zip(axes, ["ari", "nmi"]):
        ax.plot(df["resolution"], df[f"{metric}_pearson_vs_lognorm"], "-o",
                label="Pearson vs total-count")
        ax.plot(df["resolution"], df[f"{metric}_noise_floor"], "-o",
                label="noise floor (seed)")
        ax.set_xlabel("Leiden resolution"); ax.set_ylabel(metric.upper())
        ax.set_ylim(0, 1.02); ax.legend()
    axes[0].set_title("Does Pearson normalization change clustering?")
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"similarity_vs_resolution.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Side-by-side UMAPs at the reference resolution: clusters, then plasma score.
    ad_log.obs["plasma_score"] = plasma_score
    ad_pr.obs["plasma_score"] = plasma_score
    key = f"leiden_res{REF_RES:.2f}"
    fig, axes = plt.subplots(2, 2, figsize=(12, 11))
    sc.pl.umap(ad_log, color=key, ax=axes[0, 0], show=False, legend_loc="on data",
               size=6, title=f"total-count — clusters (res {REF_RES})")
    sc.pl.umap(ad_pr, color=key, ax=axes[0, 1], show=False, legend_loc="on data",
               size=6, title=f"Pearson — clusters (res {REF_RES})")
    sc.pl.umap(ad_log, color="plasma_score", ax=axes[1, 0], show=False, size=6,
               title="total-count — plasma score")
    sc.pl.umap(ad_pr, color="plasma_score", ax=axes[1, 1], show=False, size=6,
               title="Pearson — plasma score")
    fig.tight_layout()
    # Rasterize the UMAP scatter points so the PDF stays light (see 06 / commercial 05_umap_by_sample.py).
    for a in fig.axes:
        for coll in a.collections:
            coll.set_rasterized(True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG_DIR, f"umap_side_by_side_res{REF_RES}.{ext}"),
                    dpi=150, bbox_inches="tight")
    plt.close(fig)

    # Persist the Pearson-normalized, clustered object for downstream use. Raw
    # counts are kept in a layer; X holds the residuals.
    ad_pr.layers["counts"] = counts_ad.X.copy()
    ad_pr.write_h5ad(PEARSON_OUT)
    print(f"\nWrote {PEARSON_OUT}")
    print(f"Outputs → {FIG_DIR}  and  {TAB_DIR}")


if __name__ == "__main__":
    main()
