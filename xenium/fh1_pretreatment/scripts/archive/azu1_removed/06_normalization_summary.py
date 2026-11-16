"""
fh1_pretreatment AZU1-removed run — one-page normalization comparison summary.

Consolidates the total-count vs Pearson comparison onto a SINGLE page: the four
UMAP panels, the ARI/NMI-vs-resolution curves, the cluster-count curve, and a
marker heatmap for the Pearson clusters (the real-structure-vs-noise test — a
coherent cluster shows its own block of high-expression markers).

Read-only w.r.t. data; writes only a figure and a marker table.

    python xenium/fh1_pretreatment/scripts/archive/azu1_removed/06_normalization_summary.py
"""

import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

DATA_DIR = "/home/workspace/PCD_data/xenium/fh1_pretreatment/data/archive"
RUN_DIR = os.path.join(DATA_DIR, "azu1_removed")
OUT_DIR = "/home/workspace/PCD_data/xenium/fh1_pretreatment/figures/archive/azu1_removed/normalization_compare"
TAB_DIR = os.path.join(RUN_DIR, "normalization_compare")
os.makedirs(OUT_DIR, exist_ok=True)

KEY = "leiden_res0.60"
N_TOP = 3
PLASMA_IDENTITY = [
    "MZB1", "POU2AF1", "DERL3", "SDC1", "FCRL5", "TNFRSF17", "IRF4", "XBP1",
    "CD38", "PRDM1", "SLAMF7", "ITM2C", "CD79B", "CD79A", "TENT5C", "CD63",
]

pear = sc.read_h5ad(os.path.join(RUN_DIR, "combined_pearson_clustered.h5ad"))
tot = sc.read_h5ad(os.path.join(RUN_DIR, "combined_clustered.h5ad"))
sim = pd.read_csv(os.path.join(TAB_DIR, "normalization_similarity.csv"))

# Log-norm expression from raw counts — the interpretable space for markers and
# for scoring, regardless of which normalization drove each embedding.
ln = sc.AnnData(pear.layers["counts"].copy(), obs=pear.obs.copy(), var=pear.var.copy())
sc.pp.normalize_total(ln)
sc.pp.log1p(ln)
ident = [g for g in PLASMA_IDENTITY if g in ln.var_names]
sc.tl.score_genes(ln, ident, score_name="plasma_score")
plasma_score = ln.obs["plasma_score"].to_numpy()
tot.obs["plasma_score"] = plasma_score
pear.obs["plasma_score"] = plasma_score

# Markers for the Pearson clusters, on log-norm expression.
ln.obs[KEY] = pear.obs[KEY].values
sc.tl.rank_genes_groups(ln, KEY, method="wilcoxon")
groups = list(ln.obs[KEY].cat.categories)
top = {g: [ln.uns["rank_genes_groups"]["names"][g][i] for i in range(N_TOP)] for g in groups}

print("Top markers per Pearson cluster (res 0.6):")
for g in groups:
    print(f"  cluster {g:>2}: {', '.join(top[g])}")

marker_order = list(dict.fromkeys(gene for g in groups for gene in top[g]))
sub = ln[:, marker_order]
expr = pd.DataFrame(sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X),
                    columns=marker_order)
expr[KEY] = ln.obs[KEY].values
mean_by_cl = expr.groupby(KEY, observed=True)[marker_order].mean().reindex(groups)
z = (mean_by_cl - mean_by_cl.mean(0)) / (mean_by_cl.std(0) + 1e-9)
mean_by_cl.to_csv(os.path.join(TAB_DIR, "pearson_cluster_marker_means.csv"))

# ── One-page layout ────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 19))
gs = fig.add_gridspec(4, 2, height_ratios=[1.0, 1.0, 0.7, 1.5], hspace=0.35, wspace=0.2,
                      top=0.95, bottom=0.03)

ax = fig.add_subplot(gs[0, 0])
sc.pl.umap(tot, color=KEY, ax=ax, show=False, legend_loc="on data", size=5,
           title="total-count — clusters (res 0.6)")
ax = fig.add_subplot(gs[0, 1])
sc.pl.umap(pear, color=KEY, ax=ax, show=False, legend_loc="on data", size=5,
           title="Pearson — clusters (res 0.6)")
ax = fig.add_subplot(gs[1, 0])
sc.pl.umap(tot, color="plasma_score", ax=ax, show=False, size=5,
           title="total-count — plasma score")
ax = fig.add_subplot(gs[1, 1])
sc.pl.umap(pear, color="plasma_score", ax=ax, show=False, size=5,
           title="Pearson — plasma score")

ax = fig.add_subplot(gs[2, 0])
ax.plot(sim["resolution"], sim["ari_pearson_vs_lognorm"], "-o", label="Pearson vs total-count")
ax.plot(sim["resolution"], sim["ari_noise_floor"], "-o", label="noise floor (seed)")
ax.plot(sim["resolution"], sim["nmi_pearson_vs_lognorm"], "--s", label="NMI Pearson vs total")
ax.set_ylim(0, 1.02); ax.set_xlabel("Leiden resolution"); ax.set_ylabel("cluster agreement")
ax.set_title("agreement vs noise floor"); ax.legend(fontsize=8)

ax = fig.add_subplot(gs[2, 1])
ax.plot(sim["resolution"], sim["n_clusters_lognorm"], "-o", label="total-count")
ax.plot(sim["resolution"], sim["n_clusters_pearson"], "-o", label="Pearson")
ax.set_xlabel("Leiden resolution"); ax.set_ylabel("# clusters")
ax.set_title("cluster count vs resolution"); ax.legend(fontsize=8)

ax = fig.add_subplot(gs[3, :])
im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)
ax.set_xticks(range(len(marker_order)))
ax.set_xticklabels(marker_order, rotation=90, fontsize=6)
ax.set_yticks(range(len(groups)))
ax.set_yticklabels([f"c{g}" for g in groups], fontsize=7)
ax.set_title("Pearson clusters × top markers (z-scored mean log-norm) — a clean diagonal = coherent clusters")
fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01, label="z-score")

fig.suptitle("fh1_pretreatment (AZU1-removed): total-count vs Pearson normalization", fontsize=14, y=0.98)

# Rasterize the scatter points (the UMAP PathCollections) so the PDF is a light
# embedded image, not ~1.6M vector circles; axes, text, curves and the heatmap
# stay vector. Post-hoc set_rasterized per collection is the robust approach —
# sc.settings._vector_friendly and the rasterized= kwarg are version-brittle.
# See xenium/commercial/scripts/05_umap_by_sample.py.
for ax in fig.axes:
    for coll in ax.collections:
        coll.set_rasterized(True)

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT_DIR, f"normalization_summary.{ext}"), dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"\nWrote {os.path.join(OUT_DIR, 'normalization_summary.pdf')}")
