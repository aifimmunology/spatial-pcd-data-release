"""
QC diagnostics for the fh1_pretreatment Xenium pipeline. Standalone and read-only with
respect to the canonical pipeline: it reads combined_filtered.h5ad and
combined_clustered.h5ad, and writes only into qc_diagnostics/ output folders.
Nothing here modifies the pipeline objects.

Three independent tests, one section each:
  Test 1  Background / ambient-probe gene screen. A triage screen that flags
          abundant + ubiquitous + flat genes (AZU1 is the lead suspect), anchored
          to the negative-control-probe background. It flags candidates; it does
          not, on its own, prove a gene is ambient vs off-target vs segmentation
          spillover (that needs a spatial analysis, deferred).
  Test 2  Composition-bias check on plasma cells. Total-count normalization is
          compared against two references that do NOT depend on the gene counts
          the bias is suspected in: the cell-area size factor (immune to gene
          composition by construction) and analytic Pearson residuals.
  Test 3  Sensitivity of clustering to sc.pp.scale before PCA, measured across the
          resolution sweep against a same-method-different-seed noise floor.

    python xenium/fh1_pretreatment/scripts/archive/lognorm_azu1_kept/04_qc_diagnostics.py
"""

import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scanpy as sc
import harmonypy as hm
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

# ── Config (the ONLY part that differs between the commercial and fh1_pretreatment copies) ───────
CONFIG = {
    "analysis": "fh1_pretreatment",
    "data_dir": "/home/workspace/PCD_data/xenium/fh1_pretreatment/data/archive/lognorm_azu1_kept",
    "fig_dir": "/home/workspace/PCD_data/xenium/fh1_pretreatment/figures/archive/lognorm_azu1_kept/qc_diagnostics",
    "batch_key": "cohort_subject",
    "cluster_key": "leiden_res0.60",
    "n_comps": 30,
    "n_neighbors": 15,
    "resolutions": [round(0.1 * i, 2) for i in range(1, 12)],  # 0.1 … 1.1
    # Plasma-cell program, validated against this panel (identity + secretory/ER).
    # IDENTITY is also used to score/identify plasma cells; the full program is
    # used only for the descriptive "share of counts" plot — NOT to build the
    # size factor (that would be circular).
    "plasma_identity": [
        "MZB1", "POU2AF1", "DERL3", "SDC1", "FCRL5", "TNFRSF17", "IRF4", "XBP1",
        "CD38", "PRDM1", "SLAMF7", "ITM2C", "CD79B", "CD79A", "TENT5C", "CD63",
    ],
    "plasma_secretory": [
        "TXNDC5", "SEC11C", "FKBP11", "MYDGF", "MANF", "SDF2L1", "PRDX4", "PDIA6",
        "HYOU1", "SERP1", "CRELD2", "P4HB", "PRDX5", "EMC7", "NME1",
    ],
}

OBJECTS_DIR = os.path.join(CONFIG["data_dir"], "objects")
FILTERED = os.path.join(OBJECTS_DIR, "combined_filtered.h5ad")
CLUSTERED = os.path.join(OBJECTS_DIR, "combined_clustered.h5ad")
TAB_DIR = os.path.join(CONFIG["data_dir"], "qc_diagnostics")
os.makedirs(TAB_DIR, exist_ok=True)

PLASMA_PROGRAM = CONFIG["plasma_identity"] + CONFIG["plasma_secretory"]


# ── Small IO / helper utilities ───────────────────────────────────────────────
def present(genes, adata):
    return [g for g in genes if g in adata.var_names]


def savefig(fig, section, name):
    d = os.path.join(CONFIG["fig_dir"], section)
    os.makedirs(d, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(d, f"{name}.{ext}"), dpi=150, bbox_inches="tight")
    plt.close(fig)


def savecsv(df, name, index=True):
    df.to_csv(os.path.join(TAB_DIR, name), index=index)


def lognorm_from_counts(clustered):
    """A fresh log1p(normalize_total(raw counts)) view, so the two figures'
    differing filtered-X conventions don't affect scoring or Test 2."""
    ln = sc.AnnData(clustered.layers["counts"].copy(),
                    obs=clustered.obs.copy(), var=clustered.var.copy())
    sc.pp.normalize_total(ln)
    sc.pp.log1p(ln)
    return ln


def plasma_split(clustered, plasma_score):
    """Binary plasma / non-plasma split, chosen by the largest natural gap in the
    per-cluster mean plasma score. Reproducible and cutoff-free; the chosen
    clusters are reported so the split is auditable."""
    cl = clustered.obs[CONFIG["cluster_key"]]
    by_cl = pd.Series(plasma_score, index=cl.index).groupby(cl, observed=True).mean()
    by_cl = by_cl.sort_values(ascending=False)
    gaps = by_cl.values[:-1] - by_cl.values[1:]
    k = int(np.argmax(gaps)) + 1
    plasma_clusters = list(by_cl.index[:k])
    mask = cl.isin(plasma_clusters).to_numpy()
    return mask, plasma_clusters, by_cl


# ── Test 1 — background / ambient-probe screen ────────────────────────────────
def test1_background(clustered):
    print("Test 1 — background / ambient screen")
    counts = clustered.layers["counts"].tocsr()
    genes = np.asarray(clustered.var_names)
    cl = clustered.obs[CONFIG["cluster_key"]]
    clusters = list(cl.cat.categories) if hasattr(cl, "cat") else sorted(cl.unique())

    per_gene_total = np.asarray(counts.sum(0)).ravel()
    abundance = per_gene_total / per_gene_total.sum()

    # Per-cluster mean expression and detection rate.
    binar = (counts > 0).astype(np.float32)
    mean_by_cl = np.zeros((len(genes), len(clusters)))
    det_by_cl = np.zeros((len(genes), len(clusters)))
    codes = cl.cat.codes.to_numpy() if hasattr(cl, "cat") else \
        pd.Categorical(cl, categories=clusters).codes
    for j, _ in enumerate(clusters):
        rows = np.where(codes == j)[0]
        if rows.size == 0:
            continue
        mean_by_cl[:, j] = np.asarray(counts[rows].mean(0)).ravel()
        det_by_cl[:, j] = np.asarray(binar[rows].mean(0)).ravel()

    top_share = mean_by_cl.max(1) / (mean_by_cl.sum(1) + 1e-12)
    flatness = 1 - top_share                       # 1 = spread evenly across clusters
    ubiquity = (det_by_cl > 0.5).mean(1)           # fraction of clusters where detected in >50% cells
    suspect = abundance * flatness * ubiquity      # abundant AND flat AND ubiquitous

    table = pd.DataFrame({
        "gene": genes,
        "pct_of_all_counts": abundance * 100,
        "flatness": flatness,
        "ubiquity": ubiquity,
        "top_cluster": [clusters[i] for i in mean_by_cl.argmax(1)],
        "top_cluster_mean": mean_by_cl.max(1),
        "min_cluster_mean": mean_by_cl.min(1),
        "suspect_score": suspect,
    }).sort_values("suspect_score", ascending=False)
    savecsv(table, "test1_background_screen.csv", index=False)

    # Control-probe background anchor: how big is a real gene's floor vs pure
    # ambient (the negative-control-probe rate)?
    total_counts = np.asarray(counts.sum(1)).ravel()
    ctrl = clustered.obs.get("control_probe_counts")
    ctrl_rate = None
    if ctrl is not None:
        ctrl = ctrl.to_numpy().astype(float)
        ctrl_rate = np.divide(ctrl.sum(), total_counts.sum())  # panel-wide control fraction

    # Scatter: flatness vs abundance, colored by ubiquity, top suspects labeled.
    fig, ax = plt.subplots(figsize=(6, 5))
    sctr = ax.scatter(flatness, abundance * 100, c=ubiquity, s=8, cmap="viridis",
                      vmin=0, vmax=1, rasterized=True)
    ax.set_yscale("log")
    ax.set_xlabel("flatness (1 = evenly spread across clusters)")
    ax.set_ylabel("% of all counts (log)")
    ax.set_title(f"{CONFIG['analysis']}: background-gene screen")
    for _, r in table.head(6).iterrows():
        ax.annotate(r["gene"], (r["flatness"], r["pct_of_all_counts"]), fontsize=8)
    fig.colorbar(sctr, ax=ax, label="ubiquity")
    savefig(fig, "background", "screen_scatter")

    # Characterize the top suspect.
    suspect_gene = table.iloc[0]["gene"]
    ln = lognorm_from_counts(clustered)
    gi = list(genes).index(suspect_gene)

    row_floor = table.iloc[0]["min_cluster_mean"]
    with open(os.path.join(TAB_DIR, "test1_top_suspect_note.txt"), "w") as fh:
        fh.write(f"Top background suspect: {suspect_gene}\n")
        fh.write(f"  % of all counts: {table.iloc[0]['pct_of_all_counts']:.2f}\n")
        fh.write(f"  flatness: {table.iloc[0]['flatness']:.3f}  "
                 f"ubiquity: {table.iloc[0]['ubiquity']:.3f}\n")
        fh.write(f"  lowest-cluster mean count (the 'floor'): {row_floor:.3f}\n")
        if ctrl_rate is not None:
            fh.write(f"  panel-wide negative-control-probe fraction of counts: "
                     f"{ctrl_rate*100:.4f}%\n")
            fh.write("  A floor far above the control fraction argues the signal is "
                     "NOT generic ambient (off-target or segmentation spillover more "
                     "likely). This screen does not decide which; a spatial spillover "
                     "analysis is the deferred follow-up.\n")

    # Per-cluster mean + detection bars for the suspect.
    order = np.argsort(-mean_by_cl[gi])
    fig, axes = plt.subplots(2, 1, figsize=(max(6, len(clusters) * 0.4), 6), sharex=True)
    axes[0].bar(range(len(clusters)), mean_by_cl[gi][order])
    axes[0].set_ylabel("mean count")
    axes[0].set_title(f"{suspect_gene}: per-cluster expression ({CONFIG['cluster_key']})")
    axes[1].bar(range(len(clusters)), det_by_cl[gi][order] * 100)
    axes[1].set_ylabel("% cells detected")
    axes[1].set_xticks(range(len(clusters)))
    axes[1].set_xticklabels([clusters[i] for i in order], rotation=90, fontsize=7)
    axes[1].set_xlabel("cluster")
    savefig(fig, "background", f"{suspect_gene}_per_cluster")

    # UMAP + one-sample spatial map of the suspect.
    clustered.obs["_suspect_expr"] = np.asarray(
        ln[:, suspect_gene].X.todense()).ravel() if sp.issparse(ln.X) \
        else np.asarray(ln[:, suspect_gene].X).ravel()
    fig = sc.pl.umap(clustered, color="_suspect_expr", size=6, show=False,
                     return_fig=True, title=f"{suspect_gene} (log-norm)")
    savefig(fig, "background", f"{suspect_gene}_umap")

    if "spatial" in clustered.obsm and "sample" in clustered.obs:
        s0 = clustered.obs["sample"].value_counts().index[0]
        m = (clustered.obs["sample"] == s0).to_numpy()
        xy = clustered.obsm["spatial"][m]
        fig, ax = plt.subplots(figsize=(6, 6))
        sctr = ax.scatter(xy[:, 0], xy[:, 1], c=clustered.obs["_suspect_expr"].to_numpy()[m],
                          s=2, cmap="magma", rasterized=True)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_title(f"{suspect_gene} — {s0}")
        fig.colorbar(sctr, ax=ax, label="log-norm")
        savefig(fig, "background", f"{suspect_gene}_spatial_{s0}")

    del clustered.obs["_suspect_expr"]
    print(f"  top suspect: {suspect_gene} "
          f"({table.iloc[0]['pct_of_all_counts']:.1f}% of counts, "
          f"flatness {table.iloc[0]['flatness']:.2f})")


# ── Test 2 — composition bias on plasma cells ─────────────────────────────────
def test2_composition(clustered):
    print("Test 2 — composition bias (area- and Pearson-referenced)")
    counts = clustered.layers["counts"].tocsr()
    genes = np.asarray(clustered.var_names)
    total = np.asarray(counts.sum(1)).ravel().astype(float)

    ln = lognorm_from_counts(clustered)
    ident = present(CONFIG["plasma_identity"], ln)
    sc.tl.score_genes(ln, ident, score_name="plasma_score")
    plasma_score = ln.obs["plasma_score"].to_numpy()
    mask, plasma_clusters, by_cl = plasma_split(clustered, plasma_score)
    print(f"  plasma clusters ({CONFIG['cluster_key']}): {plasma_clusters}  "
          f"({mask.sum():,} cells)")

    # Size factors: total-count (suspected) vs cell-area (independent of counts).
    sf_total = total / np.median(total)
    area = clustered.obs.get("cell_area")
    if area is None:
        print("  cell_area missing — skipping area-referenced parts of Test 2")
        return
    area = area.to_numpy().astype(float)
    sf_area = area / np.median(area)
    ratio = sf_total / sf_area   # >1 in a cell = its counts are inflated per unit area

    prog = present(PLASMA_PROGRAM, clustered)
    prog_idx = [list(genes).index(g) for g in prog]
    prog_frac = np.asarray(counts[:, prog_idx].sum(1)).ravel() / np.where(total > 0, total, 1)

    # 2a — plasma-program share of counts, plasma vs rest.
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    ax.violinplot([prog_frac[mask] * 100, prog_frac[~mask] * 100], showmedians=True)
    ax.set_xticks([1, 2]); ax.set_xticklabels(["plasma", "non-plasma"])
    ax.set_ylabel("plasma-program % of a cell's counts")
    ax.set_title(f"{CONFIG['analysis']}: program share")
    savefig(fig, "composition", "2a_program_fraction")

    # 2b — size-factor ratio vs plasma score (cutoff-free), subsampled for plotting.
    rng = np.random.default_rng(0)
    take = rng.choice(len(ratio), size=min(20000, len(ratio)), replace=False)
    fig, ax = plt.subplots(figsize=(5, 4.5))
    ax.scatter(plasma_score[take], ratio[take], s=4, alpha=0.3, rasterized=True)
    ax.axhline(1.0, color="k", lw=0.8, ls="--")
    ax.set_xlabel("plasma score"); ax.set_ylabel("size-factor ratio  total / area")
    ax.set_title("count inflation vs plasma identity")
    savefig(fig, "composition", "2b_ratio_vs_score")

    # 2c — size-factor ratio, plasma vs rest.
    fig, ax = plt.subplots(figsize=(4.5, 4.5))
    ax.violinplot([ratio[mask], ratio[~mask]], showmedians=True)
    ax.axhline(1.0, color="k", lw=0.8, ls="--")
    ax.set_xticks([1, 2]); ax.set_xticklabels(["plasma", "non-plasma"])
    ax.set_ylabel("size-factor ratio  total / area")
    savefig(fig, "composition", "2c_ratio_split")

    # 2d — per non-program gene, plasma-cell deflation under total- vs area-norm.
    nonprog_idx = np.array([i for i in range(len(genes)) if i not in set(prog_idx)])

    def gene_means(idx_cells, sf):
        c = counts[idx_cells]
        inv = sp.diags(1.0 / sf[idx_cells])
        return np.asarray((inv @ c).mean(0)).ravel()

    p_idx = np.where(mask)[0]
    n_idx = np.where(~mask)[0]
    defl_plasma = gene_means(p_idx, sf_total)[nonprog_idx] / \
        (gene_means(p_idx, sf_area)[nonprog_idx] + 1e-12)
    defl_non = gene_means(n_idx, sf_total)[nonprog_idx] / \
        (gene_means(n_idx, sf_area)[nonprog_idx] + 1e-12)

    fig, ax = plt.subplots(figsize=(5, 4.5))
    bins = np.linspace(0.4, 1.6, 40)
    ax.hist(defl_plasma, bins=bins, alpha=0.6, label="plasma", density=True)
    ax.hist(defl_non, bins=bins, alpha=0.6, label="non-plasma", density=True)
    ax.axvline(1.0, color="k", lw=0.8, ls="--")
    ax.set_xlabel("total-norm / area-norm  (non-program genes)")
    ax.set_ylabel("density"); ax.legend()
    ax.set_title("plasma-cell deflation of non-program genes")
    savefig(fig, "composition", "2d_gene_deflation")

    # Pearson-residual cross-check: a second, independent method. For each
    # non-program gene, the standardized plasma-cell mean under log-norm vs under
    # Pearson residuals. Composition bias pushes the log-norm value negative; if
    # Pearson (which does not use the area factor) is less negative, two
    # independent methods agree that total-count log-norm deflates plasma cells.
    pr = sc.AnnData(counts.copy(), obs=clustered.obs.copy(), var=clustered.var.copy())
    sc.experimental.pp.normalize_pearson_residuals(pr)
    prX = np.asarray(pr.X.todense()) if sp.issparse(pr.X) else np.asarray(pr.X)
    lnX = np.asarray(ln.X.todense()) if sp.issparse(ln.X) else np.asarray(ln.X)

    def zbar_plasma(mat):
        mu = mat.mean(0); sd = mat.std(0) + 1e-9
        return (mat[p_idx].mean(0) - mu) / sd

    zb_ln = zbar_plasma(lnX)[nonprog_idx]
    zb_pr = zbar_plasma(prX)[nonprog_idx]
    fig, ax = plt.subplots(figsize=(5, 4.5))
    b = np.linspace(-1.0, 1.0, 40)
    ax.hist(zb_ln, bins=b, alpha=0.6, label="total-count log-norm", density=True)
    ax.hist(zb_pr, bins=b, alpha=0.6, label="Pearson residuals", density=True)
    ax.axvline(0.0, color="k", lw=0.8, ls="--")
    ax.set_xlabel("standardized plasma-cell mean (non-program genes)")
    ax.set_ylabel("density"); ax.legend()
    ax.set_title("plasma deflation: log-norm vs Pearson")
    savefig(fig, "composition", "2e_pearson_crosscheck")

    summary = pd.DataFrame({
        "metric": [
            "plasma cells", "non-plasma cells",
            "median program % of counts (plasma)", "median program % of counts (non)",
            "median size-factor ratio total/area (plasma)",
            "median size-factor ratio total/area (non)",
            "median non-program deflation total/area (plasma)",
            "median non-program deflation total/area (non)",
            "median standardized plasma mean, log-norm",
            "median standardized plasma mean, Pearson",
        ],
        "value": [
            int(mask.sum()), int((~mask).sum()),
            np.median(prog_frac[mask]) * 100, np.median(prog_frac[~mask]) * 100,
            np.median(ratio[mask]), np.median(ratio[~mask]),
            np.median(defl_plasma), np.median(defl_non),
            np.median(zb_ln), np.median(zb_pr),
        ],
    })
    savecsv(summary, "test2_composition_summary.csv", index=False)
    print(f"  median plasma deflation (total/area): {np.median(defl_plasma):.3f}  "
          f"(non-plasma {np.median(defl_non):.3f})")


# ── Test 3 — sensitivity to sc.pp.scale before PCA ────────────────────────────
def _run_variant(filtered, scale, seed):
    ad = filtered.copy()
    if scale:
        sc.pp.scale(ad, max_value=10)
    sc.pp.pca(ad, n_comps=CONFIG["n_comps"], use_highly_variable=False, random_state=seed)
    ho = hm.run_harmony(ad.obsm["X_pca"].astype("float64"), ad.obs,
                        CONFIG["batch_key"], random_state=seed)
    z = ho.Z_corr
    ad.obsm["X_pca_harmony"] = z.T if z.shape[0] != ad.n_obs else z
    sc.pp.neighbors(ad, use_rep="X_pca_harmony", n_neighbors=CONFIG["n_neighbors"],
                    n_pcs=CONFIG["n_comps"], random_state=seed)
    labels = {}
    for res in CONFIG["resolutions"]:
        sc.tl.leiden(ad, resolution=res, flavor="igraph", n_iterations=2,
                     directed=False, random_state=seed, key_added="_l")
        labels[res] = ad.obs["_l"].to_numpy().copy()
    return labels, ad.varm["PCs"], np.asarray(ad.var_names)


def test3_scaling(filtered):
    print("Test 3 — scaling sensitivity (3 pipeline runs; this is the slow test)")
    base, pcs_ns, vnames = _run_variant(filtered, scale=False, seed=0)
    floor, _, _ = _run_variant(filtered, scale=False, seed=1)   # noise floor
    scaled, pcs_s, _ = _run_variant(filtered, scale=True, seed=0)

    rows = []
    for res in CONFIG["resolutions"]:
        rows.append({
            "resolution": res,
            "ari_scale_vs_noscale": adjusted_rand_score(base[res], scaled[res]),
            "nmi_scale_vs_noscale": normalized_mutual_info_score(base[res], scaled[res]),
            "ari_noise_floor": adjusted_rand_score(base[res], floor[res]),
            "nmi_noise_floor": normalized_mutual_info_score(base[res], floor[res]),
            "n_clusters_noscale": len(np.unique(base[res])),
            "n_clusters_scale": len(np.unique(scaled[res])),
        })
    df = pd.DataFrame(rows)
    savecsv(df, "test3_scaling_similarity.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, metric in zip(axes, ["ari", "nmi"]):
        ax.plot(df["resolution"], df[f"{metric}_scale_vs_noscale"], "-o",
                label="scale vs no-scale")
        ax.plot(df["resolution"], df[f"{metric}_noise_floor"], "-o",
                label="noise floor (seed)")
        ax.set_xlabel("Leiden resolution"); ax.set_ylabel(metric.upper())
        ax.set_ylim(0, 1.02); ax.legend()
    axes[0].set_title(f"{CONFIG['analysis']}: does scaling change clustering?")
    savefig(fig, "scaling", "similarity_vs_resolution")

    # Split/merge contingency at the reference resolution.
    ref = 0.6 if 0.6 in base else CONFIG["resolutions"][len(CONFIG["resolutions"]) // 2]
    ct = pd.crosstab(pd.Series(base[ref], name="no_scale"),
                     pd.Series(scaled[ref], name="scale"))
    savecsv(ct, f"test3_contingency_res{ref}.csv")
    fig, ax = plt.subplots(figsize=(max(5, ct.shape[1] * 0.4),
                                    max(4, ct.shape[0] * 0.4)))
    im = ax.imshow(ct.values / ct.values.sum(1, keepdims=True), aspect="auto",
                   cmap="Blues")
    ax.set_xlabel("scale clusters"); ax.set_ylabel("no-scale clusters")
    ax.set_title(f"cluster overlap at res {ref} (row-normalized)")
    fig.colorbar(im, ax=ax)
    savefig(fig, "scaling", f"contingency_res{ref}")

    # Top PCA loadings, no-scale vs scale — do a few loud genes dominate the PCs?
    fig, axes = plt.subplots(2, 3, figsize=(13, 6))
    for pc in range(3):
        for row, (pcs, tag) in enumerate([(pcs_ns, "no-scale"), (pcs_s, "scale")]):
            load = pcs[:, pc]
            top = np.argsort(-np.abs(load))[:10][::-1]
            axes[row, pc].barh(range(10), load[top])
            axes[row, pc].set_yticks(range(10))
            axes[row, pc].set_yticklabels(vnames[top], fontsize=7)
            axes[row, pc].set_title(f"PC{pc + 1} — {tag}")
    fig.tight_layout()
    savefig(fig, "scaling", "pca_loadings")
    print("  wrote scaling similarity curve, contingency, and PCA loadings")


def main():
    print(f"=== {CONFIG['analysis']} QC diagnostics ===")
    clustered = sc.read_h5ad(CLUSTERED)
    print(f"Loaded clustered: {clustered.n_obs:,} cells x {clustered.n_vars} genes")
    test1_background(clustered)
    test2_composition(clustered)
    del clustered

    filtered = sc.read_h5ad(FILTERED)
    test3_scaling(filtered)
    print(f"\nOutputs → {CONFIG['fig_dir']}  and  {TAB_DIR}")


if __name__ == "__main__":
    main()
