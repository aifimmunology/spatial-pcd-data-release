#!/usr/bin/env python3
"""
Stage 01 — QC inspection of the RAW data. No filtering, no transformation.

Annotates QC metrics, then previews what _common.CUTOFFS would drop: a table
per sample, the distributions with each cut drawn on them, and a spatial map of
which spots go. Edit CUTOFFS in _common.py and re-run to see the consequences.

Reads:  h5ad/raw.h5ad
Writes: h5ad/qc.h5ad, tables/01_*.csv, figures/01_*.{png,pdf}

    python visium/scripts/spatial/01_qc.py
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from _common import (
    CUTOFFS,
    QC_H5AD,
    QC_LABELS,
    QC_METRICS,
    RAW_H5AD,
    SAMPLE_KEY,
    cutoffs_frame,
    dropped_mask,
    read_h5ad,
    save_fig,
    save_table,
    write_h5ad,
)

# QC metric -> the CUTOFFS key that cuts on it
METRIC_KEY = dict(zip(QC_METRICS, ["min_genes", "min_counts", "max_pct_mt"]))

adata = read_h5ad(RAW_H5AD)
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False,
                           inplace=True)

samples = list(adata.obs[SAMPLE_KEY].cat.categories)
obs = adata.obs

# --- What the cutoffs would drop -------------------------------------------
drop = dropped_mask(obs, CUTOFFS)

rows = []
for sample in samples:
    sub = obs[obs[SAMPLE_KEY] == sample]
    t = CUTOFFS[sample]
    rows.append(pd.Series({
        "n_spots": len(sub),
        "fail_genes": int((sub["n_genes_by_counts"] <= t["min_genes"]).sum()),
        "fail_counts": int((sub["total_counts"] <= t["min_counts"]).sum()),
        "fail_mt": int((sub["pct_counts_mt"] >= t["max_pct_mt"]).sum()),
        "dropped": int(drop[sub.index].sum()),
    }, name=sample))

preview = cutoffs_frame(CUTOFFS).join(pd.DataFrame(rows))
preview["kept"] = preview["n_spots"] - preview["dropped"]
preview["dropped_pct"] = (100 * preview["dropped"] / preview["n_spots"]).round(2)

print("Cutoff preview — nothing is filtered at this stage.")
print("(the three fail_* columns overlap, so they do not sum to dropped)")
print(preview.to_string())
print(f"\n{int(drop.sum()):,} of {len(drop):,} spots would drop ({100 * drop.mean():.2f}%)")

save_table(preview, "01_cutoff_preview.csv")
save_table(obs.groupby(SAMPLE_KEY, observed=True)[QC_METRICS]
           .describe().round(2), "01_qc_summary_per_sample.csv")

# --- Distributions per sample, with that sample's cut on each --------------
fig, axes = plt.subplots(len(samples), 3, figsize=(16, 3.4 * len(samples)), squeeze=False)
for row, sample in zip(axes, samples):
    sub = obs[obs[SAMPLE_KEY] == sample]
    for ax, metric, label in zip(row, QC_METRICS, QC_LABELS):
        cut = CUTOFFS[sample][METRIC_KEY[metric]]
        ax.hist(sub[metric], bins=80, color="steelblue")
        ax.axvline(cut, color="crimson", lw=2)
        ax.set_yscale("log")
        ax.set_xlim(0, max(np.percentile(sub[metric], 99), cut * 2))
        ax.set_xlabel(label, fontsize=9)
        ax.set_title(f"{sample} — cut at {cut:g}", fontsize=10)
fig.suptitle("QC distributions per sample; red = that sample's cutoff", y=1.0)
fig.tight_layout()
save_fig("01_qc_distributions")

# --- Which spots go ---------------------------------------------------------
adata.obs["qc_status"] = pd.Categorical(np.where(drop, "dropped", "kept"),
                                        categories=["kept", "dropped"])

fig, axes = plt.subplots(2, 3, figsize=(18, 12), squeeze=False)
for ax, sample in zip(axes.flatten(), samples):
    sub = adata[adata.obs[SAMPLE_KEY] == sample]
    n = int((sub.obs["qc_status"] == "dropped").sum())
    sc.pl.spatial(sub, color="qc_status", library_id=str(sub.obs["library_id"].iloc[0]),
                  palette={"kept": "#2166AC", "dropped": "#D6191C"},
                  title=f"{sample} — {n:,} dropped ({100 * n / sub.n_obs:.0f}%)",
                  ax=ax, show=False)
fig.suptitle("Spots the cutoffs would remove (red)", y=1.0)
fig.tight_layout()
save_fig("01_qc_spatial")

del adata.obs["qc_status"]  # cutoff-dependent; not part of the checkpoint
write_h5ad(adata, QC_H5AD)
