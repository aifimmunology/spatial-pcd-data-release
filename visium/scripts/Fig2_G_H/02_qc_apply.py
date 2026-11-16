#!/usr/bin/env python3
"""
Stage 02 — apply the cutoffs in _common.CUTOFFS. Spots only; genes are kept.

Same rule 01_qc.py previewed (_common.dropped_mask), so the preview and the
apply step cannot disagree. The cutoffs used are stored in uns['cutoffs'], so a
filtered object is always traceable back to the cuts that made it.

Reads:  h5ad/qc.h5ad
Writes: h5ad/filtered.h5ad, tables/02_filtering_report.csv

    python visium/scripts/spatial/02_qc_apply.py
"""

import pandas as pd

from _common import (
    CUTOFFS,
    FILTERED_H5AD,
    QC_H5AD,
    SAMPLE_KEY,
    cutoffs_frame,
    dropped_mask,
    read_h5ad,
    save_table,
    write_h5ad,
)

adata = read_h5ad(QC_H5AD)
obs = adata.obs
drop = dropped_mask(obs, CUTOFFS)

rows = []
for sample in obs[SAMPLE_KEY].cat.categories:
    sub = obs[obs[SAMPLE_KEY] == sample]
    t = CUTOFFS[sample]
    rows.append(pd.Series({
        "n_spots": len(sub),
        "fail_genes": int((sub["n_genes_by_counts"] <= t["min_genes"]).sum()),
        "fail_counts": int((sub["total_counts"] <= t["min_counts"]).sum()),
        "fail_mt": int((sub["pct_counts_mt"] >= t["max_pct_mt"]).sum()),
        "removed": int(drop[sub.index].sum()),
    }, name=sample))

report = cutoffs_frame(CUTOFFS).join(pd.DataFrame(rows))
report["kept"] = report["n_spots"] - report["removed"]

print("Per-sample cutoffs and what they removed")
print("(the three fail_* columns overlap, so they do not sum to removed):")
print(report.to_string())

adata = adata[~drop].copy()
print(f"\n{adata.n_obs:,} spots kept of {len(drop):,} ({100 * (~drop).mean():.2f}%)")

save_table(report, "02_filtering_report.csv")
adata.uns["cutoffs"] = {s: dict(t) for s, t in CUTOFFS.items()}
write_h5ad(adata, FILTERED_H5AD)
