#!/usr/bin/env python3
"""
Stage 00 — HISE zips -> per-sample Space Ranger outs/ -> combined AnnData.

A sample whose outs/ is already extracted is read straight off disk, so re-runs
are local and offline. Only a missing sample touches HISE.

Writes: h5ad/raw.h5ad  (raw counts, no QC, no transform)

    python visium/scripts/spatial/00_load.py
"""

import glob
import os
import zipfile

import scanpy as sc

from _common import INPUT_ROOT, OUTS_BASE, RAW_H5AD, SAMPLE_KEY, SAMPLES, write_h5ad

adatas = {}
for sample_id, uuid in SAMPLES.items():
    sample_outs = os.path.join(OUTS_BASE, sample_id)

    if not os.path.isdir(sample_outs):
        import hisepy as hp  # imported lazily — only needed on a cache miss

        hp.cache_files(file_ids=[uuid])
        hits = glob.glob(os.path.join(INPUT_ROOT, "**", f"*{uuid}*"), recursive=True)
        with zipfile.ZipFile(max(hits, key=os.path.getmtime)) as zf:
            for member in zf.namelist():
                if "/outs/" in member and not member.endswith("/"):
                    target = os.path.join(
                        sample_outs, member[member.index("/outs/") + len("/outs/") :]
                    )
                    os.makedirs(os.path.dirname(target), exist_ok=True)
                    with zf.open(member) as src, open(target, "wb") as dst:
                        dst.write(src.read())

    adata = sc.read_visium(sample_outs)
    adata.var_names_make_unique()
    adata.obs_names_make_unique()
    # Space Ranger's library id (the uns['spatial'] key) is not our short sample
    # id; carry it in obs so spatial plots on a subset can name the library.
    adata.obs["library_id"] = list(adata.uns["spatial"])[0]
    adatas[sample_id] = adata

combined = sc.concat(adatas, label=SAMPLE_KEY, join="outer")
combined.obs_names_make_unique()

# sc.concat drops uns; restore it so spatial plotting works on subsets
combined.uns["spatial"] = {lib: meta for adata in adatas.values()
                           for lib, meta in adata.uns["spatial"].items()}

print(combined.obs[SAMPLE_KEY].value_counts().to_string())
write_h5ad(combined, RAW_H5AD)
