#!/usr/bin/env python3
"""Assemble ONE staged Space Ranger run into the resolved GEO deliverable layout
(reference/final_visium_structure_resolved.txt) -- Task 2.

Input: a run staged by stage_run_outs.py under
    data/visium-data/staged/<RUN_ID>/outs/{...,spatial/...}
Output: the per-sample deliverable tree under
    data/visium-data/geo_build/<SAMPLE_ID>/

The combined Space Ranger matrix holds both Gene Expression and Antibody
Capture features; the layout's separate GEX / ADT h5s are DERIVED by splitting
on matrix/features/feature_type while preserving native 10x CellRanger h5 format
(matrix-group layout) so scanpy/Seurat read them directly. Everything else is a
copy (some renamed to the sample-prefixed deliverable name).

LOCAL ONLY -- writes to disk under data/. No watchfolder upload (that step is
outward-facing and gated on comps-assigned sample IDs).

Resolved layout decisions baked in here:
  - ADT_CST_bridge/  -- OMITTED: no CST/bridge source exists in any run, so the
                        bucket is resolved out of the deliverable (Task 3).
  - <sample>_H&E.tiff -- the CytAssist image (outs/spatial/cytassist_image.tiff)
                        is used as the H&E for now; EVOS microscopy sourcing
                        remains Task 6.
Both decisions are logged so they are explicit, not silent.

Naming: the long, searchable id (visium_<tis>_<exp>_geo-ready) is used ONLY for
the top folder (and the zip built from it). Files inside keep their bare Space
Ranger names; the folder disambiguates them. The comps-assigned id (Task 5) will
replace the <tis>_<exp> stem later.

Usage:
    python assemble_geo_run.py [RUN_ID] [FOLDER_NAME]
    # defaults: RUN_ID=TIS06377-001_EXP-01092,
    #           FOLDER_NAME=visium_<run_id lowercased>_geo-ready
"""

import os
import shutil
import sys

import h5py
import numpy as np
import scipy.sparse as sp

BASE = "/home/workspace/PCD_data/visium/data/visium-data"
STAGE_ROOT = os.path.join(BASE, "staged")
BUILD_ROOT = os.path.join(BASE, "geo_build")

# Per-feature datasets always carried into a split matrix. (_all_tag_keys is
# copied verbatim; target_sets handled specially -- GEX only, remapped.)
PER_FEATURE = ["id", "name", "feature_type", "genome", "read", "pattern",
               "sequence", "isotype_control", "secondary_name"]


def split_matrix(in_path, out_path, keep_feature_type):
    """Write a feature_type-subset of a 10x matrix h5, preserving native format.

    Returns (n_features_kept, nnz_kept, nnz_total)."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with h5py.File(in_path, "r") as fin:
        m = fin["matrix"]
        data, indices, indptr = m["data"][:], m["indices"][:], m["indptr"][:]
        shape = m["shape"][:]
        csc = sp.csc_matrix((data, indices, indptr), shape=tuple(shape))
        ftype = np.array([x.decode() for x in m["features/feature_type"][:]])
        mask = ftype == keep_feature_type
        if not mask.any():
            raise ValueError(f"no '{keep_feature_type}' features in {in_path}")
        old2new = -np.ones(len(ftype), dtype=np.int64)
        old2new[mask] = np.arange(mask.sum())
        sub = csc[mask, :].tocsc()
        sub.sort_indices()

        with h5py.File(out_path, "w") as fo:
            for k, v in fin.attrs.items():
                fo.attrs[k] = v
            g = fo.create_group("matrix")
            g.create_dataset("barcodes", data=m["barcodes"][:])
            g.create_dataset("data", data=sub.data.astype(m["data"].dtype))
            g.create_dataset("indices", data=sub.indices.astype(m["indices"].dtype))
            g.create_dataset("indptr", data=sub.indptr.astype(m["indptr"].dtype))
            g.create_dataset("shape",
                             data=np.array([int(mask.sum()), csc.shape[1]],
                                           dtype=m["shape"].dtype))
            ff = g.create_group("features")
            ff.create_dataset("_all_tag_keys", data=m["features/_all_tag_keys"][:])
            for col in PER_FEATURE:
                ff.create_dataset(col, data=m["features/" + col][mask])
            # Probe target_sets index into the GEX feature block only. Keep them
            # for the GEX split (remapped to new positions); drop for ADT.
            if keep_feature_type == "Gene Expression" and "target_sets" in m["features"]:
                ts = ff.create_group("target_sets")
                for tname, tds in m["features/target_sets"].items():
                    raw = tds[:]
                    # target_set indices may be stored as ints (raw matrix) or
                    # byte-strings of ints (filtered matrix, CR 3.0.0 quirk).
                    if raw.dtype.kind in ("i", "u"):
                        old = raw.astype(np.int64)
                    else:
                        old = np.array([int(x) for x in raw], dtype=np.int64)
                    new = old2new[old]
                    if (new < 0).any():
                        raise ValueError(f"target_set '{tname}' references "
                                         f"non-GEX features -- cannot remap")
                    if np.array_equal(new, old):
                        # identity remap (GEX is the contiguous index prefix):
                        # preserve the original dataset verbatim (dtype + bytes).
                        ts.create_dataset(tname, data=raw)
                    else:
                        ts.create_dataset(tname, data=new.astype(np.uint32))
        return int(mask.sum()), int(sub.nnz), int(csc.nnz)


def copy_as(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copy2(src, dst)
    return os.path.getsize(dst)


def deliverable_id(run_id):
    """run id (e.g. 'TIS06377-001_EXP-01092') -> searchable deliverable folder /
    zip name (e.g. 'visium_tis06377-001_exp-01092_geo-ready')."""
    return f"visium_{run_id.lower()}_geo-ready"


def main():
    run_id = sys.argv[1] if len(sys.argv) > 1 else "TIS06377-001_EXP-01092"
    folder = sys.argv[2] if len(sys.argv) > 2 else deliverable_id(run_id)

    outs = os.path.join(STAGE_ROOT, run_id, "outs")
    spatial = os.path.join(outs, "spatial")
    if not os.path.isdir(outs):
        print(f"ERROR: staged run not found: {outs}\n  run stage_run_outs.py first.")
        return

    # Naming: the long, searchable id lives ONLY on the top folder (and the zip
    # built from it). Files inside keep their bare Space Ranger names; the folder
    # tree disambiguates them. The two DERIVED ADT matrices keep an `ADT_` marker
    # (a content descriptor from the template) so they don't collide with the GEX
    # matrices of the same basename.
    dest = os.path.join(BUILD_ROOT, folder)
    ge = os.path.join(dest, "gene_expression")
    print(f"Assembling run '{run_id}' -> '{folder}'")
    print(f"  source: {outs}")
    print(f"  dest:   {dest}\n")

    log = []  # (relative deliverable path, bytes, note)

    def have(name, sub="outs"):
        p = os.path.join(spatial if sub == "spatial" else outs, name)
        return p if os.path.exists(p) else None

    # ---- gene_expression: GEX-only matrix splits ----
    n, nnz, tot = split_matrix(os.path.join(outs, "raw_feature_bc_matrix.h5"),
                               os.path.join(ge, "raw", "raw_feature_bc_matrix.h5"),
                               "Gene Expression")
    log.append(("gene_expression/raw/raw_feature_bc_matrix.h5",
                None, f"GEX split: {n} feats, {nnz}/{tot} nnz"))
    n, nnz, tot = split_matrix(os.path.join(outs, "filtered_feature_bc_matrix.h5"),
                               os.path.join(ge, "processed", "filtered_feature_bc_matrix.h5"),
                               "Gene Expression")
    log.append(("gene_expression/processed/filtered_feature_bc_matrix.h5",
                None, f"GEX split: {n} feats, {nnz}/{tot} nnz"))

    # molecule_info: copy as-is into processed/
    copy_as(os.path.join(outs, "molecule_info.h5"),
            os.path.join(ge, "processed", "molecule_info.h5"))
    log.append(("gene_expression/processed/molecule_info.h5", None, "copy"))

    # spatial/: bare copies (aligned_tissue_image.jpg optional -- newer SR only)
    SPATIAL_FILES = [
        "tissue_hires_image.png", "tissue_lowres_image.png", "tissue_positions.csv",
        "scalefactors_json.json", "aligned_fiducials.jpg", "detected_tissue_image.jpg",
        "aligned_tissue_image.jpg", "spatial_enrichment.csv",
    ]
    for name in SPATIAL_FILES:
        src = have(name, "spatial")
        if src:
            copy_as(src, os.path.join(ge, "spatial", name))
            log.append((f"gene_expression/spatial/{name}", None, "copy"))
        else:
            log.append((f"gene_expression/spatial/{name}", None,
                        "ABSENT in this run (optional)"))

    # probe/: probe assay artifacts
    for name in ("probe_set.csv", "raw_probe_bc_matrix.h5"):
        src = have(name)
        if src:
            copy_as(src, os.path.join(ge, "probe", name))
            log.append((f"gene_expression/probe/{name}", None, "copy"))

    # ---- ADT_conventional: ADT-only matrix splits + isotype norm ----
    adt = os.path.join(dest, "ADT_conventional")
    n, nnz, tot = split_matrix(os.path.join(outs, "raw_feature_bc_matrix.h5"),
                               os.path.join(adt, "raw", "ADT_raw_feature_bc_matrix.h5"),
                               "Antibody Capture")
    log.append(("ADT_conventional/raw/ADT_raw_feature_bc_matrix.h5",
                None, f"ADT split: {n} feats, {nnz}/{tot} nnz"))
    n, nnz, tot = split_matrix(os.path.join(outs, "filtered_feature_bc_matrix.h5"),
                               os.path.join(adt, "processed", "ADT_filtered_feature_bc_matrix.h5"),
                               "Antibody Capture")
    log.append(("ADT_conventional/processed/ADT_filtered_feature_bc_matrix.h5",
                None, f"ADT split: {n} feats, {nnz}/{tot} nnz"))
    src = have("isotype_normalization_factors.csv")
    if src:
        copy_as(src, os.path.join(adt, "isotype_normalization_factors.csv"))
        log.append(("ADT_conventional/isotype_normalization_factors.csv", None, "copy"))

    # ---- H&E: use the CytAssist image (resolved fallback; EVOS sourcing is
    #      Task 6, deferred). Placed at the folder root per the resolved layout.
    src = have("cytassist_image.tiff", "spatial")
    if src:
        copy_as(src, os.path.join(dest, "H&E.tiff"))
        log.append(("H&E.tiff", None, "copy from cytassist_image.tiff"))

    # ---- per-run QC + web summary (folder root level) ----
    src = have("metrics_summary.csv")
    if src:
        copy_as(src, os.path.join(dest, "metrics_summary.csv"))
        log.append(("metrics_summary.csv", None, "copy"))
    # web_summary: canonical name or alias -> normalized to web_summary.html
    ws = (have("web_summary.html") or have("web_summary_exp897.html")
          or have("web_summary_exp898.html"))
    if ws:
        copy_as(ws, os.path.join(dest, "web_summary.html"))
        log.append(("web_summary.html", None,
                    f"copy{' (normalized from alias)' if not ws.endswith('web_summary.html') else ''}"))

    # ---- intentionally omitted ----
    omitted = [
        ("ADT_CST_bridge/", "OMITTED -- no CST/bridge source in any run; "
                            "bucket resolved out of the deliverable (Task 3)"),
        ("H&E (EVOS)", "using cytassist_image.tiff for now; EVOS "
                       "microscopy sourcing remains Task 6"),
    ]

    # ---- report ----
    print("Built deliverable files:")
    for path, _sz, note in log:
        flag = "  !!" if "ABSENT" in note else "    "
        print(f"{flag}{path:62} [{note}]")
    print("\nIntentionally omitted:")
    for path, note in omitted:
        print(f"    {path:62} [{note}]")
    print(f"\nSample tree assembled under: {dest}")


if __name__ == "__main__":
    main()
