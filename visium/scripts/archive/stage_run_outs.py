#!/usr/bin/env python3
"""Stage ONE processed Space Ranger run's deliverable-relevant `outs/` files
locally, ahead of the GEO reorganize step (Task 2).

Downloads a single run zip from the `PCD_Spatial` HISE project store, then
extracts ONLY the members the GEO deliverable draws from -- the combined
matrices, molecule_info, web_summary, feature_reference, probe artifacts,
per-run QC csvs, and the `spatial/` images/tables. The bulky non-deliverables
(possorted_genome_bam, analysis/, deconvolution/, cloupe.cloupe, the redundant
MTX matrix dirs) are deliberately left in the archive and never extracted.

Why match the CANONICAL outs only: these zips carry duplicate copies of the
same basenames in QC / cloupe / pipeline-intermediate subtrees. We locate the
one `.../outs/` directory whose direct children include the combined matrices +
feature_reference.csv, and pull members from exactly that prefix.

hisepy stages the downloaded zip on a read-only FUSE mount under ~/input/
(NOT ~/<store>/), so we locate it via glob. zipfile reads the central directory
and extracts only the members we ask for; the rest of the archive is never
unpacked. Peak writable disk ~= one run's deliverable payload (here ~650MB,
dominated by molecule_info.h5).

Usage:
    python stage_run_outs.py            # stages the default run below
"""

import glob
import os
import zipfile

import hisepy as hp

STORE_NAME = "PCD_Spatial"

# Default target: SB008 EXP-01092 TIS06377-001 -- ADT_PRESENT, no BAM (smaller
# zip), a clean modern Space Ranger layout. One of the 4 duplicate-paired runs.
RUN_KEY = ("AIFI-2024-07-03T19:50:28.897069523Z/"
           "SB008_EXP-01092_TIS06377-001_corrected_sr_outs.zip")
RUN_ID = "TIS06377-001_EXP-01092"

INPUT_ROOT = os.path.expanduser("~/input")
STAGE_ROOT = "/home/workspace/PCD_data/visium/data/visium-data/staged"

# Deliverable-relevant members, by the level they live at. These mirror the
# resolved layout (reference/final_visium_structure_resolved.txt). Some are
# optional (not every run carries aligned_tissue_image.jpg / web_summary
# aliases) -- missing optional members are reported, not fatal.
OUTS_FILES = {
    "raw_feature_bc_matrix.h5",
    "filtered_feature_bc_matrix.h5",
    "molecule_info.h5",
    "web_summary.html",
    "feature_reference.csv",
    "metrics_summary.csv",
    "isotype_normalization_factors.csv",
    "probe_set.csv",
    "raw_probe_bc_matrix.h5",
}
# Older BM runs ship a sample-suffixed web summary instead of web_summary.html.
WEB_SUMMARY_ALIASES = {"web_summary_exp897.html", "web_summary_exp898.html"}

SPATIAL_FILES = {
    "tissue_hires_image.png",
    "tissue_lowres_image.png",
    "tissue_positions.csv",
    "scalefactors_json.json",
    "aligned_fiducials.jpg",
    "detected_tissue_image.jpg",
    "aligned_tissue_image.jpg",
    "spatial_enrichment.csv",
    "cytassist_image.tiff",
}

# Files used to identify the canonical outs directory among duplicates.
CANONICAL_MARKERS = {"raw_feature_bc_matrix.h5",
                     "filtered_feature_bc_matrix.h5",
                     "feature_reference.csv"}


def find_canonical_outs(names):
    """Return the zip-internal path prefix (ending in '/outs') whose direct
    children best match the combined-matrix markers. Picks the best overlap,
    breaking ties by child count."""
    outs_children = {}  # outs_prefix -> set(direct child basenames)
    for n in names:
        parts = n.rstrip("/").split("/")
        if "outs" not in parts:
            continue
        idx = len(parts) - 1 - parts[::-1].index("outs")  # last 'outs' segment
        # direct child of this outs dir?
        if idx == len(parts) - 2:  # parts[idx] == 'outs', parts[idx+1] == child
            prefix = "/".join(parts[: idx + 1])
            outs_children.setdefault(prefix, set()).add(parts[idx + 1])
    if not outs_children:
        return None, {}
    best = max(outs_children,
               key=lambda p: (len(outs_children[p] & CANONICAL_MARKERS),
                              len(outs_children[p])))
    return best, outs_children[best]


def main():
    zip_basename = os.path.basename(RUN_KEY)
    print(f"Staging run: {RUN_ID}")
    print(f"  store key: {RUN_KEY}\n")

    # 1. Download the whole zip to the FUSE mount (no byte-range read available).
    print("Downloading zip from project store (this fetches the full archive)...")
    ok = hp.download_from_project_store(store_name=STORE_NAME, file_name=RUN_KEY)
    if not ok:
        print("  download_from_project_store returned False -- aborting.")
        return

    # 2. Locate the staged zip under ~/input/ (read-only FUSE; not ~/<store>/).
    matches = glob.glob(os.path.join(INPUT_ROOT, "**", zip_basename), recursive=True)
    if not matches:
        print(f"  ERROR: staged zip not found under {INPUT_ROOT}")
        return
    zip_path = max(matches, key=os.path.getmtime)
    print(f"  staged at: {zip_path}\n")

    # 3. Identify the canonical outs/ prefix, then extract only what we need.
    out_dir = os.path.join(STAGE_ROOT, RUN_ID, "outs")
    spatial_dir = os.path.join(out_dir, "spatial")
    os.makedirs(spatial_dir, exist_ok=True)

    wanted_outs = OUTS_FILES | WEB_SUMMARY_ALIASES
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        prefix, children = find_canonical_outs(names)
        if prefix is None:
            print("  ERROR: no canonical outs/ directory found in zip.")
            return
        print(f"  canonical outs prefix: {prefix}")
        print(f"  outs direct children: {len(children)}\n")

        extracted, missing = [], []

        def extract_member(member, dest_dir):
            dst = os.path.join(dest_dir, os.path.basename(member))
            with zf.open(member) as src, open(dst, "wb") as out:
                # stream in chunks; molecule_info.h5 is ~450MB
                while True:
                    chunk = src.read(8 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
            return os.path.getsize(dst)

        # outs-level members (exact path = prefix/<name>)
        for name in sorted(wanted_outs):
            member = f"{prefix}/{name}"
            if member in names:
                sz = extract_member(member, out_dir)
                extracted.append((name, sz))
            elif name in OUTS_FILES and name != "web_summary.html":
                missing.append(name)
        # web_summary: required, but may be present only via alias
        if not any(n.startswith("web_summary") for n, _ in extracted):
            missing.append("web_summary.html (and aliases)")

        # spatial-level members (exact path = prefix/spatial/<name>)
        for name in sorted(SPATIAL_FILES):
            member = f"{prefix}/spatial/{name}"
            if member in names:
                sz = extract_member(member, spatial_dir)
                extracted.append((f"spatial/{name}", sz))
            else:
                missing.append(f"spatial/{name}")

    # 4. Report.
    print("Extracted:")
    total = 0
    for name, sz in extracted:
        total += sz
        print(f"  {name:48} {sz/1e6:8.1f} MB")
    print(f"  {'TOTAL':48} {total/1e6:8.1f} MB")
    if missing:
        print("\nNot present (optional / alias-only / absent in this run):")
        for m in missing:
            print(f"  {m}")
    print(f"\nStaged under: {out_dir}")


if __name__ == "__main__":
    main()
