#!/usr/bin/env python3
"""
Assemble a folder of all ADT panels.

For every processed Visium Space Ranger run, pull ONLY
`outs/feature_reference.csv` (the antibody panel definition) out of its zip in
the HISE project store `PCD_Spatial` and drop it, sample-prefixed, into one
folder:  data/adt-panels/.  Then compare the panels so we can confirm whether
every run shares the same antibody panel or any run differs (feeds the
`ADT_CST_bridge` decision in the README).

Why this works without unpacking each zip
------------------------------------------
`download_from_project_store` exposes the zip through HISE's read-only FUSE
cache (under ~/input/...); it does not consume writable disk.  `zipfile` reads
the archive's central directory and extracts only the one member we want, so we
never unpack the multi-GB run.

Source of truth for which zips to pull
--------------------------------------
crosscheck/adt_presence.csv  (rows with tree_found == 'Y')  -> distinct zip
basenames; the full HISE key for each comes from
project-stores/PCD_Spatial_contents.csv.

Idempotent: a panel already present in the output folder is skipped, so the
script can be re-run to fill in gaps.
"""

import glob
import os
import time
import zipfile

import pandas as pd

import hisepy as hp

BASE = "/home/workspace/PCD_data/visium/data"
ADT_PRESENCE = os.path.join(BASE, "crosscheck", "adt_presence.csv")
PCD_CONTENTS = os.path.join(BASE, "project-stores", "PCD_Spatial_contents.csv")
OUT_DIR = os.path.join(BASE, "adt-panels")
STORE_NAME = "PCD_Spatial"

TARGET_BASENAME = "feature_reference.csv"
CANONICAL_SUFFIX = "outs/feature_reference.csv"
SLEEP_BETWEEN = 2.0  # be gentle: download re-lists the big store internally


def build_zip_keymap():
    """Map each processed zip basename -> its full HISE key in PCD_Spatial.

    A row in adt_presence may reference a zip by basename; the full key (with
    the AIFI-<timestamp>/ prefix) lives in PCD_Spatial_contents.csv.  When a
    basename was uploaded more than once (same name, different timestamp), we
    prefer the key whose timestamp matches the tree_file recorded in
    adt_presence, else fall back to the first.
    """
    adt = pd.read_csv(ADT_PRESENCE)
    pcd = pd.read_csv(PCD_CONTENTS)
    pcd["base"] = pcd["name"].apply(lambda n: os.path.basename(str(n)))

    base_to_keys = {}
    for _, row in pcd.iterrows():
        base_to_keys.setdefault(row["base"], []).append(row["name"])

    proc = adt[adt["tree_found"] == "Y"]
    result = {}  # zip_basename -> (full_key, tree_file)
    for _, r in proc.iterrows():
        z = r.get("cellranger_output_file_name")
        if not isinstance(z, str) or not z.strip():
            continue
        if z in result:
            continue
        keys = base_to_keys.get(z, [])
        if not keys:
            print(f"  !! no PCD_Spatial key for {z}")
            continue
        chosen = keys[0]
        if len(keys) > 1:
            tf = str(r.get("tree_file", ""))
            stamp = tf.split("__", 1)[0] if "__" in tf else ""
            for k in keys:
                if stamp and stamp in k:
                    chosen = k
                    break
            print(f"  ~ {z}: {len(keys)} uploads, using {chosen}")
        result[z] = chosen
    return result


def locate_downloaded(file_key):
    """Find where hisepy placed the file under the read-only FUSE mount."""
    pattern = os.path.expanduser(os.path.join("~", "input", "**", file_key))
    matches = glob.glob(pattern, recursive=True)
    if not matches:
        return None
    return max(matches, key=os.path.getmtime)


def extract_one(file_key, out_path):
    """Download the zip and extract the canonical feature_reference.csv."""
    ok = hp.download_from_project_store(store_name=STORE_NAME, file_name=file_key)
    if not ok:
        print("    download reported failure")
        return False
    zip_path = locate_downloaded(file_key)
    if not zip_path or not os.path.exists(zip_path):
        print("    could not locate downloaded zip under ~/input/")
        return False

    with zipfile.ZipFile(zip_path) as zf:
        members = [m for m in zf.namelist() if m.endswith(TARGET_BASENAME)]
        if not members:
            print("    no feature_reference.csv in zip (GEX-only?)")
            return False
        # Prefer the canonical outs/ copy over any pipeline intermediates.
        canonical = [m for m in members if m.endswith(CANONICAL_SUFFIX)]
        member = canonical[0] if canonical else members[0]
        if not canonical:
            print(f"    note: using non-canonical member {member}")
        with zf.open(member) as src, open(out_path, "wb") as dst:
            dst.write(src.read())
    return True


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    keymap = build_zip_keymap()
    print(f"\n{len(keymap)} processed zips to pull feature_reference.csv from.\n")

    extracted, skipped, failed = [], [], []
    for i, (zip_base, key) in enumerate(sorted(keymap.items()), 1):
        stem = os.path.splitext(zip_base)[0]
        out_path = os.path.join(OUT_DIR, f"{stem}__feature_reference.csv")
        print(f"[{i}/{len(keymap)}] {zip_base}")
        if os.path.exists(out_path):
            print("    already have it; skipping")
            skipped.append(zip_base)
            continue
        try:
            if extract_one(key, out_path):
                print(f"    -> {os.path.basename(out_path)}")
                extracted.append(zip_base)
            else:
                failed.append(zip_base)
        except Exception as e:  # keep going; one bad zip shouldn't stop the run
            print(f"    ERROR: {e}")
            failed.append(zip_base)
        time.sleep(SLEEP_BETWEEN)

    print("\n=== summary ===")
    print(f"extracted: {len(extracted)}  skipped(existing): {len(skipped)}  failed: {len(failed)}")
    if failed:
        print("failed:")
        for z in failed:
            print(f"  {z}")


if __name__ == "__main__":
    main()
