#!/usr/bin/env python
"""Download a fixed set of files from the PCD_Spatial HISE project store
into the local visium-data/ directory.

hp.download_from_project_store writes each file to a hise-managed path under
'~/input/<id>/PCD_Spatial/<uuid>/<file_name>', so after each download we locate
the file and move it into TARGET_DIR, preserving the 'AIFI-.../' name prefix as
a subfolder (some files share a basename, so the prefix keeps them distinct).

Usage:
    python visium_data_download.py
"""

import glob
import os
import shutil

import hisepy as hp

# Project store holding the data (the store is named "PCD_Spatial").
STORE_NAME = "PCD_Spatial"

# Local destination directory for the downloaded files.
TARGET_DIR = "/home/workspace/PCD_data/visium/data/visium-data"

# Where hise stages downloaded files before we move them.
INPUT_ROOT = os.path.expanduser("~/input")

# Exact file names to download, as listed by list_files_in_project_store().
FILES = [
    "AIFI-2024-03-14T17:52:07.667273266Z/TIS05398-001-007-EXP-00971-SB005.zip",
    "AIFI-2024-03-14T17:51:49.757643408Z/EXP-00971 AACR 3+4-20240229T010204Z-001.zip",
    "AIFI-2024-03-14T17:52:01.813837457Z/TIS05394-001-010-EXP-00971-SB005.zip",
    "AIFI-2024-06-10T23:08:48.250692659Z/TIS05393_EXP00971_AACR_SB5.zip",
    "AIFI-2024-06-10T23:09:28.900356701Z/TIS05395_EXP00971_AACR_SB5.zip",
    "AIFI-2024-06-10T23:12:49.667054782Z/TIS05398_EXP00971_AACR_SB5.zip",
    "AIFI-2024-06-10T23:15:03.199377095Z/TIS05395_EXP00971_AACR_SB5.zip",
    "AIFI-2024-06-10T23:16:03.199656687Z/TIS05394_EXP00971_AACR_SB5.zip",
    "AIFI-2024-06-10T23:17:13.523763688Z/TIS05393_EXP00971_AACR_SB5.zip",
    "AIFI-2024-03-14T17:51:55.870095031Z/SB005_EXP00963_rerun_outs.zip",
    "AIFI-2024-06-10T23:13:27.824091334Z/TIS05397_EXP00963_AACR_SB5.zip",
    "AIFI-2024-06-10T23:14:13.664453937Z/TIS05396_EXP00963_AACR_SB5.zip",
    "AIFI-2024-06-10T23:17:45.568843181Z/TIS05392_EXP00963_AACR_SB5.zip",
    "AIFI-2024-06-10T23:18:52.069714292Z/TIS04788_EXP00963_AACR_SB5.zip",
    "AIFI-2024-06-10T23:20:01.122865514Z/EXP00937_BM.zip",
    "AIFI-2024-06-10T23:22:06.740736449Z/EXP00898_BM.zip",
    "AIFI-2024-06-10T23:22:51.980946983Z/EXP00897_BM.zip",
    "AIFI-2024-06-10T23:21:35.580329354Z/EXP00934_BM.zip",
    "AIFI-2024-06-10T23:20:26.232014326Z/EXP00936_BM.zip",
    "AIFI-2024-06-10T23:21:01.685829928Z/EXP00935_BM.zip",
    "AIFI-2024-06-10T23:23:30.477514464Z/EXP00875_Visium_BM.zip",
    "AIFI-2024-07-03T20:20:20.619605327Z/SB008_EXP-01106_TMA00144.zip",
    "AIFI-2024-07-03T20:30:34.747740525Z/SB008_EXP-01106_TMA00147_full.zip",
]


def move_to_target(file_name):
    """Find the just-downloaded file under INPUT_ROOT and move it into TARGET_DIR.

    Returns the destination path, or None if the staged file could not be found.
    """
    matches = glob.glob(os.path.join(INPUT_ROOT, "**", file_name), recursive=True)
    if not matches:
        return None
    # If multiple stages exist, take the most recently modified.
    src = max(matches, key=os.path.getmtime)
    dst = os.path.join(TARGET_DIR, file_name)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.move(src, dst)
    return dst


def main():
    os.makedirs(TARGET_DIR, exist_ok=True)
    succeeded, failed = [], []
    for i, file_name in enumerate(FILES, start=1):
        print(f"[{i}/{len(FILES)}] Downloading '{file_name}' from '{STORE_NAME}' ...")
        try:
            ok = hp.download_from_project_store(
                store_name=STORE_NAME,
                file_name=file_name,
            )
            if not ok:
                print("  download returned False")
                failed.append(file_name)
                continue
            dst = move_to_target(file_name)
            if dst is None:
                print("  ERROR: download reported success but staged file not found")
                failed.append(file_name)
            else:
                print(f"  -> {dst}")
                succeeded.append(file_name)
        except Exception as err:  # keep going if one file fails
            print(f"  ERROR downloading '{file_name}': {err}")
            failed.append(file_name)

    print(f"\nDone. {len(succeeded)} downloaded to {TARGET_DIR}, {len(failed)} failed.")
    if failed:
        print("Failed files:")
        for name in failed:
            print(f"  {name}")


if __name__ == "__main__":
    main()
