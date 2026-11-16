#!/usr/bin/env python
"""Download the EXP-00872 (TIS05396-001-001, slide V53U12-040 area A1) FASTQs
from the PCD_Spatial HISE project store, to enable a from-scratch
`spaceranger count` rerun of the only FASTQ-only / no-processed-output sample.

The FASTQs were surfaced into the readable PCD_Spatial store by the watchfolder
ingest under SB002/FC00516_done/outs/fastq_path/HNKWTDMXY/.  Both libraries are
present (matching the slide-mate EXP00875 run's library layout):
  - EXP-00872_G  (Gene Expression,  S16, L001+L002)
  - EXP-00872_P  (Antibody Capture,  S2, L001+L002)

hp.download_from_project_store stages each file in HISE's read-only FUSE cache
(~/input/...), so we locate it there and COPY (not move) it into TARGET_DIR with
its standard Illumina basename, which is what `spaceranger --libraries` expects.

Usage:  python scripts/download_exp00872_fastqs.py
"""

import glob
import os
import shutil

import pandas as pd

import hisepy as hp

STORE_NAME = "PCD_Spatial"
PCD_CONTENTS = "/home/workspace/PCD_data/visium/data/project-stores/PCD_Spatial_contents.csv"
TARGET_DIR = "/home/workspace/PCD_data/visium/data/rerun-exp00872/fastqs"
INPUT_ROOT = os.path.expanduser("~/input")


def resolve_keys():
    """Return the full PCD_Spatial keys for the EXP-00872 _G_/_P_ fastq.gz files."""
    pcd = pd.read_csv(PCD_CONTENTS)
    names = pcd["name"].astype(str)
    mask = (
        names.str.contains("EXP-00872_[GP]_", regex=True)
        & names.str.contains(r"\.fastq\.gz$", regex=True)
        & ~names.str.contains("stats|html|Reports", case=False, regex=True)
    )
    return sorted(names[mask].tolist())


def stage_locally(file_name):
    """Find the just-downloaded file under the FUSE mount and copy it to TARGET_DIR."""
    matches = glob.glob(os.path.join(INPUT_ROOT, "**", file_name), recursive=True)
    if not matches:
        return None
    src = max(matches, key=os.path.getmtime)
    dst = os.path.join(TARGET_DIR, os.path.basename(file_name))
    shutil.copy2(src, dst)  # copy: the FUSE source is read-only
    return dst


def main():
    os.makedirs(TARGET_DIR, exist_ok=True)
    keys = resolve_keys()
    print(f"{len(keys)} FASTQ files to pull into {TARGET_DIR}\n")
    ok_n, fail = 0, []
    for i, key in enumerate(keys, 1):
        base = os.path.basename(key)
        dst = os.path.join(TARGET_DIR, base)
        print(f"[{i}/{len(keys)}] {base}")
        if os.path.exists(dst) and os.path.getsize(dst) > 0:
            print("    already staged; skipping")
            ok_n += 1
            continue
        try:
            if not hp.download_from_project_store(store_name=STORE_NAME, file_name=key):
                print("    download returned False")
                fail.append(base)
                continue
            out = stage_locally(key)
            if out:
                print(f"    -> {out}  ({os.path.getsize(out)/1e6:.1f} MB)")
                ok_n += 1
            else:
                print("    ERROR: staged file not found under ~/input/")
                fail.append(base)
        except Exception as e:
            print(f"    ERROR: {e}")
            fail.append(base)

    print(f"\nDone. {ok_n} staged, {len(fail)} failed.")
    for f in fail:
        print(f"  failed: {f}")


if __name__ == "__main__":
    main()
