#!/usr/bin/env python
"""Generate zip-tree .md files for the two newly-located SB009 / EXP-01130
Visium samples, reusing the machinery in visium_zip_tree.py.

These are the samples added to 'visium data tracking.xlsx' (main sheet):

    id            cellranger output file name           project store
    TIS06381-001  SB009-TIS06381-002-EXP-01130.zip      visium_data
    TIS06381-002  SB009-TIS06381-001-EXP-01130.zip      visium_data

They live in the 'visium_data' project store (NOT PCD_Spatial), so we point
visium_zip_tree at that store. The exact in-store path (AIFI-prefixed) is
resolved at run time by listing the store and matching on basename, so we do
not hard-code a path we have not verified.

Output goes to the same visium-zip-trees/ dir as the rest, named after each
zip's basename.

Requires a working HISE session (list/download must not 401).
"""

import os

import hisepy as hp

import visium_zip_tree as vzt

STORE_NAME = "visium_data"

# Basenames we are looking for, taken verbatim from the tracking sheet.
WANTED_BASENAMES = [
    "SB009-TIS06381-001-EXP-01130.zip",
    "SB009-TIS06381-002-EXP-01130.zip",
]


def resolve_paths():
    """Return the full in-store paths for WANTED_BASENAMES.

    Lists the store once and matches each wanted basename against the 'name'
    column. Raises if the store can't be listed or a wanted file is missing.
    """
    df = hp.list_files_in_project_store(STORE_NAME)
    if df is None or "name" not in df:
        raise RuntimeError(f"could not list project store {STORE_NAME!r}")
    names = list(df["name"])
    resolved = []
    for base in WANTED_BASENAMES:
        hits = [n for n in names if os.path.basename(n) == base]
        if not hits:
            print(f"  WARNING: {base} not found in store {STORE_NAME!r}")
            continue
        if len(hits) > 1:
            print(f"  NOTE: {base} matched {len(hits)} entries; using all")
        resolved.extend(hits)
    return resolved


def main():
    os.makedirs(vzt.TREE_DIR, exist_ok=True)
    # Point the shared machinery at the visium_data store (used for both the
    # download call and the tree-file header).
    vzt.STORE_NAME = STORE_NAME

    files = resolve_paths()
    if not files:
        print("No matching files resolved; nothing to do.")
        return

    ok, bad = [], []
    for i, file_name in enumerate(files, start=1):
        print(f"[{i}/{len(files)}] {file_name}")
        try:
            (ok if vzt.process(file_name) else bad).append(file_name)
        except Exception as err:
            print(f"  ERROR: {err}")
            bad.append(file_name)
    print(f"\nDone. {len(ok)} tree(s) written to {vzt.TREE_DIR}, {len(bad)} failed.")
    for name in bad:
        print(f"  FAILED: {name}")


if __name__ == "__main__":
    main()
