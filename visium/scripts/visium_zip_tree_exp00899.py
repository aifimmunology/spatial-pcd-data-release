#!/usr/bin/env python
"""Generate zip-tree .md files for the two EXP-00899/900 SB003 zips in the
PCD_Spatial project store, reusing the machinery in visium_zip_tree.py.

These are the raw CytAssist image zip and the EXP-00899_G_L1 (FASTQ) zip that
showed up as NO_TREE / "No — raw CytAssist / FASTQ only" in the ADT and Space
Ranger checks. Generating their trees lets us see exactly what they contain.

Output goes to the same visium-zip-trees/ dir as the rest, named after each
zip's basename.
"""

import os

import visium_zip_tree as vzt

FILES = [
    "AIFI-2026-06-05T00:22:09.419430065Z/SB003/EXP-00899_G_L1.zip",
    "AIFI-2026-06-05T00:22:01.474862462Z/SB003/"
    "EXP-00899-900_CAVG10142_2023-12-07_12-03-56_V53U05-373_1701980849_"
    "CytAssist-20240111T192851Z-001.zip",
]


def main():
    os.makedirs(vzt.TREE_DIR, exist_ok=True)
    ok, bad = [], []
    for i, file_name in enumerate(FILES, start=1):
        print(f"[{i}/{len(FILES)}] {file_name}")
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
