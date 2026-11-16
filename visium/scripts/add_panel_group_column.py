#!/usr/bin/env python3
"""
Add an "ADT panel group" column to the `main` tab of the tracking sheet.

Each sample's processed run (its `cellranger output file name` zip) is matched
to the panel group that compare_adt_panels.py assigned to its
feature_reference.csv. The column is inserted right after
`cellranger output file name`.

Groups (see the "ADT markers" tab / _panel_comparison.csv):
  Group 1 (conv+CST, 65)  - conventional 35 + CST bridge 30
  Group 2 (conv, 35)      - conventional panel only
  Group 3 (CST, 30)       - CST bridge panel only

Rows with no extracted panel (the 6 unresolved "unk" samples) are left blank.
openpyxl is used so existing cell formatting/other sheets are preserved.
Idempotent: re-running refreshes the column in place rather than duplicating it.
"""

import os

import pandas as pd
from openpyxl import load_workbook

BASE = "/home/workspace/PCD_data/visium/data"
XLSX = os.path.join(BASE, "visium data tracking.xlsx")
COMP = os.path.join(BASE, "adt-panels", "_panel_comparison.csv")
SHEET = "main"
SRC_COL = "cellranger output file name"
NEW_COL = "ADT panel group"

HASH2LABEL = {
    "06a75899e2": "Group 1 (conv+CST, 65)",
    "930bee5560": "Group 2 (conv, 35)",
    "73b67cad7d": "Group 3 (CST, 30)",
}


def build_stem2label():
    comp = pd.read_csv(COMP)
    out = {}
    for _, r in comp.iterrows():
        stem = r["file"].replace("__feature_reference.csv", "")
        out[stem] = HASH2LABEL[r["panel_group"]]
    return out


def label_for(zipname, stem2label):
    if not isinstance(zipname, str) or not zipname.strip():
        return ""
    stem = zipname[:-4] if zipname.lower().endswith(".zip") else zipname
    return stem2label.get(stem, "")


def main():
    stem2label = build_stem2label()
    book = load_workbook(XLSX)
    ws = book[SHEET]

    headers = {ws.cell(row=1, column=c).value: c for c in range(1, ws.max_column + 1)}
    src_idx = headers[SRC_COL]

    # Insert fresh (or reuse if already present from a prior run).
    if NEW_COL in headers:
        new_idx = headers[NEW_COL]
    else:
        new_idx = src_idx + 1
        ws.insert_cols(new_idx)
        ws.cell(row=1, column=new_idx, value=NEW_COL)
        # src column is unaffected (left of insertion point)

    matched = 0
    for r in range(2, ws.max_row + 1):
        zipname = ws.cell(row=r, column=src_idx).value
        lbl = label_for(zipname, stem2label)
        ws.cell(row=r, column=new_idx, value=lbl or None)
        if lbl:
            matched += 1

    book.save(XLSX)
    print(f'Added/updated "{NEW_COL}" after "{SRC_COL}" (col {new_idx}).')
    print(f"Assigned a group to {matched} rows; {ws.max_row - 1 - matched} left blank.")


if __name__ == "__main__":
    main()
