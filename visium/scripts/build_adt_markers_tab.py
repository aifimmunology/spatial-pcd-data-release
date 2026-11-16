#!/usr/bin/env python3
"""
Add an "ADT markers" tab to data/visium data tracking.xlsx.

Two columns:
  markers      - every antibody marker id, taken verbatim from each detected
                 panel (whitespace trimmed)
  panel group  - which detected panel that marker row comes from

Layout: one row per marker PER detected panel group (markers repeat across
groups). The three groups are the distinct feature_reference.csv panels found
by compare_adt_panels.py:
  Group 1 (conv+CST, 65)  - conventional 35 + CST bridge 30  (12 runs)
  Group 2 (conv, 35)      - conventional panel only          (11 runs)
  Group 3 (CST, 30)       - CST bridge panel only            (4 runs)

Marker order within each group preserves the feature_reference.csv row order.
Existing sheets (main / supplemental / visiuim-geo) are left untouched; an
existing "ADT markers" sheet is replaced.
"""

import os

import pandas as pd
from openpyxl import load_workbook

BASE = "/home/workspace/PCD_data/visium/data"
XLSX = os.path.join(BASE, "visium data tracking.xlsx")
PANEL_DIR = os.path.join(BASE, "adt-panels")
SHEET = "ADT markers"

# One representative feature_reference.csv per detected panel group, with the
# label to record in the "panel group" column.
GROUPS = [
    ("Group 1 (conv+CST, 65)", "SB008_EXP-01095_TIS05681-001_sr_outs__feature_reference.csv"),
    ("Group 2 (conv, 35)", "EXP00875_Visium_BM__feature_reference.csv"),
    ("Group 3 (CST, 30)", "TIS04788_EXP00963_AACR_SB5__feature_reference.csv"),
]


def main():
    rows = []
    for label, fname in GROUPS:
        df = pd.read_csv(os.path.join(PANEL_DIR, fname))
        for marker in df["id"]:
            rows.append({"markers": str(marker).strip(), "panel group": label})
        print(f"{label}: {len(df)} markers")

    out = pd.DataFrame(rows, columns=["markers", "panel group"])
    print(f"\nTotal rows: {len(out)}")

    book = load_workbook(XLSX)
    if SHEET in book.sheetnames:
        del book[SHEET]
        print(f'Replaced existing "{SHEET}" sheet.')
    book.save(XLSX)

    with pd.ExcelWriter(XLSX, engine="openpyxl", mode="a") as writer:
        out.to_excel(writer, sheet_name=SHEET, index=False)

    final = pd.ExcelFile(XLSX).sheet_names
    print(f"Sheets now: {final}")


if __name__ == "__main__":
    main()
