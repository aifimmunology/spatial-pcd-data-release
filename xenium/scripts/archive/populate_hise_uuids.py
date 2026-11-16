#!/usr/bin/env python3
"""Populate a `hise uuid` column in the tracking sheet's `main` tab.

Run from the repo root:  python xenium/scripts/populate_hise_uuids.py

For every `main`-tab row, look up the HISE file UUID by matching the row's
`file name` (store-relative path) against a live listing of the PCD_Spatial
project store, then write it to a `hise uuid` column (appended as the last
column). Other tabs and existing formatting are preserved (openpyxl in place).
"""

import hisepy as hp
import openpyxl

TRACKING = "xenium/data/xenium data tracking.xlsx"
STORE = "PCD_Spatial"


def main():
    live = hp.list_files_in_project_store(STORE)
    name2uuid = dict(zip(live["name"], live["id"]))

    wb = openpyxl.load_workbook(TRACKING)
    ws = wb["main"]

    hdr = {c.value: c.column for c in ws[1] if c.value is not None}
    fn_col = hdr["file name"]
    uuid_col = hdr.get("hise uuid", max(hdr.values()) + 1)
    ws.cell(row=1, column=uuid_col, value="hise uuid")

    filled, missing = 0, []
    for r in range(2, ws.max_row + 1):
        fn = ws.cell(row=r, column=fn_col).value
        if not fn:
            continue
        uuid = name2uuid.get(fn)
        if uuid is None:
            missing.append((r, fn))
            continue
        ws.cell(row=r, column=uuid_col, value=uuid)
        filled += 1

    wb.save(TRACKING)
    print(f"Filled 'hise uuid' for {filled} rows; {len(missing)} unresolved.")
    for r, fn in missing:
        print(f"  [row {r}] {fn}")


if __name__ == "__main__":
    main()
