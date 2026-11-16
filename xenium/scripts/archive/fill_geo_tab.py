#!/usr/bin/env python
"""Fill the `geo` tab's file columns in the Xenium tracking sheet.

Each of the 91 GEO samples gets the same curated set of processed Xenium output
files, named `<sample>_<file>`. Every file in the set exists in all five analyzer
versions, so the list is identical for every sample. The two raw-data columns get
a placeholder — the raw-vs-processed design for Xenium is deferred.

`treatment` is filled from the `main` tab's "treatment timepoint" where present.

Run from the repo root:
    python xenium/scripts/fill_geo_tab.py
"""

import openpyxl

SHEET = "xenium/data/xenium data tracking.xlsx"

# Curated processed deliverables, in column order. Present in every analyzer
# version, so the same list applies to all 91 samples.
PROCESSED = [
    "cell_feature_matrix.h5",
    "cells.parquet",
    "cell_boundaries.parquet",
    "nucleus_boundaries.parquet",
    "transcripts.parquet",
    "gene_panel.json",
    "experiment.xenium",
    "metrics_summary.csv",
    "analysis_summary.html",
    "morphology.ome.tif",
]
RAW_PLACEHOLDER = "TBD"

# geo tab layout: metadata in cols 1-9, then processed, then 2 raw columns.
FIRST_PROCESSED_COL = 10
N_RAW = 2
COL_SAMPLE = 1
COL_TREATMENT = 6


def treatment_by_sample():
    ws = openpyxl.load_workbook(SHEET, data_only=True)["main"]
    hdr = {c.value: i for i, c in enumerate(ws[1]) if c.value}
    isid, itp = hdr["sample id"], hdr["treatment timepoint"]
    out = {}
    for row in ws.iter_rows(min_row=2, values_only=True):
        sid = row[isid]
        if sid and sid not in out and row[itp]:
            out[sid] = row[itp]
    return out


def main():
    treatments = treatment_by_sample()

    wb = openpyxl.load_workbook(SHEET)
    ws = wb["geo"]

    raw_start = FIRST_PROCESSED_COL + len(PROCESSED)
    last_col = raw_start + N_RAW - 1

    # Rewrite the header for the processed + raw block, clearing any leftover
    # columns from the narrower original layout.
    for offset, _ in enumerate(PROCESSED):
        ws.cell(row=1, column=FIRST_PROCESSED_COL + offset, value="*processed data file")
    for i in range(N_RAW):
        ws.cell(row=1, column=raw_start + i, value="*raw data file")
    for col in range(last_col + 1, ws.max_column + 1):
        ws.cell(row=1, column=col, value=None)

    for r in range(2, ws.max_row + 1):
        sample = ws.cell(row=r, column=COL_SAMPLE).value
        if not sample:
            continue
        for offset, fname in enumerate(PROCESSED):
            ws.cell(row=r, column=FIRST_PROCESSED_COL + offset, value=f"{sample}_{fname}")
        for i in range(N_RAW):
            ws.cell(row=r, column=raw_start + i, value=RAW_PLACEHOLDER)
        if not ws.cell(row=r, column=COL_TREATMENT).value and sample in treatments:
            ws.cell(row=r, column=COL_TREATMENT, value=treatments[sample])
        for col in range(last_col + 1, ws.max_column + 1):
            ws.cell(row=r, column=col, value=None)

    wb.save(SHEET)
    print(f"Filled {ws.max_row - 1} rows: {len(PROCESSED)} processed columns + "
          f"{N_RAW} raw placeholders. Treatment filled from main where available.")


if __name__ == "__main__":
    main()
