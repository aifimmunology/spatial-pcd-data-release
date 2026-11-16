#!/usr/bin/env python
"""Add two manifest columns to the tracking sheet's `main` tab.

1. `xenium analyzer software version` — the Xenium Analyzer version for each row,
   determined from the run timestamp in its output zip name (same timestamp =
   same run = same version), mapped via the version trees under
   xenium/data/xenium-zip-trees/.
2. `downstream analysis (Fig2 D-F)` — flags the six sections used in the Fig2
   D/E/F pipeline.

Run from the repo root with the HISE-authenticated environment:
    /home/workspace/environment/scanpyv3/bin/python xenium/scripts/add_main_columns.py
"""

import glob
import re
import unittest.mock as mock

import pandas as pd
pd.set_option("future.infer_string", False)

with mock.patch.dict("sys.modules", {"hisepy.training_job": mock.MagicMock(),
                                     "hisepy.ai_ray_transformer": mock.MagicMock()}):
    import hisepy as hp

import openpyxl

SHEET = "xenium/data/xenium data tracking.xlsx"
STORE = "PCD_Spatial"
TS_RE = re.compile(r"__(\d{8}__\d{6})\.zip$")

VERSION_COL = "xenium analyzer software version"
FIG2_COL = "downstream analysis (Fig2 D-F)"

FIG2_SECTIONS = {
    "TIS04788-001-010", "TIS04789-001-017", "TIS05392-001-007",
    "TIS05393-001-009", "TIS05394-001-008", "TIS05397-001-008",
}


def ts_to_version():
    out = {}
    for path in glob.glob("xenium/data/xenium-zip-trees/*.md"):
        m = re.search(r"/(\d{8}__\d{6})__(xenium-[\d.]+)\.md$", path)
        if m:
            out[m.group(1)] = m.group(2)
    return out


def main():
    id2name = dict(zip(*[hp.list_files_in_project_store(STORE)[k] for k in ("id", "name")]))
    ts2ver = ts_to_version()

    wb = openpyxl.load_workbook(SHEET)
    ws = wb["main"]
    hdr = {c.value: c.column for c in ws[1] if c.value is not None}
    iu = hdr["hise uuid"]
    isec = hdr["section id"]

    # Append after the last column that actually has a header VALUE. Do NOT use
    # ws.max_column: the sheet carries empty-but-styled cells far past the real
    # table, so max_column overshoots and would drop new columns into a gap.
    def append_col(name):
        last = max(c.column for c in ws[1] if c.value is not None)
        col = last + 1
        ws.cell(row=1, column=col, value=name)
        return col

    ver_col = hdr.get(VERSION_COL) or append_col(VERSION_COL)
    fig2_col = hdr.get(FIG2_COL) or append_col(FIG2_COL)

    n_ver = n_fig2 = 0
    for r in range(2, ws.max_row + 1):
        uuid = ws.cell(row=r, column=iu).value
        if uuid:
            name = id2name.get(uuid)
            m = TS_RE.search(name) if name else None
            if m and m.group(1) in ts2ver:
                ws.cell(row=r, column=ver_col, value=ts2ver[m.group(1)])
                n_ver += 1
        if ws.cell(row=r, column=isec).value in FIG2_SECTIONS:
            ws.cell(row=r, column=fig2_col, value="yes")
            n_fig2 += 1

    wb.save(SHEET)
    print(f"Wrote version for {n_ver} rows; flagged {n_fig2} Fig2 sections.")


if __name__ == "__main__":
    main()
