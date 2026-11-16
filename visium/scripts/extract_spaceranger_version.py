"""Extract the Space Ranger version used for each processed Visium run and write
it as a new column in the `main` tab of the tracking sheet.

For every processed-run zip in PCD_Spatial, the script mounts the zip via the
HISE FUSE cache and reads ONLY small housekeeping files from its central
directory -- the actual zip bytes are never fully fetched. The version is read
from three candidate locations inside the zip, tried in order:

  1. `_version`    -- a short text file whose entire content is the version
                      string (e.g. "2.1.1"); present in some Space Ranger
                      outputs but not all.
  2. `_cmdline`    -- the raw shell command that launched the run.  Space Ranger
                      is typically invoked as a versioned binary, so the binary
                      path contains the version (e.g. "/opt/spaceranger-2.1.1/
                      bin/spaceranger count ...").  The script looks for the
                      pattern "spaceranger-X.Y.Z" anywhere in that text.
  3. `_invocation` -- the Martian-format file recording every pipeline
                      parameter.  Space Ranger sometimes records the version
                      as `pipelines_version = "X.Y.Z"`.

If none of the three candidates yields a version string, the row is left blank
and printed as a warning.

Output: the `main` tab of the tracking sheet gains (or refreshes) a column
named "spaceranger version" inserted immediately after "cellranger output file
name".  Rows with no zip filename (the two samples not yet processed through
Space Ranger) are left blank.

Run from the repo root:
    python visium/scripts/extract_spaceranger_version.py
"""

import glob
import os
import re
import time
import zipfile

import hisepy as hp
import pandas as pd
from openpyxl import load_workbook

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = "/home/workspace/PCD_data/visium/data"
XLSX = os.path.join(BASE, "visium data tracking.xlsx")
LIVE_CSV = os.path.join(BASE, "project-stores", "PCD_Spatial_contents.csv")

STORE = "PCD_Spatial"
SHEET = "main"
SRC_COL = "cellranger output file name"
NEW_COL = "spaceranger version"

# Regex patterns used when parsing text members of the zip.
# Matches "spaceranger-2.1.1" anywhere in a path or command line.
RE_PATH_VER = re.compile(r"spaceranger[-_](\d+\.\d+\.\d+)", re.I)
# Matches `pipelines_version = "2.1.1"` in a Martian invocation file.
RE_INV_VER = re.compile(r'pipelines_version\s*=\s*"([^"]+)"')
# Matches `"Pipeline Version","spaceranger-2.1.1"` in web_summary.html.
RE_WEB_VER = re.compile(r'"Pipeline Version","spaceranger-(\d+\.\d+\.\d+)"', re.I)


# ---------------------------------------------------------------------------
# Version extraction helpers
# ---------------------------------------------------------------------------

def _read_member(zf, name):
    """Read a zip member by exact name and return its text, or None."""
    try:
        return zf.read(name).decode("utf-8", "replace")
    except KeyError:
        return None


def _find_member_ending(names, suffix):
    """Return the first member name that ends with `suffix`, shortest path first."""
    hits = sorted((n for n in names if n.endswith(suffix)), key=len)
    return hits[0] if hits else None


def extract_version_from_zip(zf):
    """Return the Space Ranger version string from an open ZipFile, or None.

    Tries four candidates in order:
      1. `_version`        -- full file content is the version string
      2. `_cmdline`        -- binary path contains "spaceranger-X.Y.Z"
      3. `_invocation`     -- contains `pipelines_version = "X.Y.Z"`
      4. `web_summary.html`-- contains `"Pipeline Version","spaceranger-X.Y.Z"`
    """
    names = zf.namelist()

    # --- Strategy 1: _version file ---
    # The file contains just the version string, possibly with whitespace.
    target = _find_member_ending(names, "_version")
    if target:
        txt = _read_member(zf, target)
        if txt:
            ver = txt.strip()
            # Accept only if it looks like a version (digits and dots, no spaces).
            if re.match(r"^\d+\.\d+[\.\d]*$", ver):
                return ver, "_version"

    # --- Strategy 2: _cmdline file ---
    # The binary path usually contains "spaceranger-X.Y.Z" or similar.
    target = _find_member_ending(names, "_cmdline")
    if target:
        txt = _read_member(zf, target)
        if txt:
            m = RE_PATH_VER.search(txt)
            if m:
                return m.group(1), "_cmdline"

    # --- Strategy 3: _invocation file ---
    # Martian invocation may record pipelines_version = "X.Y.Z".
    target = _find_member_ending(names, "_invocation")
    if target:
        txt = _read_member(zf, target)
        if txt:
            m = RE_INV_VER.search(txt)
            if m:
                return m.group(1), "_invocation"

    # --- Strategy 4: web_summary.html ---
    # Outputs-only zips (no housekeeping files) still contain the web summary,
    # which records: "Pipeline Version","spaceranger-X.Y.Z"
    target = _find_member_ending(names, "web_summary.html")
    if target:
        txt = _read_member(zf, target)
        if txt:
            m = RE_WEB_VER.search(txt)
            if m:
                return m.group(1), "web_summary.html"

    return None, None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # Build a map from zip basename -> store object (first hit wins on dups).
    live = pd.read_csv(LIVE_CSV)
    live["base"] = live["name"].str.split("/").str[-1]
    by_base = {}
    for _, r in live.iterrows():
        by_base.setdefault(r["base"], r)

    # Load the tracking sheet to find which zips to inspect.
    main_df = pd.read_excel(XLSX, sheet_name=SHEET)
    rows = [
        (str(r[SRC_COL]))
        for _, r in main_df.iterrows()
        if isinstance(r[SRC_COL], str) and r[SRC_COL].strip()
    ]
    # Deduplicate while preserving order (the same zip won't appear twice, but
    # be defensive).
    seen = set()
    unique_zips = []
    for z in rows:
        if z not in seen:
            seen.add(z)
            unique_zips.append(z)

    matched_in_store = [z for z in unique_zips if z in by_base]
    unmatched = [z for z in unique_zips if z not in by_base]
    if unmatched:
        print(f"WARNING: {len(unmatched)} zip(s) not found in store listing:")
        for z in unmatched:
            print(f"  {z}")

    print(f"{len(matched_in_store)} zips to inspect")

    # Query each zip for its version.
    zip_to_version = {}
    for i, base in enumerate(matched_in_store, 1):
        obj = by_base[base]
        name = obj["name"]
        try:
            cands = glob.glob(os.path.expanduser(f"~/input/**/{base}"), recursive=True)
            if cands:
                zp = cands[0]
            else:
                hp.download_from_project_store(store_name=STORE, file_name=name)
                cands = glob.glob(os.path.expanduser(f"~/input/**/{base}"), recursive=True)
                if not cands:
                    raise FileNotFoundError("not found in FUSE mount after download")
                zp = cands[0]
            with zipfile.ZipFile(zp) as zf:
                ver, source = extract_version_from_zip(zf)
            zip_to_version[base] = ver or ""
            status = f"{ver} (from {source})" if ver else "NOT FOUND"
            print(f"[{i}/{len(matched_in_store)}] {base}: {status}")
        except Exception as e:
            print(f"[{i}/{len(matched_in_store)}] {base}: ERR {type(e).__name__} {str(e)[:120]}")
            zip_to_version[base] = ""
        time.sleep(2)  # avoid hammering the HISE listing endpoint

    # Write the new column to the tracking sheet, preserving all other sheets
    # and formatting.
    book = load_workbook(XLSX)
    ws = book[SHEET]

    headers = {ws.cell(row=1, column=c).value: c for c in range(1, ws.max_column + 1)}
    src_idx = headers[SRC_COL]

    # Idempotent: reuse the column if already present, otherwise insert.
    if NEW_COL in headers:
        new_idx = headers[NEW_COL]
        print(f'\nRefreshing existing "{NEW_COL}" column (col {new_idx}).')
    else:
        new_idx = src_idx + 1
        ws.insert_cols(new_idx)
        ws.cell(row=1, column=new_idx, value=NEW_COL)
        print(f'\nInserted "{NEW_COL}" after "{SRC_COL}" (col {new_idx}).')

    filled = 0
    for r in range(2, ws.max_row + 1):
        zipname = ws.cell(row=r, column=src_idx).value
        if isinstance(zipname, str) and zipname.strip():
            ver = zip_to_version.get(zipname.strip(), "")
        else:
            ver = ""
        ws.cell(row=r, column=new_idx, value=ver or None)
        if ver:
            filled += 1

    book.save(XLSX)
    print(f"Wrote version for {filled} rows; {ws.max_row - 1 - filled} left blank.")


if __name__ == "__main__":
    main()
