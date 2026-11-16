#!/usr/bin/env python3
"""Populate the 'xenium analyzer software version' column in the tracking sheet.

Run from the repo root with a HISE-authenticated Python environment:
    python xenium/scripts/populate_analyzer_versions.py

Strategy:
  - The version is stored in `experiment.xenium` (a JSON file) inside each
    output zip. Every sample run by the Xenium instrument in the same batch
    shares the same run timestamp (the YYYYMMDD__HHMMSS token at the end of
    the file name) and therefore the same analyzer version. So we only need to
    download one zip per unique run timestamp, not one per row.

  - The tracking sheet has 95 rows and 88 unique zip files, but only 20 unique
    run timestamps. We download one representative zip per timestamp, extract
    experiment.xenium from it, read the version, clean up, and then write that
    version to every row sharing the same timestamp.
"""

import json
import os
import re
import zipfile

import hisepy as hp
import openpyxl

TRACKING = "xenium/data/xenium data tracking.xlsx"
STORE = "PCD_Spatial"
VERSION_FIELD = "xenium analyzer software version"
FILE_NAME_FIELD = "file name"
UUID_FIELD = "hise uuid"

# `experiment.xenium` may have changed key names across analyzer versions.
# We check these in order and use the first one found.
VERSION_KEYS = ["analysis_sw_version", "software_version", "xenium_analyzer_version"]


def extract_run_timestamp(file_name: str) -> str:
    """Return the YYYYMMDD__HHMMSS run-timestamp token from a Xenium zip name.

    The file name format is:
      <ingest-prefix>/output-<instrument>__<slide>__<sample>__<YYYYMMDD>__<HHMMSS>.zip
    The last two double-underscore-delimited tokens before .zip are the run
    date and run time — together they uniquely identify one Xenium Analyzer run.
    """
    m = re.search(r"__(\d{8}__\d{6})\.zip$", file_name)
    return m.group(1) if m else "UNKNOWN"


def get_version_from_zip(zip_path: str) -> str | None:
    """Open a Xenium output zip and read the analyzer version from experiment.xenium.

    experiment.xenium is a small JSON file at the root of the output directory
    inside the zip. We extract just that file without unpacking the whole archive.
    Returns the version string, or None if the file or field is missing.
    """
    try:
        with zipfile.ZipFile(zip_path) as zf:
            # Find experiment.xenium — it may sit at the archive root or one
            # level down inside an output-... folder.
            candidates = [n for n in zf.namelist() if n.endswith("experiment.xenium")]
            if not candidates:
                print(f"  WARNING: experiment.xenium not found in {zip_path}")
                return None
            # Prefer the shortest path (most likely root-level).
            target = min(candidates, key=len)
            with zf.open(target) as fh:
                data = json.load(fh)
            for key in VERSION_KEYS:
                if key in data:
                    return str(data[key])
            print(f"  WARNING: none of {VERSION_KEYS} found in {target}")
            print(f"  Available keys: {list(data.keys())}")
            return None
    except zipfile.BadZipFile as e:
        print(f"  ERROR: bad zip file {zip_path}: {e}")
        return None


def main():
    wb = openpyxl.load_workbook(TRACKING)
    ws = wb["main"]

    hdr = {c.value: c.column for c in ws[1] if c.value is not None}
    fn_col = hdr[FILE_NAME_FIELD]
    uuid_col = hdr[UUID_FIELD]
    ver_col = hdr[VERSION_FIELD]

    # --- Step 1: collect rows grouped by run timestamp ---
    # ts_to_rows: timestamp -> list of (row_number, file_name, uuid)
    ts_to_rows: dict[str, list[tuple[int, str, str]]] = {}
    for r in range(2, ws.max_row + 1):
        fn = ws.cell(row=r, column=fn_col).value
        if not fn:
            continue
        uuid = ws.cell(row=r, column=uuid_col).value
        ts = extract_run_timestamp(fn)
        ts_to_rows.setdefault(ts, []).append((r, fn, uuid))

    print(f"Found {len(ts_to_rows)} unique run timestamps across all rows.")
    print()

    # --- Step 2: for each timestamp, download one representative zip and
    #             extract the analyzer version ---
    ts_to_version: dict[str, str] = {}

    for ts, rows in sorted(ts_to_rows.items()):
        # Pick the first row that has a UUID.
        representative = next(((fn, uuid) for _, fn, uuid in rows if uuid), None)
        if representative is None:
            print(f"[{ts}] No UUID available — skipping.")
            continue

        fn, uuid = representative
        print(f"[{ts}] Downloading via UUID {uuid}")
        print(f"      file: {fn}")

        zip_paths = hp.cache_files(file_ids=[uuid])
        if not zip_paths:
            print(f"  ERROR: cache_files returned no paths — skipping this timestamp.")
            continue
        zip_path = zip_paths[0]
        print(f"  Downloaded to: {zip_path}")

        version = get_version_from_zip(zip_path)
        if version:
            print(f"  Version: {version}")
            ts_to_version[ts] = version
        else:
            print(f"  Could not extract version.")

        # Clean up — Xenium output zips are large; remove after use.
        try:
            os.remove(zip_path)
            print(f"  Cleaned up {zip_path}")
        except OSError as e:
            print(f"  WARNING: could not remove {zip_path}: {e}")

        print()

    # --- Step 3: write versions back to every row with the matching timestamp ---
    filled = 0
    skipped_ts = []
    for ts, rows in ts_to_rows.items():
        version = ts_to_version.get(ts)
        if version is None:
            skipped_ts.append(ts)
            continue
        for row_num, _, _ in rows:
            ws.cell(row=row_num, column=ver_col, value=version)
            filled += 1

    wb.save(TRACKING)

    print(f"Wrote version to {filled} rows across {len(ts_to_version)} timestamps.")
    if skipped_ts:
        print(f"Could not resolve version for {len(skipped_ts)} timestamps:")
        for ts in skipped_ts:
            print(f"  {ts}")


if __name__ == "__main__":
    main()
