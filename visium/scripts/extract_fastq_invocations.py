"""Extract each processed Visium run's FASTQ provenance from its Space Ranger
`_invocation` (the authoritative record of what reads the run consumed).

For every processed-run zip in PCD_Spatial (matched to the tracking sheet's
`cellranger output file name`), mount it via the HISE FUSE cache and read ONLY the
top-level `_invocation` member -- zipfile pulls the central directory + that one small
member, so the full (often tens-of-GB) archive is never fetched. Parse, per library:
`library_type`, `read_path` (the flowcell / fastq dir), `sample_names` (the FASTQ prefix),
plus `sample_id` and slide/area.

Output: visium/data/crosscheck/fastq_invocations.csv  (one row per run x library).
Run from repo root: python visium/scripts/extract_fastq_invocations.py
"""
import glob
import os
import re
import time
import zipfile

import hisepy as hp
import pandas as pd

STORE = "PCD_Spatial"
LIVE = "/tmp/claude-0/-home-workspace-PCD-data-visium/c7d389ac-c071-4667-b42c-81868dbebf75/scratchpad/pcd_live.csv"
OUT = "visium/data/crosscheck/fastq_invocations.csv"

live = pd.read_csv(LIVE)
live["base"] = live["name"].str.split("/").str[-1]
main = pd.read_excel("visium/data/visium data tracking.xlsx", sheet_name="main")

# map each run's zip basename -> store object (dedup: take first hit)
by_base = {}
for _, r in live.iterrows():
    by_base.setdefault(r["base"], r)

# regex over the invocation's sample_def blocks
BLOCK = re.compile(
    r'"library_type":\s*"([^"]+)".*?"read_path":\s*"([^"]+)".*?"sample_names":\s*\[([^\]]*)\]',
    re.S,
)


def parse_invocation(txt):
    sid = re.search(r'sample_id\s*=\s*"([^"]+)"', txt)
    area = re.search(r'slide_serial_capture_area\s*=\s*"([^"]*)"', txt)
    rows = []
    for m in BLOCK.finditer(txt):
        rows.append({
            "sample_id": sid.group(1) if sid else "",
            "slide_area": area.group(1) if area else "",
            "library_type": m.group(1),
            "read_path": m.group(2),
            "sample_names": m.group(3).replace('"', "").strip(),
        })
    return rows


records = []
runs = [(r["experiment"], r["SB"], str(r["cellranger output file name"]))
        for _, r in main.iterrows()
        if str(r["cellranger output file name"]) in by_base]

print(f"{len(runs)} run zips to inspect")
for i, (exp, sb, base) in enumerate(runs, 1):
    obj = by_base[base]
    name = obj["name"]
    try:
        hp.download_from_project_store(store_name=STORE, file_name=name)
        cands = glob.glob(os.path.expanduser(f"~/input/**/{base}"), recursive=True)
        if not cands:
            raise FileNotFoundError("not found in FUSE mount")
        zp = cands[0]
        with zipfile.ZipFile(zp) as zf:
            invs = sorted((m for m in zf.namelist() if m.endswith("_invocation")), key=len)
            if not invs:
                raise ValueError("no _invocation member")
            txt = zf.read(invs[0]).decode("utf-8", "replace")
        rows = parse_invocation(txt)
        for row in rows:
            row.update({"experiment": exp, "SB": sb, "zip": base, "zip_uuid": obj["id"]})
            records.append(row)
        print(f"[{i}/{len(runs)}] {base}: {len(rows)} libs -> "
              f"{', '.join(sorted(set(r['sample_names'] for r in rows)))}")
    except Exception as e:
        print(f"[{i}/{len(runs)}] {base}: ERR {type(e).__name__} {str(e)[:120]}")
        records.append({"experiment": exp, "SB": sb, "zip": base,
                        "zip_uuid": obj["id"], "library_type": "", "read_path": "",
                        "sample_names": "", "sample_id": "", "slide_area": "",
                        "error": f"{type(e).__name__}: {str(e)[:120]}"})
    time.sleep(2)  # be gentle on the listing endpoint

df = pd.DataFrame(records)
cols = ["experiment", "SB", "zip", "sample_id", "slide_area", "library_type",
        "sample_names", "read_path", "zip_uuid"]
df = df.reindex(columns=cols + [c for c in df.columns if c not in cols])
df.to_csv(OUT, index=False)
print(f"\nwrote {len(df)} rows -> {OUT}")
