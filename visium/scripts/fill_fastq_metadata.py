"""Fill FASTQ location metadata + HISE UUIDs into the Visium tracking sheet.

Populates `fastq filename` and `fastq UUID` on the `main` tab, appends FASTQ data-flags
to `Notes`, and (re)builds a `FASTQ files` tab (one row per run x library).

Authoritative provenance comes from each run's Space Ranger `_invocation`
(visium/data/crosscheck/fastq_invocations.csv): per library, the FASTQ `sample_names`
prefix and the compute-node `read_path` the reads were consumed from. Cross-checked
against all 42 HISE project stores + the file ledger (see fastq_hise_locations*): the
ONLY per-file FASTQs archived in HISE are SB002's (flowcell HNKWTDMXY) -- EXP-00872 and
EXP-00875 -- in the PCD_Spatial store; their per-file UUIDs come from the live listing.
Every other batch's raw FASTQs were never ingested to HISE (only the processed output
zip is), so there is no FASTQ UUID to record -- only the provenance path.

Run from repo root: python visium/scripts/fill_fastq_metadata.py
"""
import re
import shutil
from datetime import datetime

import openpyxl
import pandas as pd

XLSX = "visium/data/visium data tracking.xlsx"
INV = "visium/data/crosscheck/fastq_invocations.csv"
LIVE = "visium/data/project-stores/PCD_Spatial_live_2026-07-13.csv"
PERFILE_OUT = "visium/data/crosscheck/fastq_sb002_perfile.csv"

HISE_FLOWCELL = "HNKWTDMXY"  # the only flowcell whose FASTQs are archived in HISE
LIB_SHORT = {"Gene Expression": "GEX", "Antibody Capture": "ADT"}

UUID_NOTE_HISE = "see 'FASTQ files' tab — 16 files (8 GEX + 8 ADT) in PCD_Spatial"
UUID_NOTE_ABSENT = "not in HISE — raw FASTQs never ingested (only processed output zip archived)"
ABSENT_INTERNAL = "n/a — AIFI Internal (dropped from release)"

# FASTQ data-flags (keyed by tracking-sheet id) -> appended to main `Notes` + tab notes
FLAGS = {
    "TIS06392-002": "FASTQ FLAG: FASTQs deliberately switched (sample_id "
                    "'TIS06392-002_modified_switched_fastq') — GEX sourced from TIS06377-001, "
                    "ADT from TIS06392-002.",
    "TIS06377-001": "FASTQ FLAG: FASTQs deliberately switched (sample_id "
                    "'TIS06377-001_modified_switched_fastq') — GEX sourced from TIS06392-002, "
                    "ADT from TIS06377-001.",
    "TMA00147": "FASTQ FLAG: output zip carries TMA00144's sample_id + area (V53S25-332-A1) and "
                "uses TMA00144 FASTQs — TMA00147 output appears to duplicate/mislabel TMA00144.",
    "TMA00144": "FASTQ FLAG: the TMA00147 output zip duplicates this run's sample_id/FASTQs "
                "(see TMA00147).",
}


def clean_source(read_path, sb):
    s = str(read_path).replace("/home/jupyter/Spatial/", "").replace("/home/jupyter/", "")
    if s and not s.startswith("SB") and isinstance(sb, str) and sb.startswith("SB"):
        s = f"{sb}/{s}"
    return s


# ---- per-file SB002 UUIDs (prefix -> [(filename, uuid), ...]) --------------------
live = pd.read_csv(LIVE)
live_fq = live[live["name"].str.endswith(".fastq.gz", na=False)].copy()
live_fq["fn"] = live_fq["name"].str.split("/").str[-1]
live_fq["prefix"] = live_fq["fn"].str.replace(r"_S\d+_L\d+_[IR][12]_001\.fastq\.gz$", "", regex=True)

perfile_rows = []
uuids_by_prefix = {}
for pref in ["EXP-00872_G", "EXP-00872_P", "EXP-00875_G", "EXP-00875_P"]:
    sub = live_fq[live_fq["prefix"] == pref].sort_values("fn")
    uuids_by_prefix[pref] = list(sub["id"])
    for _, r in sub.iterrows():
        perfile_rows.append({"prefix": pref, "filename": r["fn"], "hise_uuid": r["id"],
                             "store": "PCD_Spatial", "full_path": r["name"]})
pd.DataFrame(perfile_rows).to_csv(PERFILE_OUT, index=False)
assert all(len(v) == 8 for v in uuids_by_prefix.values()), uuids_by_prefix

# ---- per run x library provenance from invocations -------------------------------
inv = pd.read_csv(INV)
# EXP-00872 has no _invocation (reprocessed run); supply its known provenance.
EXP00872 = {
    "Gene Expression": ("EXP-00872_G", "SB002/FC00516_done/outs/fastq_path/HNKWTDMXY/"),
    "Antibody Capture": ("EXP-00872_P", "SB002/FC00516_done/outs/fastq_path/HNKWTDMXY/"),
}

# build: zip -> {library_type: (prefix, source, flowcell)}
by_zip = {}
for _, r in inv.iterrows():
    if pd.isna(r["sample_names"]) or not str(r["sample_names"]).strip():
        continue
    pref = str(r["sample_names"]).replace('"', "").strip()
    src = clean_source(r["read_path"], r["SB"])
    flow = re.search(r"(H[A-Z0-9]{8}|HNKWTDMXY)", str(r["read_path"]))
    by_zip.setdefault(r["zip"], {})[r["library_type"]] = (pref, src, flow.group(0) if flow else "")


def libs_for(zip_name):
    if zip_name == "EXP00872_TIS05396_A1.zip":
        return {lt: (p, s, HISE_FLOWCELL) for lt, (p, s) in EXP00872.items()}
    return by_zip.get(zip_name, {})


# ---- back up then edit -----------------------------------------------------------
bak = f"{XLSX}.bak-{datetime.now():%Y%m%d-%H%M%S}"
shutil.copy2(XLSX, bak)
print("backup:", bak)

wb = openpyxl.load_workbook(XLSX)
ws = wb["main"]
H = {c.value: c.column for c in ws[1]}
c_id, c_exp, c_env, c_sb = H["id"], H["experiment"], H["HISE environment"], H["SB"]
c_zip, c_fn, c_uuid, c_notes = (H["cellranger output file name"], H["fastq filename"],
                                H["fastq UUID"], H["Notes"])

fastq_tab_rows = []
updated = 0
for row in range(2, ws.max_row + 1):
    sid = ws.cell(row, c_id).value
    if not sid:
        continue
    updated += 1
    env = str(ws.cell(row, c_env).value or "")
    sb = ws.cell(row, c_exp).value and ws.cell(row, c_sb).value
    zip_name = str(ws.cell(row, c_zip).value)
    exp = ws.cell(row, c_exp).value
    libs = libs_for(zip_name)

    if "Internal" in env:  # SB009, dropped
        ws.cell(row, c_fn).value = ABSENT_INTERNAL
        ws.cell(row, c_uuid).value = ABSENT_INTERNAL
    elif libs:
        gex = libs.get("Gene Expression")
        adt = libs.get("Antibody Capture")
        in_hise = any(v[2] == HISE_FLOWCELL for v in libs.values())
        src = (gex or adt)[1]
        parts = []
        if gex:
            parts.append(f"{gex[0]} (GEX)")
        if adt:
            parts.append(f"{adt[0]} (ADT)")
        ws.cell(row, c_fn).value = f"{' + '.join(parts)} — {src}"
        ws.cell(row, c_uuid).value = UUID_NOTE_HISE if in_hise else UUID_NOTE_ABSENT
        # FASTQ files tab rows (one per library)
        for lt in ("Gene Expression", "Antibody Capture"):
            if lt not in libs:
                continue
            pref, lsrc, flow = libs[lt]
            hit = flow == HISE_FLOWCELL
            uu = "; ".join(uuids_by_prefix.get(pref, [])) if hit else "N/A — not in HISE"
            fastq_tab_rows.append([sid, exp, ws.cell(row, c_sb).value, lt, pref, lsrc,
                                   "Y" if hit else "N",
                                   "PCD_Spatial" if hit else "",
                                   uu, FLAGS.get(str(sid), "")])
    else:
        ws.cell(row, c_fn).value = "unresolved (no _invocation)"
        ws.cell(row, c_uuid).value = UUID_NOTE_ABSENT

    # append FASTQ flag to Notes without clobbering
    if str(sid) in FLAGS:
        cur = ws.cell(row, c_notes).value
        flag = FLAGS[str(sid)]
        ws.cell(row, c_notes).value = f"{cur} | {flag}" if cur and str(cur).strip() else flag

# ---- (re)build FASTQ files tab ---------------------------------------------------
if "FASTQ files" in wb.sheetnames:
    del wb["FASTQ files"]
fs = wb.create_sheet("FASTQ files")
cols = ["sample id", "experiment", "SB", "library", "fastq sample prefix",
        "source (compute-node read_path)", "in_HISE", "store", "HISE UUIDs", "notes"]
fs.append(cols)
for r in fastq_tab_rows:
    fs.append(r)

wb.save(XLSX)
print(f"updated {updated} main rows; wrote {len(fastq_tab_rows)} rows to 'FASTQ files'")
print(f"per-file SB002 UUIDs -> {PERFILE_OUT} ({len(perfile_rows)} files)")
