#!/usr/bin/env python
"""Survey the internal structure of Xenium output bundles, one per analyzer version.

Every Xenium output bundle is a zip in the PCD_Spatial project store, referenced
by the `hise uuid` on the tracking sheet's `main` tab. Inside a bundle the file
names are generic (cell_feature_matrix.h5, transcripts.parquet, ...) — the sample
identity lives only in the outer zip name. So every bundle produced by the same
Xenium Analyzer run shares an identical internal file list.

Runs are grouped by the YYYYMMDD__HHMMSS timestamp token at the end of the zip
name: same timestamp -> same run -> same analyzer version -> same file list. There
are 20 unique timestamps across the 97 rows, so we download one representative
bundle per timestamp (20 downloads) and that fully characterises the file
convention for all samples.

For each representative we:
    download zip by UUID -> read its internal structure + analyzer version
    -> write a Markdown tree file -> delete the zip -> next.
Peak disk usage is a single bundle. hisepy only offers a full download (no
byte-range read), so each zip is fully fetched once before we can list it.

Output: one Markdown tree per timestamp under TREE_DIR, named
`<timestamp>__<version>.md`.

Run from the repo root with the HISE-authenticated environment:
    /home/workspace/environment/scanpyv3/bin/python xenium/scripts/xenium_zip_tree.py
"""

import json
import os
import re
import sys
import unittest.mock as mock
import zipfile
from collections import defaultdict

import pandas as pd

# hisepy builds a DataFrame from the store listing; this pandas build chokes on
# Arrow-backed string columns, so force plain-Python strings before importing.
pd.set_option("future.infer_string", False)

# Two hisepy submodules fail to import in this env and aren't needed here.
with mock.patch.dict("sys.modules", {"hisepy.training_job": mock.MagicMock(),
                                     "hisepy.ai_ray_transformer": mock.MagicMock()}):
    import hisepy as hp

import openpyxl

STORE = "PCD_Spatial"
SHEET = "xenium/data/xenium data tracking.xlsx"
TREE_DIR = "xenium/data/xenium-zip-trees"

# experiment.xenium key names for the analyzer version have drifted across
# versions; check these in order and use the first present.
VERSION_KEYS = ["analysis_sw_version", "software_version", "xenium_analyzer_version"]

TS_RE = re.compile(r"__(\d{8}__\d{6})\.zip$")


def human_size(num_bytes):
    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{int(size)}B" if unit == "B" else f"{size:.1f}{unit}"
        size /= 1024


def build_tree(infolist):
    """Build an ASCII tree (list of lines) from a zip's ZipInfo list."""
    root = {}
    for info in infolist:
        parts = [p for p in info.filename.split("/") if p]
        if not parts:
            continue
        node = root
        for part in parts[:-1]:
            node = node.setdefault(part, {})
        last = parts[-1]
        if info.is_dir():
            node.setdefault(last, {})
        else:
            node.setdefault(last, {})["__size__"] = info.file_size

    lines = []

    def walk(node, prefix):
        entries = sorted((k for k in node if k != "__size__"),
                         key=lambda k: ("__size__" in node[k], k.lower()))
        for i, name in enumerate(entries):
            last = i == len(entries) - 1
            connector = "└── " if last else "├── "
            child = node[name]
            if "__size__" in child:
                lines.append(f"{prefix}{connector}{name}  ({human_size(child['__size__'])})")
            else:
                lines.append(f"{prefix}{connector}{name}/")
                walk(child, prefix + ("    " if last else "│   "))

    walk(root, "")
    return lines


def read_version(zf):
    names = [n for n in zf.namelist() if n.endswith("experiment.xenium")]
    if not names:
        return "unknown"
    with zf.open(min(names, key=len)) as fh:
        data = json.load(fh)
    for key in VERSION_KEYS:
        if key in data:
            return str(data[key])
    return "unknown"


def representatives():
    """One (timestamp, sample_id, store_name, uuid) per unique run timestamp."""
    live = hp.list_files_in_project_store(STORE)
    id2name = dict(zip(live["id"], live["name"]))

    ws = openpyxl.load_workbook(SHEET, data_only=True)["main"]
    hdr = {c.value: i for i, c in enumerate(ws[1]) if c.value}
    iu, isid = hdr["hise uuid"], hdr["sample id"]

    groups = defaultdict(list)
    for row in ws.iter_rows(min_row=2, values_only=True):
        uuid = row[iu]
        if not uuid:
            continue
        name = id2name.get(uuid)
        m = TS_RE.search(name) if name else None
        if m:
            groups[m.group(1)].append((row[isid], name, uuid))

    reps = []
    for ts in sorted(groups):
        sid, name, uuid = sorted(groups[ts])[0]
        reps.append((ts, sid, name, uuid))
    return reps


def write_tree_file(ts, sample_id, store_name, version, instrument, infolist):
    n_files = sum(1 for i in infolist if not i.is_dir())
    n_dirs = sum(1 for i in infolist if i.is_dir())
    total = sum(i.file_size for i in infolist if not i.is_dir())
    dst = os.path.join(TREE_DIR, f"{ts}__{version}.md")
    with open(dst, "w") as fh:
        fh.write(f"# Xenium run {ts}\n\n")
        fh.write(f"- Analyzer version: `{version}`\n")
        fh.write(f"- Instrument: `{instrument}`\n")
        fh.write(f"- Representative sample: `{sample_id}`\n")
        fh.write(f"- Store file: `{store_name}`\n")
        fh.write(f"- Entries: {n_files} files, {n_dirs} directories\n")
        fh.write(f"- Uncompressed total: {human_size(total)}\n\n")
        fh.write("```\n")
        fh.write("\n".join(build_tree(infolist)))
        fh.write("\n```\n")
    return dst


def main():
    os.makedirs(TREE_DIR, exist_ok=True)
    reps = representatives()
    print(f"{len(reps)} unique run timestamps to survey.\n")

    ok, failed = 0, []
    for i, (ts, sid, name, uuid) in enumerate(reps, 1):
        instr = name.split("output-")[-1].split("__")[0]
        print(f"[{i}/{len(reps)}] {ts}  {sid}  ({instr})")
        zip_path = None
        try:
            paths = hp.cache_files(file_ids=[uuid])
            if not paths:
                print("  cache_files returned no paths"); failed.append(ts); continue
            zip_path = paths[0]
            with zipfile.ZipFile(zip_path) as zf:
                version = read_version(zf)
                infolist = zf.infolist()
            dst = write_tree_file(ts, sid, name, version, instr, infolist)
            print(f"  version {version} -> {dst}")
            ok += 1
        except Exception as err:
            print(f"  ERROR: {err}")
            failed.append(ts)
        finally:
            try:
                if zip_path and os.path.exists(zip_path):
                    os.remove(zip_path)
            except OSError as err:
                print(f"  WARNING: could not delete zip: {err}")

    print(f"\nDone. {ok} trees in {TREE_DIR}, {len(failed)} failed.")
    for ts in failed:
        print(f"  failed: {ts}")


if __name__ == "__main__":
    main()
