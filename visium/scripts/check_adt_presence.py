#!/usr/bin/env python3
"""
Cross-check every Visium sample in "visium data tracking.xlsx" (sheet `main`)
against its zip-tree listing under visium-zip-trees/ and mark which samples
have ADT (antibody-capture) data present.

WHY a single keyword is not enough
-----------------------------------
Space Ranger always emits an `analysis/` *skeleton* (graphclust + kmeans_2..10
for both gene_expression_* and antibody_capture_*). So merely finding the word
"antibody" marks every tree positive and would HIDE a real gene-expression-only
run. Instead we look for the artifacts that only exist when an antibody panel
was actually run and produced real output:

  * outs/feature_reference.csv .......... the antibody panel definition
  * analysis/.../antibody_capture_graphclust/clusters.csv  with real size
                                          (per-spot cluster calls; empty/absent
                                           on a GEX-only run)

Verdict per sample:
  ADT_PRESENT    feature_reference.csv AND a non-trivial antibody clusters.csv
  SCAFFOLD_ONLY  antibody_capture dirs exist but clusters.csv missing/placeholder
  GEX_ONLY       a real tree with NO antibody_capture analysis at all
  NO_TREE        no zip tree on disk for this sample (unprocessed / fastq-only)

Outputs crosscheck/adt_presence.csv and prints a summary. Read-only; no HISE.
"""

import os
import re
import glob
import shutil
import datetime
import pandas as pd
import openpyxl

BASE = "/home/workspace/PCD_data/visium/data"
XLSX = os.path.join(BASE, "visium data tracking.xlsx")
TREES = os.path.join(BASE, "visium-zip-trees")
OUT = os.path.join(BASE, "crosscheck")
os.makedirs(OUT, exist_ok=True)

# Header for the column written back into the xlsx `main` tab.
SR_COL = "Space Ranger output"

# A clusters.csv smaller than this is treated as a placeholder, not real calls.
MIN_CLUSTERS_KB = 5.0

SIZE_RE = re.compile(r"\(([0-9.]+)\s*([KMG]?)B\)")


def size_to_kb(token):
    """'(91.6KB)' / '(1.6MB)' -> float KB, or None."""
    if not token:
        return None
    m = SIZE_RE.search(token)
    if not m:
        return None
    val, unit = float(m.group(1)), m.group(2)
    return {"": 1 / 1024.0, "K": 1.0, "M": 1024.0, "G": 1024.0 * 1024}[unit] * val


def index_trees():
    """Map a zip filename -> tree .md path.

    Tree files are named '<zipname>.zip.md', sometimes with an
    'AIFI-<timestamp>__' prefix, so index both the full name and the part
    after '__'.
    """
    idx = {}
    for path in glob.glob(os.path.join(TREES, "*.zip.md")):
        zip_name = os.path.basename(path)[: -len(".md")]  # strip trailing .md
        idx.setdefault(zip_name, path)
        if "__" in zip_name:
            idx.setdefault(zip_name.split("__", 1)[1], path)
    return idx


def find_trees_for_cell(cell, idx):
    """Resolve the 'cellranger output file name' cell (may list several zips
    separated by ';') to tree paths. Returns (zip_names, tree_paths)."""
    if not isinstance(cell, str) or not cell.strip():
        return [], []
    names = [n.strip() for n in cell.split(";") if n.strip()]
    paths = []
    for n in names:
        p = idx.get(n)
        if p is None:
            # last-resort suffix match (handles stray prefixes)
            for k, v in idx.items():
                if k.endswith(n) or n.endswith(k):
                    p = v
                    break
        if p:
            paths.append(p)
    return names, paths


def assess_tree(path):
    """Parse one tree .md and return ADT evidence for it."""
    with open(path) as f:
        text = f.read()
    lines = text.splitlines()

    has_ab_dirs = "antibody_capture" in text
    has_ab_stage = ("SUMMARIZE_ANTIBODY_ANALYSIS" in text) or ("antibody_analysis/" in text)

    # feature_reference.csv (+ size of first occurrence)
    fref_kb = None
    has_fref = False
    for ln in lines:
        if "feature_reference.csv" in ln:
            has_fref = True
            fref_kb = size_to_kb(ln)
            break

    # antibody_capture_graphclust/clusters.csv size (real per-spot calls)
    clust_kb = None
    armed = False
    for ln in lines:
        if "antibody_capture_graphclust/" in ln:
            armed = True
            continue
        if armed and "clusters.csv" in ln:
            clust_kb = size_to_kb(ln)
            break

    real_clusters = clust_kb is not None and clust_kb >= MIN_CLUSTERS_KB

    if has_fref and real_clusters:
        verdict = "ADT_PRESENT"
    elif has_ab_dirs and not real_clusters:
        verdict = "SCAFFOLD_ONLY"
    elif not has_ab_dirs:
        verdict = "GEX_ONLY"
    else:
        verdict = "ADT_PRESENT"

    evidence = []
    if has_fref:
        evidence.append(f"feature_reference.csv({fref_kb:.1f}KB)" if fref_kb else "feature_reference.csv")
    if clust_kb is not None:
        evidence.append(f"antibody_graphclust_clusters.csv({clust_kb:.1f}KB)")
    if has_ab_stage:
        evidence.append("antibody_analysis_stage")
    if not evidence:
        evidence.append("no antibody artifacts")

    return {
        "verdict": verdict,
        "feature_reference_kb": round(fref_kb, 1) if fref_kb else "",
        "antibody_clusters_kb": round(clust_kb, 1) if clust_kb else "",
        "evidence": "; ".join(evidence),
    }


def sr_status_for_row(r, idx):
    """One Space Ranger-output status string per spreadsheet row.

    'Yes' when a Space Ranger output tree exists for the sample; otherwise
    'No — <what is available instead>'.
    """
    cell = r.get("cellranger output file name")
    _, paths = find_trees_for_cell(cell, idx)
    if paths:
        return "Yes"
    fastq = r.get("fastq filename")
    if isinstance(cell, str) and "cytassist" in cell.lower():
        return "No — raw CytAssist zip only (not Space Ranger outs)"
    if isinstance(fastq, str) and fastq.strip():
        return "No — FASTQ only"
    return "No — not yet processed / no spatial output"


def update_xlsx(df, idx):
    """Write the SR_COL column into the `main` tab in row order, preserving
    all other sheets/formatting. Idempotent: reuses the column if it exists.
    Backs up the workbook first."""
    stamp = datetime.date.today().strftime("%Y%m%d")
    backup = os.path.join(BASE, f"visium data tracking.backup-{stamp}.xlsx")
    if not os.path.exists(backup):
        shutil.copy2(XLSX, backup)
        print(f"Backed up workbook -> {os.path.basename(backup)}")

    wb = openpyxl.load_workbook(XLSX)
    ws = wb["main"]
    header = [ws.cell(row=1, column=c).value for c in range(1, ws.max_column + 1)]

    if SR_COL in header:
        col = header.index(SR_COL) + 1
    elif None in header:
        col = header.index(None) + 1          # first empty header cell
    else:
        col = ws.max_column + 1
    ws.cell(row=1, column=col, value=SR_COL)

    for i, (_, r) in enumerate(df.iterrows()):
        ws.cell(row=2 + i, column=col, value=sr_status_for_row(r, idx))

    wb.save(XLSX)
    print(f"Wrote column '{SR_COL}' (col {col}) into main tab for {len(df)} rows.")


def main():
    idx = index_trees()
    df = pd.read_excel(XLSX, sheet_name="main")

    rows = []
    for _, r in df.iterrows():
        cell = r.get("cellranger output file name")
        names, paths = find_trees_for_cell(cell, idx)

        if not paths:
            rows.append({
                "id": r["id"], "experiment": r["experiment"], "SB": r.get("SB"),
                "project store": r.get("project store"),
                "cellranger_output_file_name": cell if isinstance(cell, str) else "",
                "tree_file": "", "tree_found": "N", "adt_verdict": "NO_TREE",
                "feature_reference_kb": "", "antibody_clusters_kb": "",
                "evidence": "no zip tree on disk (unprocessed / fastq-only / raw CytAssist)",
            })
            continue

        # A row may map to several zips; assess each.
        for name, path in zip(names + [None] * (len(paths) - len(names)), paths):
            a = assess_tree(path)
            rows.append({
                "id": r["id"], "experiment": r["experiment"], "SB": r.get("SB"),
                "project store": r.get("project store"),
                "cellranger_output_file_name": name or os.path.basename(path)[:-3],
                "tree_file": os.path.basename(path), "tree_found": "Y",
                "adt_verdict": a["verdict"],
                "feature_reference_kb": a["feature_reference_kb"],
                "antibody_clusters_kb": a["antibody_clusters_kb"],
                "evidence": a["evidence"],
            })

    out = pd.DataFrame(rows)
    out_path = os.path.join(OUT, "adt_presence.csv")
    out.to_csv(out_path, index=False)

    print(f"Wrote {out_path}  ({len(out)} rows)\n")
    print("=== verdict counts ===")
    print(out["adt_verdict"].value_counts().to_string())
    print("\n=== per-sample ===")
    show = out[["id", "experiment", "adt_verdict", "antibody_clusters_kb",
                "feature_reference_kb", "tree_file"]]
    print(show.to_string(index=False))

    no_tree = out[out["adt_verdict"] == "NO_TREE"]
    if not no_tree.empty:
        print("\n=== samples with NO processed tree (cannot confirm ADT) ===")
        print(no_tree[["id", "experiment", "project store", "cellranger_output_file_name"]].to_string(index=False))

    # Write the Space Ranger-output status back into the xlsx main tab.
    print()
    update_xlsx(df, idx)


if __name__ == "__main__":
    main()
