#!/usr/bin/env python3
"""
Resolve differences between the GEO deliverable template and the actual Space
Ranger output, ahead of reorganizing the 31 processed runs.

`reference/final_visium_structure_example.txt` defines the target per-sample
layout (gene_expression raw/processed/spatial, ADT_conventional, ADT_CST_bridge,
web_summary, H&E).  Before we can transform real runs into that layout we need to
know, for every run:

  1. Which template-required SOURCE files are MISSING from `outs/`
     (e.g. a run that lacks molecule_info.h5, or uses the old
     `tissue_positions_list.csv` name instead of `tissue_positions.csv`).
  2. Which files the actual `outs/` carries that the template does NOT mention
     but arguably SHOULD carry into the deliverable (e.g. cytassist_image.tiff,
     metrics_summary.csv).

This is LOCAL ONLY: it parses the per-zip content trees already captured in
data/visium-zip-trees/*.md (no HISE).  The set of 31 processed runs and the
tree file for each comes from crosscheck/adt_presence.csv (tree_found == 'Y').

We only inspect the canonical `<run>/outs/` directory and its `outs/spatial/`
subdirectory -- the two levels the template draws from.  Deeper analysis/ and
pipeline-intermediate trees are intentionally ignored.

Outputs (data/crosscheck/):
  template_vs_output_matrix.csv  one row per run x expected source file (present?)
  template_vs_output_extras.csv  files in outs/ + outs/spatial/ not in template
  template_vs_output_report.txt  human-readable summary of the differences
"""

import glob
import os
import re

import pandas as pd

BASE = "/home/workspace/PCD_data/visium/data"
TREE_DIR = os.path.join(BASE, "visium-zip-trees")
ADT_PRESENCE = os.path.join(BASE, "crosscheck", "adt_presence.csv")
OUT_DIR = os.path.join(BASE, "crosscheck")

# ---------------------------------------------------------------------------
# What the template (reference/final_visium_structure_example.txt) requires.
# Keyed by template bucket -> the SOURCE file(s) that must exist in Space Ranger
# `outs/` to build it.  ADT_conventional is DERIVED from the same combined
# matrices (split on feature_type) so its source is the GEX matrices + the
# feature_reference.csv that proves an antibody panel exists.  ADT_CST_bridge and
# H&E have no in-tree source (see README key findings) and are not checked here.
# ---------------------------------------------------------------------------
# expected file basename -> (where it lives: "outs" or "spatial", template bucket)
EXPECTED = {
    "raw_feature_bc_matrix.h5":      ("outs",    "gene_expression/raw"),
    "filtered_feature_bc_matrix.h5": ("outs",    "gene_expression/processed"),
    "molecule_info.h5":              ("outs",    "gene_expression/processed"),
    "web_summary.html":              ("outs",    "web_summary.html"),
    # feature_reference.csv: not a literal deliverable file, but required to
    # split out ADT_conventional; track it as a source dependency.
    "feature_reference.csv":         ("outs",    "ADT_conventional (split source)"),
    "tissue_hires_image.png":        ("spatial", "gene_expression/spatial"),
    "tissue_lowres_image.png":       ("spatial", "gene_expression/spatial"),
    "tissue_positions.csv":          ("spatial", "gene_expression/spatial"),
    "scalefactors_json.json":        ("spatial", "gene_expression/spatial"),
    "aligned_fiducials.jpg":         ("spatial", "gene_expression/spatial"),
}

# Accepted alternate names (older Space Ranger). Maps the canonical template
# name -> set of names that satisfy it.
ALIASES = {
    "tissue_positions.csv": {"tissue_positions.csv", "tissue_positions_list.csv"},
    # Two BM runs ship a sample-suffixed web summary instead of the canonical name.
    "web_summary.html": {"web_summary.html", "web_summary_exp897.html",
                         "web_summary_exp898.html"},
}

CONNECTOR = re.compile(r"^(?P<indent>(?:[│ ]   )*)(?:├──|└──) (?P<name>.+)$")
SIZE_SUFFIX = re.compile(r"  \([0-9.]+[KMGT]?B\)\s*$")


def parse_outs(md_path):
    """Return (outs_children, spatial_children) basename sets for a tree md.

    `outs/` is NOT at a fixed depth -- some zips wrap it under several wrapper
    directories.  So we reconstruct the full path of every entry via a depth
    stack, collect the direct children of every directory whose path ends in
    `/outs`, and pick the canonical `outs` (the one that best matches the
    expected Space Ranger files).  Its `spatial/` direct children are returned
    too.  Deeper garbled / duplicate subtrees thus get ignored automatically.
    """
    expected_outs = {n for n, (lvl, _) in EXPECTED.items() if lvl == "outs"}
    outs_candidates = {}     # outs_path -> set(child basenames)
    spatial_candidates = {}  # outs_path -> set(spatial child basenames)
    stack = {}  # depth -> name
    in_fence = False

    with open(md_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.strip() == "```":
                in_fence = not in_fence
                continue
            if not in_fence:
                continue
            m = CONNECTOR.match(line)
            if not m:
                continue
            depth = len(m.group("indent")) // 4
            name = SIZE_SUFFIX.sub("", m.group("name")).rstrip()
            clean = name.rstrip("/")
            stack[depth] = clean
            for d in list(stack):  # prune deeper entries
                if d > depth:
                    del stack[d]
            if depth == 0:
                continue
            parent_path = "/".join(stack[d] for d in range(depth))
            if parent_path.endswith("/outs") or parent_path == "outs":
                outs_candidates.setdefault(parent_path, set()).add(clean)
            elif parent_path.endswith("/outs/spatial") or parent_path == "outs/spatial":
                outs_path = parent_path[: -len("/spatial")]
                spatial_candidates.setdefault(outs_path, set()).add(clean)

    if not outs_candidates:
        return set(), set()
    # canonical outs = best overlap with expected files, then most children
    best = max(outs_candidates,
               key=lambda p: (len(outs_candidates[p] & expected_outs),
                              len(outs_candidates[p])))
    return outs_candidates[best], spatial_candidates.get(best, set())


def satisfied(canonical, present_set):
    """Is the canonical expected file present, allowing known aliases?"""
    names = ALIASES.get(canonical, {canonical})
    hit = names & present_set
    return (bool(hit), next(iter(hit)) if hit else "")


def main():
    adt = pd.read_csv(ADT_PRESENCE)
    runs = adt[adt["tree_found"] == "Y"].copy()
    # dedup by tree_file (some ids share a run / appear twice)
    runs = runs.drop_duplicates(subset=["tree_file"]).reset_index(drop=True)

    matrix_rows = []
    extras_rows = []
    # canonical set of names the template knows about, per level
    template_outs = {n for n, (lvl, _) in EXPECTED.items() if lvl == "outs"}
    template_spatial = {n for n, (lvl, _) in EXPECTED.items() if lvl == "spatial"}
    # extras we expect to see but are deliberately non-deliverable noise
    alias_names = {a for s in ALIASES.values() for a in s}

    missing_count = {}  # tree_file -> list of missing canonical files
    for _, r in runs.iterrows():
        tf = str(r["tree_file"])
        md_path = os.path.join(TREE_DIR, tf)
        if not os.path.exists(md_path):
            print(f"  !! tree md not found: {tf}")
            continue
        outs_c, spatial_c = parse_outs(md_path)
        run_id = f"{r['id']} / {r['experiment']}"
        missing = []
        for canonical, (lvl, bucket) in EXPECTED.items():
            present_set = outs_c if lvl == "outs" else spatial_c
            ok, found_as = satisfied(canonical, present_set)
            matrix_rows.append({
                "id": r["id"], "experiment": r["experiment"], "tree_file": tf,
                "level": lvl, "template_bucket": bucket,
                "expected_file": canonical, "present": "Y" if ok else "N",
                "found_as": found_as if found_as != canonical else "",
            })
            if not ok:
                missing.append(f"{lvl}/{canonical}")
        if missing:
            missing_count[tf] = missing

        # extras: in outs/ or outs/spatial/ but not part of the template
        for nm in sorted(outs_c - template_outs - alias_names):
            extras_rows.append({"id": r["id"], "experiment": r["experiment"],
                                "tree_file": tf, "level": "outs", "extra": nm})
        for nm in sorted(spatial_c - template_spatial - alias_names):
            extras_rows.append({"id": r["id"], "experiment": r["experiment"],
                                "tree_file": tf, "level": "spatial", "extra": nm})

    mdf = pd.DataFrame(matrix_rows)
    edf = pd.DataFrame(extras_rows)
    os.makedirs(OUT_DIR, exist_ok=True)
    mdf.to_csv(os.path.join(OUT_DIR, "template_vs_output_matrix.csv"), index=False)
    edf.to_csv(os.path.join(OUT_DIR, "template_vs_output_extras.csv"), index=False)

    # ---- report ----
    n_runs = runs["tree_file"].nunique()
    lines = [f"Template vs Space Ranger output -- {n_runs} processed runs", ""]

    lines.append("== 1. Template-required source files: presence across runs ==")
    for canonical, (lvl, bucket) in EXPECTED.items():
        sub = mdf[mdf["expected_file"] == canonical]
        n_present = (sub["present"] == "Y").sum()
        alt = sub[(sub["present"] == "Y") & (sub["found_as"] != "")]["found_as"].unique()
        note = f"  (some via alias: {', '.join(alt)})" if len(alt) else ""
        flag = "" if n_present == n_runs else "   <-- NOT IN ALL RUNS"
        lines.append(f"  {lvl:8} {canonical:32} {n_present}/{n_runs} present{note}{flag}")
    lines.append("")

    lines.append("== 2. Runs missing one or more required source files ==")
    if missing_count:
        for tf, miss in sorted(missing_count.items()):
            lines.append(f"  {tf}")
            for m in miss:
                lines.append(f"      missing: {m}")
    else:
        lines.append("  (none -- every run has every template-required source file)")
    lines.append("")

    lines.append("== 3. Files present in outs/ or outs/spatial/ but NOT in the template ==")
    lines.append("   (candidates the deliverable layout may want to carry / map)")
    if not edf.empty:
        tally = (edf.groupby(["level", "extra"]).size()
                 .reset_index(name="n_runs").sort_values(["level", "n_runs"], ascending=[True, False]))
        for _, t in tally.iterrows():
            lines.append(f"  {t['level']:8} {t['extra']:36} in {t['n_runs']}/{n_runs} runs")
    else:
        lines.append("  (none)")
    lines.append("")

    report = os.path.join(OUT_DIR, "template_vs_output_report.txt")
    with open(report, "w") as f:
        f.write("\n".join(lines) + "\n")

    print("\n".join(lines))
    print(f"\nWrote:\n  {os.path.join(OUT_DIR, 'template_vs_output_matrix.csv')}"
          f"\n  {os.path.join(OUT_DIR, 'template_vs_output_extras.csv')}"
          f"\n  {report}")


if __name__ == "__main__":
    main()
