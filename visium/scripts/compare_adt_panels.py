#!/usr/bin/env python3
"""
Compare the ADT panels assembled by extract_feature_references.py.

Reads every `*__feature_reference.csv` in data/adt-panels/, normalizes each to
its antibody-marker set, and groups the runs by identical panel.  Answers the
question the README's "Key finding" raises: does every processed run share the
same Antibody Capture panel, or does any run differ?

A panel is identified by the set of antibody rows it defines.  Two runs are the
"same panel" iff they define the same set of (id, name, sequence, feature_type)
tuples — column order, row order, and any pipeline-added formatting are ignored.

Outputs
-------
  data/adt-panels/_panel_comparison.csv   one row per run: file, marker counts,
                                          panel group id, isotype-control count
  data/adt-panels/_panel_groups.txt       human-readable: each distinct panel,
                                          its members, and pairwise diffs

Local-only; no HISE.  Run after extract_feature_references.py.
"""

import glob
import hashlib
import os

import pandas as pd

BASE = "/home/workspace/PCD_data/visium/data"
PANEL_DIR = os.path.join(BASE, "adt-panels")

# 10x feature_reference.csv key columns; "sequence" is the barcode that pins a
# marker, "feature_type" separates Antibody Capture from any Gene Expression rows.
KEY_COLS = ["id", "name", "sequence", "feature_type"]
ISOTYPE_HINTS = ("igg", "isotype", "control")


def load_panel(path):
    df = pd.read_csv(path)
    df.columns = [c.strip().lower() for c in df.columns]
    # Keep only antibody-capture rows when feature_type is present.
    if "feature_type" in df.columns:
        ab = df[df["feature_type"].str.contains("antibody", case=False, na=False)]
        if not ab.empty:
            df = ab
    cols = [c for c in KEY_COLS if c in df.columns]
    rows = set()
    for _, r in df.iterrows():
        rows.add(tuple(str(r.get(c, "")).strip() for c in cols))
    return df, rows, cols


def panel_hash(rows):
    canon = "\n".join(sorted("\t".join(t) for t in rows))
    return hashlib.sha1(canon.encode()).hexdigest()[:10]


def is_isotype(row_tuple, cols):
    text = " ".join(row_tuple).lower()
    return any(h in text for h in ISOTYPE_HINTS)


def main():
    paths = sorted(glob.glob(os.path.join(PANEL_DIR, "*__feature_reference.csv")))
    if not paths:
        print(f"No panels found in {PANEL_DIR}. Run extract_feature_references.py first.")
        return

    panels = {}   # path -> (rows, cols)
    summary = []
    groups = {}   # hash -> list of paths
    for p in paths:
        _, rows, cols = load_panel(p)
        panels[p] = (rows, cols)
        h = panel_hash(rows)
        groups.setdefault(h, []).append(p)
        n_iso = sum(1 for t in rows if is_isotype(t, cols))
        summary.append({
            "file": os.path.basename(p),
            "n_markers": len(rows),
            "n_isotype_controls": n_iso,
            "n_real_markers": len(rows) - n_iso,
            "panel_group": h,
        })

    sdf = pd.DataFrame(summary).sort_values(["panel_group", "file"])
    out_csv = os.path.join(PANEL_DIR, "_panel_comparison.csv")
    sdf.to_csv(out_csv, index=False)

    print(f"{len(paths)} panels across {len(groups)} distinct panel group(s).\n")
    print(sdf.to_string(index=False))

    # Human-readable group report with pairwise diffs.
    lines = [f"{len(paths)} panels; {len(groups)} distinct panel group(s).", ""]
    ordered = sorted(groups.items(), key=lambda kv: -len(kv[1]))
    for gi, (h, members) in enumerate(ordered, 1):
        rows0, cols0 = panels[members[0]]
        lines.append(f"=== Panel group {gi}  [{h}]  ({len(members)} run(s), {len(rows0)} markers) ===")
        for m in members:
            lines.append(f"   {os.path.basename(m)}")
        lines.append("")

    if len(ordered) > 1:
        lines.append("=== Pairwise differences between groups ===")
        ref_h, ref_members = ordered[0]
        ref_rows = panels[ref_members[0]][0]
        ref_ids = {t[0] for t in ref_rows}
        for gi, (h, members) in enumerate(ordered[1:], 2):
            rows, _ = panels[members[0]]
            ids = {t[0] for t in rows}
            only_ref = sorted(ref_ids - ids)
            only_grp = sorted(ids - ref_ids)
            lines.append(f"\nGroup {gi} [{h}] vs Group 1 [{ref_h}] (by marker id):")
            lines.append(f"   only in Group 1: {only_ref or '(none)'}")
            lines.append(f"   only in Group {gi}: {only_grp or '(none)'}")
    else:
        lines.append("All runs share an IDENTICAL antibody panel.")

    out_txt = os.path.join(PANEL_DIR, "_panel_groups.txt")
    with open(out_txt, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"\nWrote {out_csv}")
    print(f"Wrote {out_txt}")
    if len(groups) == 1:
        print("\n==> All processed runs share an IDENTICAL antibody panel.")
    else:
        print(f"\n==> {len(groups)} DISTINCT panels detected — see _panel_groups.txt for diffs.")


if __name__ == "__main__":
    main()
