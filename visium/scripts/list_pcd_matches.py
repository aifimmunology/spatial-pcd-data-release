#!/usr/bin/env python3
"""List files in PCD_Spatial matching either sample's TIS/EXP id (or variants)."""
import pandas as pd, re

CSV = "/home/workspace/PCD_data/visium/data/project-stores/PCD_Spatial_contents.csv"

# samples to search: (TIS id, EXP id)
SAMPLES = [
    ("TIS05396-001-001", "EXP-00872"),
    ("TIS05397-001-001", "EXP-00876"),
]

def variants(tis, exp):
    """All formatting variants to match for a TIS id and EXP id."""
    td = re.match(r"TIS(\d+)", tis).group(1)   # e.g. 05396
    ed = re.search(r"(\d+)", exp).group(1)      # e.g. 00872
    return [
        tis, tis.replace("-", "_"), tis.replace("-", ""), "TIS" + td, td,  # TIS variants
        exp, "EXP" + ed, "EXP_" + ed, ed,                                  # EXP variants
    ]

df = pd.read_csv(CSV, dtype=str).fillna("")
for tis, exp in SAMPLES:
    pats = [v.lower() for v in variants(tis, exp)]
    hits = df[df["name"].str.lower().apply(lambda n: any(p in n for p in pats))]
    print(f"\n=== {tis} / {exp} : {len(hits)} file(s) ===")
    for name in hits["name"]:
        print("  ", name)
