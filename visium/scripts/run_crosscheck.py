#!/usr/bin/env python3
"""
Cross-check the 6 "unk" Visium rows against all 33 project-store contents.
LOCAL ONLY -- reads the *_contents.csv files already on disk + _store_index.csv.
No HISE calls. Writes 3 paper-trail CSVs under crosscheck/.
"""
import pandas as pd, glob, os, re

BASE = "/home/workspace/PCD_data/visium/data"
PS   = os.path.join(BASE, "project-stores")
OUT  = os.path.join(BASE, "crosscheck")
os.makedirs(OUT, exist_ok=True)

# --- the 6 unk rows: (TIS id, EXP id) ---
UNK = [
    ("TIS05396-001-001", "EXP-00872"),
    ("TIS05744-001-013", "EXP-00876"),
    ("TIS06381-001",     "EXP-01130"),
    ("TIS06381-002",     "EXP-01130"),
]

def tis_variants(tis):
    """All formatting variants for a TIS id, weakest->strongest is encoded via label."""
    digits = re.match(r"TIS(\d+)", tis).group(1)        # e.g. 05396
    prefix = "TIS" + digits                              # TIS05396
    return [
        ("dashed",      tis),                            # TIS05396-001-001
        ("underscored", tis.replace("-", "_")),          # TIS05396_001_001
        ("nosep",       tis.replace("-", "")),           # TIS05396001001
        ("bare_prefix", prefix),                         # TIS05396
        ("bare_digits", digits),                         # 05396  (catch-all superset)
    ]

def exp_variants(exp):
    digits = re.search(r"(\d+)", exp).group(1)           # 00872
    return [
        ("dashed",      exp),                            # EXP-00872
        ("nosep",       "EXP" + digits),                 # EXP00872
        ("underscored", "EXP_" + digits),                # EXP_00872
        ("bare_digits", digits),                         # 00872  (catch-all superset)
    ]

# specificity ranking for "strongest_variant_matched"
STRENGTH = {"dashed": 4, "underscored": 3, "nosep": 2, "bare_prefix": 1, "bare_digits": 0}

def platform_guess(name):
    n = name.lower()
    if "xetg" in n or "xenium" in n: return "Xenium"
    if "visium" in n:                return "Visium"
    return "?"

# --- load store index so we record ALL 33 stores (incl. 6 failed/empty) ---
idx = pd.read_csv(os.path.join(PS, "_store_index.csv"))
all_stores = list(idx["store_name"])
status_by_store = dict(zip(idx["store_name"], idx["status"]))

# --- preload per-store name lists (chunked for the big PLMyeloma file) ---
store_names = {}   # store -> list of (file_id, name, type)
store_count = {}
for store in all_stores:
    path = os.path.join(PS, f"{store}_contents.csv")
    if not os.path.exists(path):
        store_names[store] = []
        store_count[store] = 0
        continue
    rows = []
    for chunk in pd.read_csv(path, usecols=lambda c: c in ("id", "name", "type"),
                             chunksize=100000, dtype=str):
        chunk = chunk.fillna("")
        rows.extend(zip(chunk.get("id", ""), chunk["name"], chunk.get("type", "")))
    store_names[store] = rows
    store_count[store] = len(rows)

# --- run the search ---
long_rows, hit_rows = [], []
for tis, exp in UNK:
    searches = [("TIS", lbl, pat) for lbl, pat in tis_variants(tis)] + \
               [("EXP", lbl, pat) for lbl, pat in exp_variants(exp)]
    paired = {"TIS": exp, "EXP": tis}
    for id_type, lbl, pat in searches:
        patl = pat.lower()
        for store in all_stores:
            files = store_names[store]
            matches = [(fid, name, ftype) for (fid, name, ftype) in files if patl in name.lower()]
            long_rows.append({
                "unk_id": tis, "paired_id": paired[id_type], "id_type": id_type,
                "variant_label": lbl, "search_pattern": pat, "store_name": store,
                "store_status": status_by_store.get(store, "?"),
                "store_file_count": store_count[store],
                "match_count": len(matches), "matched": "Y" if matches else "N",
            })
            for fid, name, ftype in matches:
                hit_rows.append({
                    "unk_id": tis, "id_type": id_type, "variant_label": lbl,
                    "search_pattern": pat, "store_name": store,
                    "file_id": fid, "file_name": name, "file_type": ftype,
                    "platform_guess": platform_guess(name),
                })

long_df = pd.DataFrame(long_rows)
hits_df = pd.DataFrame(hit_rows)

# --- summary: one row per unk_id ---
summ = []
for tis, exp in UNK:
    h = hits_df[hits_df["unk_id"] == tis] if not hits_df.empty else hits_df
    # dedupe hits by (store, file) for counting
    if not h.empty:
        hd = h.drop_duplicates(["store_name", "file_id", "file_name"])
        stores_hit = sorted(hd["store_name"].unique())
        other = sorted(s for s in stores_hit if s != "PCD_Spatial")
        strongest = max(hd["variant_label"], key=lambda l: STRENGTH[l])
        plats = sorted(hd["platform_guess"].unique())
        total = len(hd)
    else:
        stores_hit, other, strongest, plats, total = [], [], "", [], 0
    summ.append({
        "unk_id": tis, "paired_exp": exp, "total_hits": total,
        "stores_with_hits": ";".join(stores_hit),
        "stores_other_than_PCD_Spatial_with_hits": ";".join(other),
        "strongest_variant_matched": strongest,
        "platform_guesses": ";".join(plats),
    })
summ_df = pd.DataFrame(summ)

long_df.to_csv(os.path.join(OUT, "crosscheck_long.csv"), index=False)
hits_df.to_csv(os.path.join(OUT, "crosscheck_hits.csv"), index=False)
summ_df.to_csv(os.path.join(OUT, "crosscheck_summary.csv"), index=False)

# --- console report ---
print(f"crosscheck_long.csv  rows: {len(long_df)}  (expected {len(UNK)}*9*{len(all_stores)} = {len(UNK)*9*len(all_stores)})")
print(f"crosscheck_hits.csv  rows: {len(hits_df)}")
print(f"crosscheck_summary.csv rows: {len(summ_df)}")
print("\n=== SUMMARY ===")
print(summ_df.to_string(index=False))
if not hits_df.empty:
    other_hits = hits_df[hits_df["store_name"] != "PCD_Spatial"].drop_duplicates(["store_name","file_id","file_name"])
    print(f"\nHits in stores OTHER than PCD_Spatial: {len(other_hits)}")
    if not other_hits.empty:
        print(other_hits[["unk_id","store_name","variant_label","platform_guess","file_name"]].to_string(index=False))
