#!/usr/bin/env python3
"""Reproducible pull of HISE project-store metadata listings (READ-ONLY).

For each of 27 target stores, call hp.list_files_in_project_store(store) and
write the resulting DataFrame to <out_dir>/<store>_contents.csv. Failures are
caught per-store and recorded; the run continues. A paper-trail index covering
all 33 stores (27 pulled + 6 reused) is written to _store_index.csv.
"""
import os
import datetime
import traceback

import pandas as pd
import hisepy as hp

OUT_DIR = "/home/workspace/PCD_data/visium/data/project-stores"
EXPECTED_COLS = ["id", "name", "availability", "type", "metadata", "store_name"]

# 27 stores to pull fresh
TO_PULL = [
    "cohorts", "fh1_bmp", "rds", "fh1_data", "CD123_B", "br1_br2",
    "up2_prelim", "PedvsSenior", "fh1_treat", "olink_cohort", "RA-2021-09",
    "CompBio_demo", "UCSDCU_preRA", "Mono_Sub_A", "brictyof", "teaseqpert",
    "Alt_splice", "preRA_Bcell", "myeloid_Tseq", "UCSDCU_Y4", "imm_age",
    "ifns-immune", "StopRA", "scatac_label", "U01T1D", "Dyna_IHandA",
    "U01T1Dfastqs",
]

# 6 stores already on disk (do NOT re-pull)
REUSED = [
    "PCD_Spatial", "PCDSpatiFHMM", "PLMyeloma", "Renji_MM",
    "MM_Bone_Marr", "FH_MM_Year4",
]


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    index_rows = []

    # ---- TASK 1: pull the 27 stores ----
    for store in TO_PULL:
        ts = utc_now()
        try:
            df = hp.list_files_in_project_store(store)
            if df is None:
                df = pd.DataFrame(columns=EXPECTED_COLS)
            out_path = os.path.join(OUT_DIR, f"{store}_contents.csv")
            df.to_csv(out_path, index=False)
            n = len(df)
            print(f"[ok]   {store}: {n} rows -> {out_path}")
            index_rows.append({
                "store_name": store, "row_count": n, "source": "pulled",
                "status": "ok", "error": "", "pulled_at_utc": ts,
            })
        except Exception as e:  # noqa: BLE001 - catch-all per requirement
            err = f"{type(e).__name__}: {e}".strip().replace("\n", " ")
            print(f"[FAIL] {store}: {err}")
            traceback.print_exc()
            index_rows.append({
                "store_name": store, "row_count": "", "source": "pulled",
                "status": "failed", "error": err, "pulled_at_utc": ts,
            })

    # ---- TASK 2 (reused part): read row counts from existing CSVs ----
    for store in REUSED:
        path = os.path.join(OUT_DIR, f"{store}_contents.csv")
        try:
            n = len(pd.read_csv(path))
            status, err = "ok", ""
        except Exception as e:  # noqa: BLE001
            n = ""
            status = "failed"
            err = f"{type(e).__name__}: {e}".strip().replace("\n", " ")
        index_rows.append({
            "store_name": store, "row_count": n, "source": "reused",
            "status": status, "error": err, "pulled_at_utc": "preexisting",
        })

    # ---- write the index ----
    index_df = pd.DataFrame(index_rows, columns=[
        "store_name", "row_count", "source", "status", "error", "pulled_at_utc",
    ])
    index_path = os.path.join(OUT_DIR, "_store_index.csv")
    index_df.to_csv(index_path, index=False)

    # ---- summary ----
    pulled = index_df[index_df["source"] == "pulled"]
    n_ok = int((pulled["status"] == "ok").sum())
    n_fail = int((pulled["status"] == "failed").sum())
    print("\n===== SUMMARY =====")
    print(f"Pulled stores: {n_ok} ok / {n_fail} failed (of {len(pulled)})")
    if n_fail:
        print("Failures:")
        for _, r in pulled[pulled["status"] == "failed"].iterrows():
            print(f"  - {r['store_name']}: {r['error']}")
    print(f"Index rows: {len(index_df)} (-> {index_path})")
    on_disk = [f for f in os.listdir(OUT_DIR) if f.endswith("_contents.csv")]
    print(f"_contents.csv files on disk: {len(on_disk)}")
    print("\nNewly pulled row_counts:")
    for _, r in pulled.iterrows():
        print(f"  {r['store_name']}: {r['row_count']} ({r['status']})")


if __name__ == "__main__":
    main()
