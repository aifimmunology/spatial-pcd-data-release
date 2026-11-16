#!/usr/bin/env python3
"""Match Google-bucket Xenium outputs to tracking-sheet samples not yet in HISE.

Run from the repo root:  python xenium/scripts/match_gbucket_to_hise.py

For every `main`-tab row with a blank `file name` (i.e. not ingested into HISE),
find the matching whole-output folder(s) in gs://temp-xenium-hise-transfer/ using:
    slide id (7-digit zero-padded)  == filename slide token
    AND  section id == filename section token          (TMA-named runs)
         OR cohort  == / is a prefix of section token  (FH-subject-named runs)
"""

import re
import subprocess
import pandas as pd

BUCKET = "gs://temp-xenium-hise-transfer/"
TRACKING = "xenium/data/xenium data tracking.xlsx"

OUTPUT_RE = re.compile(
    r"output-(XETG\d+)__(\d+)__(.+)__(\d{8})__(\d{6})/?$"
)


def list_bucket_outputs():
    """Return list of (folder_name, slide_token, section_token) for top-level outputs."""
    out = subprocess.run(
        ["gsutil", "ls", BUCKET], capture_output=True, text=True, check=True
    ).stdout.splitlines()
    parsed = []
    for line in out:
        name = line.rstrip("/").split("/")[-1]
        m = OUTPUT_RE.match(name)
        if m:
            _, slide, section, _, _ = m.groups()
            parsed.append((name, slide, section))
    return parsed


def matches(row, slide_tok, section_tok):
    """True if a bucket output (slide_tok, section_tok) belongs to this main-tab row."""
    slide = row["slide id"]
    if pd.isna(slide) or str(int(slide)).zfill(7) != slide_tok:
        return False
    section_id = str(row["section id"])
    cohort = str(row["cohort"])
    return (
        section_tok == section_id
        or section_tok == cohort
        or section_tok.startswith(cohort + "_")
    )


def main():
    outputs = list_bucket_outputs()
    df = pd.read_excel(TRACKING, sheet_name="main")
    missing = df[df["file name"].isna()]

    print(f"{len(outputs)} output folders in bucket; "
          f"{len(missing)} un-ingested rows to match.\n")

    matched_folders = set()
    for _, row in missing.iterrows():
        hits = [name for (name, s, sec) in outputs if matches(row, s, sec)]
        matched_folders.update(hits)
        tag = f"{row['block id']} / {row['section id']} / slide {int(row['slide id'])} / {row['cohort']}"
        if hits:
            for h in hits:
                print(f"[MATCH] {tag}\n        -> {h}")
        else:
            print(f"[NONE ] {tag}  -- no bucket output found")

    print(f"\n{len(matched_folders)} distinct bucket folder(s) matched:")
    for f in sorted(matched_folders):
        print(f"  {f}")


if __name__ == "__main__":
    main()
