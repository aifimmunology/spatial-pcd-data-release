"""
Download all Visium Space Ranger output zips from the PCD_Spatial HISE project store.

Run from the repository root:
    python3 visium/scripts/download_spaceranger_outputs.py
"""

import os
import sys
import pandas as pd
import unittest.mock as mock

# Patch broken import before loading hisepy
with mock.patch.dict('sys.modules', {
    'hisepy.training_job': mock.MagicMock(),
    'hisepy.ai_ray_transformer': mock.MagicMock()
}):
    import hisepy as hp

STORE = 'PCD_Spatial'

SCRIPT_DIR     = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT      = os.path.dirname(SCRIPT_DIR)          # visium/scripts/ -> visium/ -> ..
TRACKING_SHEET = os.path.join(REPO_ROOT, 'data', 'visium data tracking.xlsx')
STORE_CSV      = os.path.join(REPO_ROOT, 'data', 'project-stores', 'PCD_Spatial_contents.csv')

# ── Load tracking sheet ──────────────────────────────────────────────────────

xl = pd.read_excel(TRACKING_SHEET, sheet_name='main')
samples = xl[
    xl['Space Ranger output'].str.startswith('Yes', na=False) &
    xl['spaceranger UUID'].notna()
][['id', 'spaceranger UUID', 'project store']].copy()

print(f"Samples to download: {len(samples)}")

# ── Match UUIDs to file names in the project store CSV ──────────────────────

csv = pd.read_csv(STORE_CSV)
uuid_to_name = csv.set_index('id')['name'].to_dict()

downloads = []
missing   = []

for _, row in samples.iterrows():
    uuid = row['spaceranger UUID']
    store = row['project store']
    name = uuid_to_name.get(uuid)
    if name:
        downloads.append({'sample': row['id'], 'store': store, 'name': name, 'uuid': uuid})
    else:
        missing.append({'sample': row['id'], 'uuid': uuid})

print(f"  Matched in CSV : {len(downloads)}")
if missing:
    print(f"  NOT in CSV     : {len(missing)}")
    for m in missing:
        print(f"    {m['sample']}  uuid={m['uuid']}")

# ── Download ─────────────────────────────────────────────────────────────────

succeeded = []
failed    = []

for i, d in enumerate(downloads, 1):
    print(f"\n[{i}/{len(downloads)}] {d['sample']}")
    print(f"  store : {d['store']}")
    print(f"  file  : {d['name']}")
    try:
        result = hp.download_from_project_store(
            store_name=d['store'],
            file_name=d['name']
        )
        print(f"  result: {result}")
        succeeded.append(d)
    except Exception as e:
        print(f"  ERROR : {e}")
        failed.append({**d, 'error': str(e)})

# ── Summary ──────────────────────────────────────────────────────────────────

print("\n" + "=" * 60)
print(f"Done.  Succeeded: {len(succeeded)}  Failed: {len(failed)}")
if failed:
    print("\nFailed downloads:")
    for f in failed:
        print(f"  {f['sample']}  —  {f['error']}")
