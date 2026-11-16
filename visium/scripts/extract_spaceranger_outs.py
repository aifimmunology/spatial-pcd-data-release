"""
Extract only the outs/ directory from each Visium Space Ranger zip, skipping the BAM.

Reads the tracking sheet and project-store CSV to locate each zip, then extracts
into visium/data/spaceranger_outputs/<sample_id>/outs/.

Run from anywhere — paths are resolved relative to this script's location.

Usage:
    python3 extract_spaceranger_outs.py
"""

import os
import zipfile
import pandas as pd

SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT    = os.path.dirname(SCRIPT_DIR)          # visium/
TRACKING_SHEET = os.path.join(REPO_ROOT, 'data', 'visium data tracking.xlsx')
STORE_CSV      = os.path.join(REPO_ROOT, 'data', 'project-stores', 'PCD_Spatial_contents.csv')
STORE_BASE     = '/home/workspace/input/815275457/PCD_Spatial'
OUT_BASE       = os.path.join(REPO_ROOT, 'data', 'spaceranger_outputs')

# Files inside outs/ to skip (large, not needed for scanpy/squidpy analysis)
SKIP_SUFFIXES = ('possorted_genome_bam.bam',)

# ── Build sample → zip path mapping ─────────────────────────────────────────

xl  = pd.read_excel(TRACKING_SHEET, sheet_name='main')
csv = pd.read_csv(STORE_CSV)

uuid_to_storename = csv.set_index('id')['name'].to_dict()

samples = xl[
    xl['Space Ranger output'].str.startswith('Yes', na=False) &
    xl['spaceranger UUID'].notna()
][['id', 'spaceranger UUID']].copy()

jobs = []
missing_zip = []

for _, row in samples.iterrows():
    sample_id = row['id']
    uuid      = row['spaceranger UUID']
    store_name = uuid_to_storename.get(uuid)
    if not store_name:
        missing_zip.append(sample_id)
        continue
    zip_filename = os.path.basename(store_name)
    zip_path = os.path.join(STORE_BASE, uuid, store_name)
    if not os.path.exists(zip_path):
        missing_zip.append(f"{sample_id} (zip not found: {zip_path})")
        continue
    jobs.append({'sample_id': sample_id, 'zip_path': zip_path})

print(f"Samples to extract : {len(jobs)}")
if missing_zip:
    print(f"Skipped (no zip)   : {missing_zip}")

# ── Extract outs/ from each zip ──────────────────────────────────────────────

def extract_outs(zip_path, dest_dir):
    """
    Extract all members under .../outs/ (except the BAM) into dest_dir.
    Strips everything before 'outs/' in the internal path, so the result is
    dest_dir/outs/<file>.
    """
    skipped_files = []
    extracted = 0

    with zipfile.ZipFile(zip_path, 'r') as zf:
        members = zf.infolist()

        for member in members:
            name = member.filename

            # Only process paths that contain /outs/ or are the outs/ entry itself
            if '/outs/' not in name and not name.endswith('/outs/'):
                continue

            # Derive the relative path starting from outs/
            if '/outs/' in name:
                rel = 'outs/' + name.split('/outs/', 1)[1]
            else:
                rel = 'outs/'

            # Skip the BAM (too large, not needed for scanpy/squidpy)
            if any(rel.endswith(s) for s in SKIP_SUFFIXES):
                skipped_files.append(rel)
                continue

            target = os.path.join(dest_dir, rel)

            if member.is_dir():
                os.makedirs(target, exist_ok=True)
                continue

            os.makedirs(os.path.dirname(target), exist_ok=True)

            # Skip if already extracted
            if os.path.exists(target) and os.path.getsize(target) > 0:
                continue

            with zf.open(member) as src, open(target, 'wb') as dst:
                dst.write(src.read())

            extracted += 1

    return extracted, skipped_files


succeeded = []
failed    = []

for i, job in enumerate(jobs, 1):
    sample_id = job['sample_id']
    zip_path  = job['zip_path']
    dest_dir  = os.path.join(OUT_BASE, sample_id)

    print(f"\n[{i}/{len(jobs)}] {sample_id}")
    print(f"  zip  : {os.path.basename(zip_path)}")
    print(f"  dest : {dest_dir}")

    os.makedirs(dest_dir, exist_ok=True)

    try:
        n, skipped = extract_outs(zip_path, dest_dir)
        print(f"  extracted : {n} files")
        if skipped:
            print(f"  skipped   : {skipped}")
        succeeded.append(sample_id)
    except Exception as e:
        print(f"  ERROR : {e}")
        failed.append({'sample_id': sample_id, 'error': str(e)})

# ── Summary ──────────────────────────────────────────────────────────────────

print("\n" + "=" * 60)
print(f"Done.  Succeeded: {len(succeeded)}  Failed: {len(failed)}")
if failed:
    print("\nFailed:")
    for f in failed:
        print(f"  {f['sample_id']} — {f['error']}")

print(f"\nOutputs in: {OUT_BASE}/")
