import os
import shutil
import tempfile
import zipfile

import pandas as pd
import scanpy as sc
from spatialdata_io import xenium

DATA_DIR    = '/home/workspace/PCD_data/xenium/fh1_pretreatment/data'
OUTS_DIR    = os.path.join(DATA_DIR, 'outs')
SAMPLES_DIR = os.path.join(DATA_DIR, 'samples')
OBJECTS_DIR = os.path.join(DATA_DIR, 'archive', 'lognorm_azu1_kept', 'objects')

os.makedirs(OUTS_DIR, exist_ok=True)
os.makedirs(SAMPLES_DIR, exist_ok=True)
os.makedirs(OBJECTS_DIR, exist_ok=True)

# 14 single-patient MM Pre-Treatment outputs flagged `downstream analysis (Fig3)`
# in the tracking sheet. Each entry: sample_id -> (hise_uuid, cohort_subject, timepoint).
# The first 10 already have a per-sample h5ad copied from the earlier fh1 run
# (loaded from disk below); the last 4 are downloaded fresh from HISE.
SAMPLES = {
    'TIS04399-001': ('e841cbfc-851f-4781-a707-5cf06e3e6068', 'FH1001', 'MM Pre-Treatment'),
    'TIS02434-001': ('f2957fe1-aa58-4d75-adfe-1cabd4f9bbe0', 'FH1009', 'MM Pre-Treatment'),
    'TIS02437-001': ('cba78e0f-62d1-4592-9581-a02a5cea0fae', 'FH1011', 'MM Pre-Treatment'),
    'TIS04398-001': ('b8e8fe26-5d7c-4aee-9850-2651a7073ceb', 'FH1015', 'MM Pre-Treatment'),
    'TIS02611-001': ('8460b141-a598-47c1-a79e-d30658698302', 'FH1016', 'MM Pre-Treatment'),
    'TIS02603-001': ('d14fff16-f2e1-4dcd-bf11-fcbfe364dada', 'FH1017', 'MM Pre-Treatment'),
    'TIS02604-001': ('fee36c64-5a50-4b7f-868e-616088d9701f', 'FH1018', 'MM Pre-Treatment'),
    'TIS02613-001': ('69ea81f2-1751-43ad-baf0-1952f593f4e6', 'FH1021', 'MM Pre-Treatment'),
    'TIS02421-001': ('91211de2-f801-4f2b-af3a-18017dc44005', 'FH1003', 'MM Pre-Treatment'),
    'TIS02430-001': ('c3e68735-26be-471c-834e-24833fbde0dd', 'FH1007', 'MM Pre-Treatment'),
    'TIS02424-001': ('40a0ef40-6b5a-43db-845e-d5bc58faa685', 'FH1004', 'MM Pre-Treatment'),
    'TIS02439-001': ('807357ef-5733-427f-b6ce-b828754f180c', 'FH1012', 'MM Pre-Treatment'),
    'TIS02594-001': ('f3c5af40-e58c-4143-86ce-0194ca9eab26', 'FH1006', 'MM Pre-Treatment'),
    # TIS04397-001 (FH1010) dropped: it is on a 377-gene panel, not the 480-gene
    # panel of the other 13. Including it collapsed the join='inner' gene set to
    # 258; without it the 13 share the full 480.
}


def _read_table_direct(dest):
    adata = sc.read_10x_h5(os.path.join(dest, 'cell_feature_matrix.h5'))
    adata.var_names_make_unique()
    cells = pd.read_parquet(os.path.join(dest, 'cells.parquet'))
    cells['cell_id'] = cells['cell_id'].astype(str)
    cells = cells.set_index('cell_id').reindex(adata.obs_names)
    adata.obsm['spatial'] = cells[['x_centroid', 'y_centroid']].to_numpy()
    adata.obs = adata.obs.join(cells.drop(columns=['x_centroid', 'y_centroid']))
    return adata


def download_and_read(name, uuid):
    dest = os.path.join(OUTS_DIR, name)
    if not os.path.isdir(dest):
        # Imported here, not at module load, so a run that only reuses cached
        # per-sample h5ads never triggers hisepy's account/auth check.
        import hisepy as hp
        print(f'{name}: downloading from HISE (UUID {uuid})…')
        (zip_path,) = hp.cache_files(file_ids=[uuid])
        with tempfile.TemporaryDirectory() as tmp:
            with zipfile.ZipFile(zip_path) as zf:
                zf.extractall(tmp)
            inner = next(
                os.path.join(tmp, d) for d in os.listdir(tmp)
                if os.path.isdir(os.path.join(tmp, d))
            )
            shutil.move(inner, dest)
        print(f'{name}: extracted → {dest}')

    # We only use the cell table (counts + centroid coords in obsm['spatial']),
    # so skip the images, boundaries, labels, and transcripts. This also avoids a
    # corrupt morphology_focus TIFF that crashes the image reader on some samples.
    try:
        sdata = xenium(
            dest,
            cells_boundaries=False, nucleus_boundaries=False,
            cells_labels=False, nucleus_labels=False,
            transcripts=False, morphology_mip=False, morphology_focus=False,
            aligned_images=False,
        )
        adata = sdata.tables['table']
    except Exception as err:
        # Some bundles ship a 0-byte experiment.xenium (and other stripped files);
        # spatialdata_io parses it for version detection and dies before reading
        # any data. The count matrix and cell centroids are intact, so read them
        # directly instead.
        print(f'{name}: spatialdata_io reader failed ({type(err).__name__}); '
              f'falling back to direct cell_feature_matrix + cells.parquet read')
        adata = _read_table_direct(dest)

    adata.obs_names_make_unique()
    return adata, dest


adatas = {}
for name, (uuid, subject, timepoint) in SAMPLES.items():
    out_path = os.path.join(SAMPLES_DIR, f'{name}.h5ad')
    if os.path.exists(out_path):
        print(f'{name}: {out_path} exists — loading from disk')
        adatas[name] = sc.read_h5ad(out_path)
        continue

    adata, extracted = download_and_read(name, uuid)
    adata.obs['sample'] = name
    adata.obs['cohort_subject'] = subject
    adata.obs['treatment_timepoint'] = timepoint
    adata.write_h5ad(out_path)
    print(f'{name}: {adata.n_obs:,} cells, {adata.n_vars} genes → {out_path}')
    adatas[name] = adata

    # Extracted Xenium outputs are 2-8 GB each; delete after reading so peak
    # disk holds only one extracted output at a time.
    shutil.rmtree(extracted)

combined = sc.concat(list(adatas.values()), join='inner')
combined.obs_names_make_unique()
print(f'\nCombined: {combined.n_obs:,} cells x {combined.n_vars} genes')
print('\nCells per timepoint:')
print(combined.obs['treatment_timepoint'].value_counts().to_string())
print('\nCells per patient:')
print(combined.obs['cohort_subject'].value_counts().to_string())

combined_path = os.path.join(OBJECTS_DIR, 'combined_raw.h5ad')
combined.write_h5ad(combined_path)
print(f'\nWrote {combined_path}')
