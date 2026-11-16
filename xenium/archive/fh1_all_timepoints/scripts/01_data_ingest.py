import os
import shutil
import tempfile
import zipfile

import pandas as pd
import scanpy as sc
import hisepy as hp
from spatialdata_io import xenium

DATA_DIR    = '/home/workspace/PCD_data/xenium/archive/fh1_all_timepoints/data'
OBJECTS_DIR = os.path.join(DATA_DIR, 'objects')
OUTS_DIR    = os.path.join(DATA_DIR, 'outs')
SAMPLES_DIR = os.path.join(DATA_DIR, 'samples')

os.makedirs(OUTS_DIR, exist_ok=True)
os.makedirs(SAMPLES_DIR, exist_ok=True)
os.makedirs(OBJECTS_DIR, exist_ok=True)

# 44 single-patient FH1 outputs (multi-patient TMA outputs are excluded).
# Each entry: sample_id -> (hise_uuid, cohort_subject, treatment_timepoint).
# FH1002 at Post Transplant 1 year has two sections (TIS02087-001, TIS06674-001);
# both are kept.
FH1_SAMPLES = {
    'TIS02422-001': ('3a57b697-a54b-4c11-a6c0-dd5abdd0499b', 'FH1003', 'MM End Induction 1st Draw'),
    'TIS02428-001': ('af7a7a30-7933-48b3-bfd9-a7dd5d7f16ca', 'FH1005', 'MM End Induction 1st Draw'),
    'TIS02595-001': ('1cf7b648-1750-42b9-ae1b-260e0eb1a6da', 'FH1006', 'MM End Induction 1st Draw'),
    'TIS02433-001': ('b07364e6-37be-4bf5-8d21-be1053b0de67', 'FH1008', 'MM End Induction 1st Draw'),
    'TIS02438-001': ('572b5020-4c89-4ad6-8c9d-dff525814702', 'FH1011', 'MM End Induction 1st Draw'),
    'TIS02597-001': ('df21c8fb-23b7-45e7-87e1-bd5347450318', 'FH1014', 'MM End Induction 1st Draw'),
    'TIS02605-001': ('b4114a94-fa64-401a-801e-38ae1429457a', 'FH1017', 'MM End Induction 1st Draw'),
    'TIS04963-001': ('606b86b7-b269-45f6-954d-05e4c1ecd333', 'FH1022', 'MM End Induction 1st Draw'),
    'TIS06485-001': ('810d07c2-c746-4383-9ec2-662b47fae7f3', 'FH1028', 'MM End Induction 1st Draw'),
    'TIS02087-001': ('6a6977af-5771-4d3c-8d9d-934885cab56f', 'FH1002', 'MM Post Transplant 1 year'),
    'TIS06674-001': ('5c4d1cc4-3d35-4331-9a11-8d190fa15b63', 'FH1002', 'MM Post Transplant 1 year'),
    'TIS04748-001': ('ae541f6a-2262-4e50-acc7-af8ccc158954', 'FH1003', 'MM Post Transplant 1 year'),
    'TIS04750-001': ('1ee867b6-326c-4dec-ac5b-867edc0f4f04', 'FH1006', 'MM Post Transplant 1 year'),
    'TIS04751-001': ('b544b8cb-3185-4683-8999-d8b8abfc1c38', 'FH1007', 'MM Post Transplant 1 year'),
    'TIS04755-001': ('8a338d33-b392-45e4-98eb-ad3429cbb2e8', 'FH1008', 'MM Post Transplant 1 year'),
    'TIS04434-001': ('dcb529b6-803d-40c5-a0e1-5d5a0fd1334b', 'FH1009', 'MM Post Transplant 1 year'),
    'TIS04969-001': ('e58081d9-a024-450d-9278-bc1bcd79f875', 'FH1017', 'MM Post Transplant 1 year'),
    'TIS05005-001': ('386e5819-da69-44f7-9b3f-e9d12c6a5836', 'FH1021', 'MM Post Transplant 1 year'),
    'TIS04967-001': ('ecb5134f-8031-4012-bc31-1817bdc3d56b', 'FH1005', 'MM Post Transplant 2 year'),
    'TIS04757-001': ('2e723938-18e7-4e6e-ae8a-2cd831d1e7b0', 'FH1007', 'MM Post Transplant 2 year'),
    'TIS04994-001': ('35e6a706-8036-4e74-89af-e1114b344793', 'FH1009', 'MM Post Transplant 2 year'),
    'TIS02423-001': ('06cc928b-01b3-448a-9b8c-07712810c218', 'FH1003', 'MM Post Transplant 90 Days'),
    'TIS02426-001': ('34650d6c-8aff-49c3-9c79-abc4b1251094', 'FH1004', 'MM Post Transplant 90 Days'),
    'TIS02596-001': ('b33213f7-9b2f-4b4e-a0fe-f3b6ead2f4b4', 'FH1006', 'MM Post Transplant 90 Days'),
    'TIS02431-001': ('8d888814-7a1a-43f0-bc35-487c4d4b4d5a', 'FH1007', 'MM Post Transplant 90 Days'),
    'TIS02601-001': ('520d2517-caa7-41e2-bf22-3f482355b364', 'FH1014', 'MM Post Transplant 90 Days'),
    'TIS02617-001': ('57ee56c2-9834-43a8-bd62-342c56fc2044', 'FH1018', 'MM Post Transplant 90 Days'),
    'TIS04996-001': ('02031a9c-fbbf-49cc-a513-db8771146dc7', 'FH1022', 'MM Post Transplant 90 Days'),
    'TIS05004-001': ('a6e7b0cf-0978-4c6b-b970-f5c1c913b78b', 'FH1024', 'MM Post Transplant 90 Days'),
    'TIS06484-001': ('65e15824-5a99-41f4-bc6a-742f74cd25a1', 'FH1027', 'MM Post Transplant 90 Days'),
    'TIS04399-001': ('e841cbfc-851f-4781-a707-5cf06e3e6068', 'FH1001', 'MM Pre-Treatment'),
    'TIS02421-001': ('91211de2-f801-4f2b-af3a-18017dc44005', 'FH1003', 'MM Pre-Treatment'),
    'TIS02430-001': ('c3e68735-26be-471c-834e-24833fbde0dd', 'FH1007', 'MM Pre-Treatment'),
    'TIS02434-001': ('f2957fe1-aa58-4d75-adfe-1cabd4f9bbe0', 'FH1009', 'MM Pre-Treatment'),
    'TIS02437-001': ('cba78e0f-62d1-4592-9581-a02a5cea0fae', 'FH1011', 'MM Pre-Treatment'),
    'TIS04398-001': ('b8e8fe26-5d7c-4aee-9850-2651a7073ceb', 'FH1015', 'MM Pre-Treatment'),
    'TIS02611-001': ('8460b141-a598-47c1-a79e-d30658698302', 'FH1016', 'MM Pre-Treatment'),
    'TIS02603-001': ('d14fff16-f2e1-4dcd-bf11-fcbfe364dada', 'FH1017', 'MM Pre-Treatment'),
    'TIS02604-001': ('fee36c64-5a50-4b7f-868e-616088d9701f', 'FH1018', 'MM Pre-Treatment'),
    'TIS02613-001': ('69ea81f2-1751-43ad-baf0-1952f593f4e6', 'FH1021', 'MM Pre-Treatment'),
    'TIS04431-001': ('fe50bfe9-ea10-4811-bd6d-bf2937a6eefb', 'FH1022', 'MM Pre-Treatment'),
    'TIS04707-001': ('60162f9e-2d50-405a-bcbc-89aa38c06b2b', 'FH1023', 'MM Pre-Treatment'),
    'TIS04962-001': ('b9dd1998-442b-4d0b-906c-24a17b4cd782', 'FH1024', 'MM Pre-Treatment'),
    'TIS05002-001': ('fc84eb04-1c5c-4881-b8ff-bf45f8a7e3e9', 'FH1028', 'MM Pre-Treatment'),
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
    # corrupt morphology_focus TIFF that crashes the image reader on some samples
    # (e.g. TIS02597-001).
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
        # any data (e.g. TIS04757-001). The count matrix and cell centroids are
        # intact, so read them directly instead.
        print(f'{name}: spatialdata_io reader failed ({type(err).__name__}); '
              f'falling back to direct cell_feature_matrix + cells.parquet read')
        adata = _read_table_direct(dest)

    adata.obs_names_make_unique()
    return adata, dest


adatas = {}
for name, (uuid, subject, timepoint) in FH1_SAMPLES.items():
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
    # disk holds only one extracted output at a time, not all 44.
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
