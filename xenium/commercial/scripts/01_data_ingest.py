import logging
import os
import shutil
import tempfile
import zipfile

import scanpy as sc
import hisepy as hp
from spatialdata_io import xenium

DATA_DIR    = '/home/workspace/PCD_data/xenium/commercial/data'
OBJECTS_DIR = os.path.join(DATA_DIR, 'objects')
OUTS_DIR    = os.path.join(DATA_DIR, 'outs')
SAMPLES_DIR = os.path.join(DATA_DIR, 'samples')
LOGS_DIR    = os.path.join(DATA_DIR, 'logs')

for directory in (OBJECTS_DIR, OUTS_DIR, SAMPLES_DIR, LOGS_DIR):
    os.makedirs(directory, exist_ok=True)

logging.basicConfig(
    filename=os.path.join(LOGS_DIR, '01_data_ingest.log'),
    level=logging.INFO,
    format='%(asctime)s %(message)s',
)

XENIUM_SAMPLES = {
    'TIS05393':     '70f7c0d6-8d7e-4b86-acf3-ca80137fd849',  # Multiple Myeloma
    'TIS04788':     '38fe47a0-c09f-4373-9a8a-8a957c40aa9f',  # Multiple Myeloma
    'TIS05392':     '0c879d26-970f-47fc-a705-88b103b4de80',  # Multiple Myeloma
    'TIS05394':     'd7b4b14e-9856-4169-a363-213b4e29a99e',  # Multiple Myeloma
    'TIS04789':     '7c3113c3-23bf-4840-8f45-61505048322c',  # Multiple Myeloma
    'TIS05397_008': '16752455-bd23-4f88-8d05-2108e7778496',  # Normal BM (2024)
}


def download_and_extract(uuid, dest):
    (zip_path,) = hp.cache_files(file_ids=[uuid])
    with tempfile.TemporaryDirectory() as tmp:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(tmp)
        extracted_folder = next(
            os.path.join(tmp, d) for d in os.listdir(tmp)
            if os.path.isdir(os.path.join(tmp, d))
        )
        shutil.move(extracted_folder, dest)


adatas = {}
for name, uuid in XENIUM_SAMPLES.items():
    sample_dir = os.path.join(OUTS_DIR, name)
    if not os.path.isdir(sample_dir):
        download_and_extract(uuid, sample_dir)
        logging.info('%s: downloaded from HISE (UUID %s) to %s', name, uuid, sample_dir)

    # Creates the AnnData from the raw Xenium output: .X = raw integer counts (cells x 377 genes),
    # obs = per-cell metadata from the vendor, obsm['spatial'] = each cell's x/y position in microns.
    adata = xenium(sample_dir).tables['table']
    # Rewrites any duplicated cell IDs within the sample by appending a suffix.
    adata.obs_names_make_unique()
    sample_path = os.path.join(SAMPLES_DIR, f'{name}.h5ad')
    # Saves the unmodified per-sample object to disk.
    adata.write_h5ad(sample_path)
    logging.info('%s: %d cells, %d genes -> %s', name, adata.n_obs, adata.n_vars, sample_path)
    adatas[name] = adata

# Stacks the six samples into one object, keeping only genes present in all of them,
# and adds an obs column 'sample' recording each cell's source sample.
combined = sc.concat(adatas, label='sample', join='inner')
# Cell IDs repeat across samples, so this suffixes the duplicates to make them unique.
combined.obs_names_make_unique()
combined_path = os.path.join(OBJECTS_DIR, 'combined_raw.h5ad')
combined.write_h5ad(combined_path)
logging.info(
    'combined: %d cells, %d genes -> %s\n%s',
    combined.n_obs, combined.n_vars, combined_path,
    combined.obs['sample'].value_counts().to_string(),
)
