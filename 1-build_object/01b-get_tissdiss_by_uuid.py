import os
import pandas as pd
import hisepy as hp

metadata = pd.read_csv('../0-metadata_fixer/updated_metadata.csv')
wanted_files = [os.path.basename(p) for p in metadata[metadata['flag'].isna()]['file_path']]

store = hp.list_files_in_project_store('PCD_Spatial')
store['basename'] = store['name'].fillna('').map(os.path.basename)
file_ids = store[store['basename'].isin(wanted_files)]['id'].tolist()

downloaded_paths = hp.cache_files(file_ids=file_ids)
print(f"{len(downloaded_paths)} of {len(wanted_files)} files downloaded")
