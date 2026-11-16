import csv

path = "/home/workspace/PCD_data/visium/data/visium_data_tracking .csv"

# (id, experiment) -> zip file name, derived from visium-zip-trees/*.md
mapping = {
    ("TIS05397-001-001", "EXP-00875"): "EXP00875_Visium_BM.zip",
    ("TIS06392-002", "EXP-01092"): "SB008_EXP-01092_TIS06392-002_corrected_sr_outs.zip",
    ("TIS06377-001", "EXP-01092"): "SB008_EXP-01092_TIS06377-001_corrected_sr_outs.zip",
    ("TIS05681-001", "EXP-01095"): "SB008_EXP-01095_TIS05681-001_sr_outs.zip",
    ("TIS05683-001", "EXP-01095"): "SB008_EXP-01095_TIS05683-001_sr_outs.zip",
    ("TIS05684-001-003", "EXP-01106"): "SB008_EXP-01106_TIS05684-001-003.zip",
    ("TIS06402-001-004", "EXP-01106"): "SB008_EXP-01106_TIS06402-001-004.zip",
    ("TMA00144", "EXP-01106"): "SB008_EXP-01106_TMA00144.zip",
    ("TMA00147", "EXP-01106"): "SB008_EXP-01106_TMA00147_full.zip",
    ("TIS04788-001-004", "EXP-00947"): "TIS04788_EXP00947_SB5.zip",
    ("TIS04789-001-013", "EXP-00948"): "TIS04789_EXP00948_SB5.zip",
    ("TIS05671-001-001", "EXP-00949"): "TIS05671_EXP00949_SB5.zip",
    ("TIS05672-001-006", "EXP-00950"): "TIS05672_EXP00950_SB5.zip",
    ("TIS05392-001", "EXP-00963"): "TIS05392_EXP00963_AACR_SB5.zip",
    ("TIS05396-001", "EXP-00963"): "TIS05396_EXP00963_AACR_SB5.zip",
    ("TIS05397-001", "EXP-00963"): "TIS05397_EXP00963_AACR_SB5.zip",
    ("TIS04788-001", "EXP-00963"): "TIS04788_EXP00963_AACR_SB5.zip",
    ("TIS05393-001-004", "EXP-00934"): "EXP00934_BM.zip",
    ("TIS05395-001-004", "EXP-00935"): "EXP00935_BM.zip",
    ("TIS05394-001-003", "EXP-00936"): "EXP00936_BM.zip",
    ("TIS05398-001-001", "EXP-00937"): "EXP00937_BM.zip",
    ("TIS05392-001-001", "EXP-00897"): "EXP00897_BM.zip",
    ("TIS05391-001-001", "EXP-00898"): "EXP00898_BM.zip",
}

with open(path, newline="") as f:
    rows = list(csv.reader(f))

header = rows[0]
fn_idx = header.index("file name")
id_idx = header.index("id")
exp_idx = header.index("experiment")

filled = 0
for row in rows[1:]:
    if len(row) <= fn_idx:
        continue
    key = (row[id_idx], row[exp_idx])
    if not row[fn_idx].strip() and key in mapping:
        row[fn_idx] = mapping[key]
        filled += 1

with open(path, "w", newline="") as f:
    csv.writer(f).writerows(rows)

print(f"Filled {filled} cells")
