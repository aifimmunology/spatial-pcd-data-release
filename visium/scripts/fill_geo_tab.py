import openpyxl, csv

XLSX = "visium/data/visium data tracking.xlsx"

# Build the sample list from adt_presence.csv: every distinct (id, experiment)
# with a processed Space Ranger tree on disk. Distinct full IDs are kept as
# separate GEO samples; only literal (id, experiment) repeats are collapsed.
EXP = {}  # id -> experiment
seen = set()
INCLUDE = []
for r in csv.DictReader(open("visium/data/crosscheck/adt_presence.csv")):
    if r["tree_found"] != "Y":
        continue
    key = (r["id"], r["experiment"])
    if key in seen:
        continue
    seen.add(key)
    INCLUDE.append(r["id"])
    EXP[r["id"]] = r["experiment"]

# Standardized processed deliverables (per final_visium_structure_example.txt)
PROC_SUFFIXES = [
    "raw_feature_bc_matrix.h5",
    "filtered_feature_bc_matrix.h5",
    "molecule_info.h5",
    "tissue_hires_image.png",
    "tissue_lowres_image.png",
    "tissue_positions.csv",
    "scalefactors_json.json",
    "aligned_fiducials.jpg",
    "ADT_raw_feature_bc_matrix.h5",
    "ADT_filtered_feature_bc_matrix.h5",
]
RAW_SUFFIXES = ["R1.fastq.gz", "R2.fastq.gz"]  # placeholder convention - confirm

# ---- pull metadata from main tab ----
wb_d = openpyxl.load_workbook(XLSX, data_only=True)
main = wb_d["main"]
meta = {}  # id -> dict
for row in main.iter_rows(min_row=2, values_only=True):
    sid = row[0]
    if not sid or sid in meta:
        continue
    meta[sid] = dict(tissue=row[10], disease=row[11], location=row[12],
                     tx=row[13], source=row[14], notes=row[15])

def tx_map(v):
    if not v:
        return ""
    v = str(v).strip()
    return {"Pre-treatment": "pre-treatment", "Post-treatment": "post-treatment",
            "Unknown": "unknown", "N/A": "N/A"}.get(v, v.lower())

def build_title(sid, m):
    exp = EXP.get(sid)
    suffix = f"{sid} ({exp})" if exp else sid
    if sid.startswith("TMA"):
        return f"Bone Marrow TMA Visium - {suffix}"
    d = m.get("disease")
    if d:
        return f"Bone Marrow Visium - {d} - {suffix}"
    return f"Bone Marrow Visium - {suffix}"

def build_desc(m):
    parts = []
    if m.get("disease"): parts.append(str(m["disease"]))
    if m.get("location"): parts.append(f"Site: {m['location']}")
    if m.get("source"): parts.append(f"Source: {m['source']}")
    if m.get("notes"): parts.append(str(m["notes"]))
    return "; ".join(parts)

# ---- write GEO tab ----
wb = openpyxl.load_workbook(XLSX)
ws = wb["visiuim-geo"]

META_HEADERS = ["*Sample name", "*title", "*organism", "*tissue",
                "*tissue preservation method", "treatment", "genotype",
                "*molecule", "description"]
header = list(META_HEADERS)
header += ["*processed data file"] * len(PROC_SUFFIXES)
header += ["*raw data file"] * len(RAW_SUFFIXES)

# clear existing used range (rows 1..max) across enough columns
ncol = len(header)
for r in range(1, ws.max_row + 1):
    for c in range(1, max(ws.max_column, ncol) + 1):
        ws.cell(row=r, column=c).value = None

# header
for c, h in enumerate(header, 1):
    ws.cell(row=1, column=c, value=h)

# keep the original example row as row 2
example = ["Sample 1", "Bone Marrow Visium", "Human", "Bone Marrow", "FFPE",
           "pre-treatment", "N/A", "total RNA", "Example row - format reference"]
for c, v in enumerate(example, 1):
    ws.cell(row=2, column=c, value=v)

# real samples from row 3
r = 3
for sid in INCLUDE:
    m = meta.get(sid, {})
    tissue = m.get("tissue") or ("Bone Marrow" if sid.startswith("TIS") else "")
    vals = [
        sid,                      # *Sample name
        build_title(sid, m),      # *title
        "Human",                  # *organism
        tissue,                   # *tissue
        "FFPE",                   # *tissue preservation method
        tx_map(m.get("tx")),      # treatment
        "N/A",                    # genotype
        "total RNA",              # *molecule
        build_desc(m),            # description
    ]
    vals += [f"{sid}_{suf}" for suf in PROC_SUFFIXES]
    vals += [f"{sid}_{suf}" for suf in RAW_SUFFIXES]
    for c, v in enumerate(vals, 1):
        ws.cell(row=r, column=c, value=v)
    r += 1

wb.save(XLSX)
print(f"Wrote {len(INCLUDE)} samples (rows 3-{r-1}); {ncol} columns; header row + example row 2.")
