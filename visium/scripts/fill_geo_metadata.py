"""Populate the 9 GEO sample-attribute columns on the `visiuim-geo` tab.

Fills columns 1-9 (*Sample name .. description) for every GEO-bound sample,
sourced from the `main` tab. Leaves the filename columns (10-21) untouched.
Sample set = all 30 `main` rows. (TIS06381-001 / -002, EXP-01130, will be
processed through Space Ranger later; their metadata is included now, file
columns filled once output exists.) Descriptions are kept public-clean:
internal QC/tracking tokens are stripped.
"""
import openpyxl

XLSX = "visium/data/visium data tracking.xlsx"
DROP = set()  # no samples excluded; all 30 main-tab rows are GEO-bound

# Notes tokens that are internal-only and must not leak into the public GEO description.
DROP_PATTERNS = ["fastq flag", "pair with", "visium not confirmed",
                 "visium status unclear", "st-000", "same patient", "shared w"]


def tx_map(v):
    if not v:
        return "unknown"
    v = str(v).strip()
    return {"Pre-treatment": "pre-treatment", "Post-treatment": "post-treatment",
            "Unknown": "unknown", "N/A": "unknown"}.get(v, v.lower())


def clean_notes(notes):
    if not notes:
        return []
    toks = []
    for chunk in str(notes).split("|"):
        for tok in chunk.split(";"):
            tok = tok.strip()
            if not tok:
                continue
            low = tok.lower()
            if any(p in low for p in DROP_PATTERNS):
                continue
            toks.append(tok)
    return toks


def build_title(sid, exp, dis):
    d = dis if dis else "Bone Marrow"
    return f"{d} Visium — {sid} ({exp})"


def build_desc(dis, site, src, notes):
    parts = []
    if dis:
        parts.append(str(dis))
    if site:
        parts.append(f"Site: {site}")
    if src:
        parts.append(f"Source: {src}")
    parts += clean_notes(notes)
    return "; ".join(parts)


wb = openpyxl.load_workbook(XLSX)
main = wb["main"]
ws = wb["visiuim-geo"]

r = 2  # row 1 is the header; real samples start at row 2
n = 0
for row in main.iter_rows(min_row=2, values_only=True):
    sid = row[0]
    if not sid or sid in DROP:
        continue
    exp, tissue, dis, site, tx, src, notes = row[1], row[12], row[13], row[14], row[15], row[16], row[17]
    tissue = tissue or "Bone Marrow"
    vals = [
        sid,                                 # *Sample name
        build_title(sid, exp, dis),          # *title
        "Homo sapiens",                      # *organism
        tissue,                              # *tissue
        "FFPE",                              # *tissue preservation method
        tx_map(tx),                          # treatment
        "N/A",                               # genotype
        "total RNA",                         # *molecule
        build_desc(dis, site, src, notes),   # description
    ]
    for c, v in enumerate(vals, 1):
        ws.cell(row=r, column=c, value=v)
    r += 1
    n += 1

wb.save(XLSX)
print(f"Wrote GEO metadata for {n} samples (rows 2-{r-1}), columns 1-9. Filename columns untouched.")
