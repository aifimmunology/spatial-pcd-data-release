# PCD Visium Spatial Data — Inventory, Validation & GEO Prep

> **Status:** In progress (working repository). This README is a living document for the
> AIFI data team. Sections marked _TODO_ or `<!-- FILL -->` need owner input.

## What this project is

This repo organizes and validates **10x Genomics Visium** spatial transcriptomics data for the
plasma-cell-disorder / multiple-myeloma (PCD/MM) study, sourced from **HISE** project stores at
AIFI, and prepares it for a **public NCBI GEO submission**.

It drives three intertwined goals:

1. **Internal data inventory / tracking** — maintain an accurate catalog of which Visium samples
   exist, which HISE store holds them, and what each archive actually contains.
2. **Validation / QC** — confirm each sample has a real processed Space Ranger output and determine
   whether antibody-capture (ADT) data is genuinely present.
3. **Assemble the GEO tracking sheet** — gather all the metadata and file locations for every
   processed run into one complete tracking sheet (the deliverable). The Visium side of this
   project is **inventory + metadata assembly only — we do not reformat, split, or move/upload the
   data files themselves.**

The deliverable — and the single source of truth for sample status — is
**`visium/data/visium data tracking.xlsx`** (sheets: `main`, `supplemental`, `visiuim-geo`, plus
`ADT markers`). Everything else in this repo exists to populate, verify, or extend that sheet.

## Scope at a glance

- **30 in-study Visium samples** tracked across multiple experiments/cohorts (`EXP-008xx`–`EXP-011xx`,
  batches `SB002`–`SB008`). (TIS05744 was removed — not part of this study.)
- **All 30 are GEO-bound** and populated on the `visiuim-geo` sheet (sample metadata). **28** have
  processed Space Ranger output; the remaining **2** (`TIS06381-001` / `-002`, EXP-01130, AIFI
  Internal) are to be processed later — file columns pending. See
  [Outstanding needs](#outstanding-needs--todo).
- Primary HISE store: **`PCD_Spatial`**. Cross-checked against **33** project stores total.

## Repository layout

This modality is self-contained under `visium/` at the repo root (see the top-level `README.md` for
how Visium and Xenium fit together):

```
visium/
├── README.md                          # this file
├── notebooks/
│   └── visium_data.ipynb              # scratch / HISE exploration
├── scripts/                           # the workflows (see below)
│   └── archive/                       # superseded reformat/stage scripts (out of scope)
├── reference/
│   ├── final_visium_structure_example.txt   # reference only — a GEO layout example, NOT a target we conform to
│   └── hisepy.txt                     # hisepy API notes
└── data/                              # NOT committed (.gitignore) — large / generated
    ├── visium data tracking.xlsx      # ★ the deliverable / source of truth (main / supplemental / visiuim-geo / ADT markers)
    ├── project-stores/                # *_contents.csv listings + _store_index.csv
    ├── visium-zip-trees/              # one .md content tree per zip
    ├── adt-panels/                    # per-run feature_reference.csv + _panel_comparison.csv
    ├── visium-data/                   # staged downloads (transient)
    └── crosscheck/                    # adt_presence.csv, crosscheck_*.csv
```

> **Note:** `visium/data/` is git-ignored (via `**/data/`; project-store contents, zip trees,
> tracking sheet are large and/or regenerable). Only code and reference docs are committed.

## Workflows & scripts

All scripts assume the working directory `/home/workspace/PCD_data` (the repo root) and reference
their inputs/outputs via the modality-prefixed path `visium/data/...`. Where noted, they also need a
`hisepy` session with access to the relevant HISE project stores.

### 1. Discover & pull store metadata
- **`pull_stores.py`** — for 27 target stores, calls `hp.list_files_in_project_store(store)` and
  writes `<store>_contents.csv`; 6 stores already on disk are reused. Writes a paper-trail
  `_store_index.csv` covering all 33. Read-only, fault-tolerant per store.
- **`list_pcd_matches.py`** — list files in `PCD_Spatial` matching a sample's TIS/EXP id (or
  formatting variants).

### 2. Targeted file download
- **`visium_data_download.py`** — download a fixed set of files from `PCD_Spatial` into
  `visium/data/visium-data/`. HISE stages files under a read-only FUSE path
  (`~/input/<id>/<store>/<uuid>/<file>`); the script locates and moves them, preserving the
  `AIFI-…/` prefix as a subfolder to disambiguate shared basenames.

### 3. Lightweight zip inventory (no large zips kept)
- **`visium_zip_tree.py`** — the core inventory tool. For each zip in `PCD_Spatial`:
  download → read internal structure → write a `.md` content tree → delete → next. Peak disk = one
  zip. (hisepy offers no byte-range read, so each zip is fully fetched once.)
- **`visium_zip_tree_sb009.py`**, **`visium_zip_tree_exp00899.py`** — batch-specific variants.

### 4. Cross-check unresolved samples
- **`run_crosscheck.py`** — LOCAL ONLY (no HISE). Cross-checks the unresolved ("unk") rows against
  all 33 `*_contents.csv` using TIS/EXP id formatting variants (dashed / underscored / nosep /
  bare-prefix / bare-digits). Writes `crosscheck_hits.csv`, `crosscheck_long.csv`,
  `crosscheck_summary.csv`. (This drove locating the now-resolved samples; the remaining 3 sit in an
  account this session can't read — see Outstanding needs.)

### 5. ADT (antibody-capture) detection
- **`check_adt_presence.py`** — cross-checks every sample in `main` against its zip tree and
  classifies ADT status. **Key decision:** the word "antibody" alone is unreliable — Space Ranger
  always emits an `analysis/` skeleton for antibody_capture. Real ADT is confirmed only by
  `outs/feature_reference.csv` **plus** a non-trivial `antibody_capture_graphclust/clusters.csv`.
  Verdicts: `ADT_PRESENT` / `SCAFFOLD_ONLY` / `GEX_ONLY` / `NO_TREE`. Output: `adt_presence.csv`.

### 6. Single-artifact extraction
- **`download_cytassist_image.py`** — pull only `outs/spatial/cytassist_image.tiff` (the canonical
  copy) from a zip via the central directory, without unpacking the whole archive.
- **`download_web_summary.py`** — pull only `web_summary.html` from a zip the same way.

### 7. ADT panel assembly (metadata)
- **`extract_feature_references.py`** — pull only `outs/feature_reference.csv` (the antibody panel
  definition) from every processed run's zip into `visium/data/adt-panels/`, sample-prefixed. Reads the
  member from the zip's central directory via the FUSE mount; never unpacks the full archive.
- **`compare_adt_panels.py`** — LOCAL ONLY. Normalize each `feature_reference.csv` to its
  antibody-marker set and group runs by identical panel. Writes `_panel_comparison.csv` and
  `_panel_groups.txt`.
- **`build_adt_markers_tab.py`** — write the `ADT markers` tab into the tracking sheet (one row per
  marker per detected panel group, preserving `feature_reference.csv` order).
- **`add_panel_group_column.py`** — add the `ADT panel group` column to the `main` tab, matching
  each sample's run to its panel group from `_panel_comparison.csv`. Idempotent.

### 8. Tracking sheet & GEO metadata
- **`fill_filenames.py`** — fill the cellranger/zip output filename column from a curated
  `(TIS, EXP) → zip name` mapping derived from the zip trees.
- **`fill_geo_tab.py`** — build the `visiuim-geo` sheet from `adt_presence.csv` (every distinct
  `(id, experiment)` with a processed tree), recording the GEO metadata and per-run file
  references. (No files are produced or moved — this populates the tracking sheet only.)

### Archived — out of scope
`scripts/archive/` holds the earlier reformat-and-stage tooling (`stage_run_outs.py`,
`assemble_geo_run.py`, `build_deliverable_tree_md.py`, `compare_template_vs_output.py`). These
restructured runs into the deliverable layout and previewed/staged the result — **no longer in
scope** now that the Visium deliverable is the tracking sheet, not reformatted data files. Kept for
reference.

## Data provenance & key decisions

- **HISE FUSE behavior** — `download_from_project_store` does **not** copy onto writable disk; it
  exposes files via a read-only FUSE cache mount under `~/input/`. Locate downloads by glob; they
  need not (and cannot) be deleted from there. (Scripts that "delete" zips operate on copies staged
  into `visium/data/visium-data/`.)
- **HISE endpoint throttling** — `download_from_project_store` re-calls
  `list_files_in_project_store` internally on *every* call; listing the large `PCD_Spatial` store in
  a tight loop can trigger HTTP 500s (`/hydration/analysis/projectStore/files`). Space out calls /
  add backoff; single calls succeed.
- **Combined matrix** — Space Ranger emits ONE combined matrix per run holding both Gene Expression
  and Antibody Capture features (they are distinguishable by `feature_type`). Noted as a property of
  the data; we are **not** splitting or re-deriving matrices in this project.
- **ID formatting variants** — TIS/EXP ids appear inconsistently across stores (dashes vs
  underscores vs none vs bare digits). Matching must try all variants; bare-digits is the
  catch-all superset. See `run_crosscheck.py` / `list_pcd_matches.py`.
- **ADT detection** — see script 5 above; do not trust the "antibody" keyword.
- **Platform ambiguity** — some ids that look like Visium samples resolve to **Xenium** or other
  assays in other stores; cross-check flags `platform_guesses` accordingly.
- **HISE environment / cross-account access** — the `main` sheet's `HISE environment` column records
  which account/bucket holds each sample. Samples in the **AIFI Internal / AIFI Internal Account**
  environment live in a separate bucket that **this HISE account cannot read**, so their metadata and
  file locations can't be captured from this session (see Outstanding needs).
- **`final_visium_structure_example.txt` is reference only** — it documents an example per-sample
  GEO directory structure (gene_expression raw/processed/spatial, ADT_conventional, ADT_CST_bridge,
  web_summary, H&E). It is **not** a layout we conform files to; the deliverable here is the tracking
  sheet, not a restructured data package.

## Outstanding needs / TODO

- [x] ~~**Access the 3 samples in the AIFI Internal account**~~ — **closed; internal-account access
      no longer needed.** These three lived in the unreadable **AIFI Internal / AIFI Internal Account**
      bucket, but the blocker is now moot:
      - `TIS06381-001 / EXP-01130` (AIFI Internal) — **processed; Space Ranger 3.0.0.** Output zip
        `SB009-TIS06381-001-EXP-01130.zip` lives in the AIFI Internal account's **`visium_data`**
        project store (readable only from an AIFI-Internal HISE session, not the PCD session that
        holds `PCD_Spatial`). Version read from the run's `_cmdline`
        (`/home/jupyter/spaceranger-3.0.0/spaceranger count`). Metadata is in the GEO tab.
      - `TIS06381-002 / EXP-01130` (AIFI Internal) — **processed; Space Ranger 3.0.0.** Output zip
        `SB009-TIS06381-002-EXP-01130.zip`, same `visium_data` store, same version and `_cmdline`
        source. Metadata is in the GEO tab.

        *(Note: these two ran under Space Ranger 3.0.0, not the 2.1.1 used across `PCD_Spatial` —
        they were processed separately and later. `spaceranger version` values pending write into
        the tracking sheet, which is only reachable from the PCD session.)*
      - `TIS05396-001-001 / EXP-00872` — **reprocessed independently** through Space Ranger (see
        "Locate FASTQs" below), so the AIFI-Internal processed output was never needed. Remaining
        step is uploading the rerun output — see new TODO below.

      *(The former TIS05744 "unk" samples — `EXP-00876` / `EXP-00899` / `EXP-00900` — were likewise
      removed: not part of this study.)*
- [x] ~~**Upload the reprocessed EXP-00872 output to the PCD_Spatial project store.**~~ — **Done
      (2026-07-01).** Ingested via the watchfolder as `EXP00872_TIS05396_A1.zip` (HISE UUID
      `e94a7fae-67c3-49ed-917e-854414cb07c1`) in `PCD_Spatial`; confirmed present in the live store.
- [ ] **Locate FASTQs + inputs to rerun Space Ranger** — find the raw sequencing FASTQs (and any
      other data-specific inputs) needed to rerun `spaceranger count` from scratch. This is
      **CytAssist Gene & Protein Expression** (antibody capture).

      **Update (2026-06-30): FASTQs are not the universal blocker we first thought.** The PCD
      watchfolder has ingested per-experiment GEX+Antibody FASTQs into the **readable PCD_Spatial
      store**, under `SB00x/FC*_done/outs/fastq_path/<flowcell>/EXP-XXXXX_{G,P}_S#_L00#_{I1,I2,R1,R2}_001.fastq.gz`.
      Confirmed present for SB001/SB002: EXP-00826/827/848/849/851/852/875/876 and **EXP-00872**.
      (Where FASTQs are not surfaced, the prior fallback still holds: 23/27 runs carry a
      `possorted_genome_bam.bam` → regenerate via 10x `bamtofastq`; the 4 BAM-less runs cannot.)

      **Worked example — EXP-00872 / TIS05396-001-001 (the lone FASTQ-only sample): reprocessed ✅.**
      It is area **A1** of slide **V53U12-040**; its slide-mate **EXP-00875** (area B1) has a complete
      processed run, whose `_cmdline`/`_invocation` (spaceranger **2.1.1**) we recovered to lock in
      identical inputs (saved to `visium/data/crosscheck/exp00875-run-params/`). All inputs sourced:
      - **FASTQs** (`--libraries`, GEX `EXP-00872_G` + Antibody `EXP-00872_P`) — ✅ in PCD_Spatial; staged to `visium/data/rerun-exp00872/fastqs/`.
      - **CytAssist image** (`--cytaimage`) — ✅ `..._A_tis05396-001-001.tif` (uploaded CytAssist zip).
      - **Slide layout** (`--slide V53U12-040 --area A1` + `--loupe-alignment V53U12-040-A1.json`) — ✅ uploaded CytAssist zip.
      - **Transcriptome** (`--transcriptome`) — ✅ `refdata-gex-GRCh38-2020-A` on disk (exact name match).
      - **Probe set** (`--probe-set`) — ✅ `Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv`, byte-identical to the bundled copy.
      - **Feature reference** (`--feature-ref`) — ✅ stock 10x `Visium_Human_Immune_Cell_Profiling_Panel_v1.0.csv` (bundled). Antibody panel **assumed identical** to the slide-mate (same run/batch, stock panel); not independently confirmed for EXP-00872.

      **Rerun completed (2026-06-30)** under **Space Ranger 2.1.1** (matching the original; 4.0.1 is
      also installed but would not be bit-identical). Output:
      `visium/data/rerun-exp00872/output/EXP00872_TIS05396_A1/outs/` (launcher: `visium/data/rerun-exp00872/run_spaceranger_2.1.1.sh`).
      QC: 3,857 spots under tissue, 17,978 genes; ADT real (30M antibody reads, median 3,222 UMIs/spot,
      8 antibody clusters, 3% unrecognized → panel assumption validated). Data flag: *Estimated UMIs from
      genomic DNA = 0.49* (high gDNA background — a property of this Normal-BM sample, not a run error).

      **➡️ Uploaded (2026-07-01) ✅** — the regenerated run is in the **`PCD_Spatial`** project store
      as `EXP00872_TIS05396_A1.zip` (HISE UUID `e94a7fae-67c3-49ed-917e-854414cb07c1`), ingested via
      the watchfolder. Nothing further outstanding for this sample.
- [x] ~~**Verify the ADT panels against the tracking sheet**~~ — **Done.** `extract_feature_references.py`
      + `compare_adt_panels.py` assembled and grouped the per-run `feature_reference.csv` panels, and
      each sample's panel is recorded in the sheet (the `ADT panel group` column in `main` + the
      `ADT markers` tab). Knowing which sample has which panel was the goal; that mapping is tracked.
- [ ] **Record H&E image locations** — capture the EVOS microscopy image location from the `H&E url`
      column (and note `cytassist_image.tiff` presence) in the tracking sheet. We record locations,
      not the images themselves — no images are pulled or moved.
- [ ] **Add a "HISE output UUID" column to the Visium tracking sheet.** Each sample's Space Ranger
      output zip corresponds to one UUID in the HISE database; record that UUID per sample so each
      row links back to its concrete deliverable in HISE.
- [ ] **Populate GEO sample metadata on the `visiuim-geo` tab** — fill the sample-attribute
      columns for all samples: *Sample name, *title, *organism, *tissue, *tissue preservation
      method, treatment, genotype, *molecule, and description. **Not** the filename columns.
- [x] ~~**Compile all FASTQ location metadata + HISE UUIDs.**~~ — **Done (2026-07-14).**
      FASTQ provenance for all 30 samples was recovered from each run's Space Ranger
      `_invocation` (sample-name prefix + source flowcell/dir per library;
      `visium/data/crosscheck/fastq_invocations.csv`). **Only SB002's raw FASTQs
      (EXP-00872, EXP-00875; flowcell HNKWTDMXY) are archived in HISE** — in `PCD_Spatial`,
      with per-file UUIDs (32 files; `fastq_sb002_perfile.csv`). Exhaustive search of all
      **42** HISE project stores + the file ledger confirmed the SB003/SB004/SB005/SB008
      raw FASTQs were **never ingested to HISE** (only the processed output zips are); they
      exist only on the compute-node paths recorded in the invocations. Recorded in the
      tracking sheet: `main.fastq filename` (prefix + source, all rows), `main.fastq UUID`
      (tab pointer for SB002, "not in HISE" otherwise), and a new **`FASTQ files` tab**
      (one row per run×library, with the 32 SB002 UUIDs). Scripts:
      `extract_fastq_invocations.py`, `fill_fastq_metadata.py`. Locations/UUIDs only — no
      FASTQs moved.
- [ ] **Resolve two FASTQ data-integrity flags surfaced during FASTQ compilation** (also
      noted in the sheet's `Notes` / `FASTQ files` tab):
      - **EXP-01092 (SB008) FASTQs deliberately switched** — the two runs' `_invocation`
        sample_ids are `TIS06377-001_modified_switched_fastq` / `TIS06392-002_modified_switched_fastq`;
        GEX and Antibody were cross-sourced between the two tissues (hence the
        "corrected_sr_outs" zips). Confirm this was intentional and correct.
      - **TMA00147 output = TMA00144 data** — `SB008_EXP-01106_TMA00147_full.zip` carries
        TMA00144's `sample_id` + capture area (V53S25-332-A1) and uses the `SB008TMA00144`
        FASTQs, i.e. it appears to duplicate/mislabel TMA00144. Verify TMA00147 has a
        genuine distinct run.

## Setup & usage

- **Python deps:** `hisepy`, `pandas`, `openpyxl` (plus stdlib `zipfile`, `glob`).
- **HISE access:** an authenticated `hisepy` session with read access to the PCD project stores is
  required for any script that touches HISE (pull / download / zip-tree). The cross-check and ADT
  scripts are **local-only** and need no HISE.
- **Typical order:** `pull_stores.py` → `visium_zip_tree.py` → `check_adt_presence.py` →
  `run_crosscheck.py` → `extract_feature_references.py` → `compare_adt_panels.py` →
  `build_adt_markers_tab.py` / `add_panel_group_column.py` → `fill_filenames.py` / `fill_geo_tab.py`.

<!-- FILL: study/PI context, GEO series accession (once assigned), submission deadline, data owner -->
```
Study / cohort context:   <!-- FILL -->
GEO series accession:     <!-- FILL (pending) -->
Submission target date:   <!-- FILL -->
Data owner / contact:     <!-- FILL -->
```
