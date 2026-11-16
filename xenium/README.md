# PCD Xenium Spatial Data

> **Status:** Scaffolding — not yet started. This README is a stub to be filled in. Sections marked
> _TODO_ or `<!-- FILL -->` need owner input.

## What this project is

The Xenium side of the PCD/MM spatial project: organize and validate **10x Genomics Xenium** data
sourced from **HISE** project stores at AIFI, toward the same **NCBI GEO submission** as the Visium
data. See the top-level [`README.md`](../README.md) for how the modalities fit together, and
[`visium/README.md`](../visium/README.md) for the established pattern this side will mirror.

<!-- FILL: goals for the Xenium side — inventory / QC / GEO metadata assembly, and any differences
     from the Visium workflow. -->

## Scope at a glance

<!-- FILL: number of Xenium samples, experiments/cohorts, HISE store(s), how many are GEO-bound. -->

## Repository layout

Each analysis is named after the samples it uses, not the paper figure it feeds, and is
self-contained: its scripts, data, figures and notebooks all sit under one folder.

```
xenium/
├── README.md
├── commercial/            # 6 commercial-biobank sections (paper Figure 2 D–F) — see commercial/README.md
│   ├── README.md
│   ├── scripts/           # numbered pipeline 01 … 08 (+ archive/)
│   ├── data/              # NOT committed — objects/, samples/, resolution_sweep/, qc_diagnostics/, logs/, outs/
│   ├── figures/           # NOT committed — qc_calibration/, processing/, umap_by_sample/, umap_celltypes/,
│   │                      #   spatial_celltypes/, qc_diagnostics/, archive/
│   └── notebooks/archive/ # the original exploratory notebooks for these 6 samples
├── fh1_pretreatment/      # 13 FH1 MM Pre-Treatment samples (paper Figure 3) — see fh1_pretreatment/README.md
│   ├── README.md
│   ├── scripts/           # config.py + _pipeline.py + 01_ingest / 02_filter / 03_process (+ archive/)
│   ├── data/              # NOT committed — samples/, outs/, archive/ (earlier runs)
│   ├── figures/           # NOT committed — archive/ (earlier runs)
│   └── notebooks/
├── archive/
│   └── fh1_all_timepoints/   # retired: all 44 FH1 samples, every timepoint — scripts/, data/, figures/
├── scripts/               # shared GEO-metadata / tracking-sheet helpers (see below)
├── notebooks/             # shared exploration (explore_multitissue_outputs.ipynb)
├── reference/             # reference docs
└── data/                  # NOT committed — shared files used by no single analysis
    ├── xenium-zip-trees/      # file tree per analyzer version
    ├── xenium data tracking.xlsx
    ├── mmc2.xlsx
    ├── Content_*.xlsx
    └── README.md
```

Every `data/` and `figures/` folder, at any depth, is git-ignored. Inside each analysis, `archive/`
holds earlier runs and retired scripts: kept so past results stay traceable, never written to again.

## Workflows & scripts

GEO-metadata workflow (run from the repo root):
- `scripts/xenium_zip_tree.py` — survey one output bundle per analyzer version (grouped by run
  timestamp); writes a file tree per version to `data/xenium-zip-trees/`. Needs the HISE-authenticated
  environment (`/home/workspace/environment/scanpyv3/bin/python`).
- `scripts/fill_geo_tab.py` — fill the `geo` tab's 10 processed-file columns (sample-id-prefixed) and
  the raw placeholders; also fills `treatment` from the `main` tab. No HISE needed.
- `scripts/add_main_columns.py` — add `xenium analyzer software version` and the
  `downstream analysis (Fig2 D-F)` flag to the `main` tab. Needs HISE.

## Data provenance & key decisions

### GEO submission — deposit file convention

The GEO submission has **91 samples** (sample-level, the `geo`-tab rows), drawn from the 97
section-level rows on the `main` tab. Six samples were run twice; we keep the **newer** output
(latest run timestamp): TIS02439-001, TIS04436-001, TIS05395-001, TIS05396-001, TIS05397-001,
TIS05398-001.

**Analyzer-version survey.** Xenium output bundles are zips in the `PCD_Spatial` HISE project store,
one per `hise uuid`. Their internal file names are generic (the sample identity is only in the outer
zip name), so every bundle from the same Xenium Analyzer run shares an identical file list. Runs are
keyed by the `YYYYMMDD__HHMMSS` timestamp in the zip name — 20 unique timestamps across the 97 rows,
spanning **five analyzer versions**: 1.7.1.0, 2.0.0.10, 3.1.0.4, 3.2.0.6, 3.2.0.7.
`scripts/xenium_zip_tree.py` downloaded one bundle per timestamp and wrote a file tree per version to
`data/xenium-zip-trees/`. The file set is uniform across versions except: the morphology image layout
(a flat `morphology_focus.ome.tif` in 1.7.1.0 vs. a `morphology_focus/` directory from 2.0.0.10 on),
and `transcripts.csv.gz` (present only in 1.7.1.0 and 2.0.0.10).

**Curated processed files (per sample).** Ten files, all present in every analyzer version, listed in
the `geo` tab as `<sample>_<file>`:

    cell_feature_matrix.h5, cells.parquet, cell_boundaries.parquet, nucleus_boundaries.parquet,
    transcripts.parquet, gene_panel.json, experiment.xenium, metrics_summary.csv,
    analysis_summary.html, morphology.ome.tif

We use `morphology.ome.tif` (the full pyramidal tissue image) because it is uniform across all five
versions, sidestepping the version-dependent `morphology_focus` naming. Excluded: all `.zarr.zip`
archives (redundant; the interactive 10x Xenium Explorer is out of scope), the redundant alternate
encodings (`cell_feature_matrix/` MEX trio, and the `.csv.gz` twins of the parquet tables), and the
vendor `analysis/` clustering subtree.

**Raw files: deferred.** The two `*raw data file` columns hold a `TBD` placeholder. Xenium has no
FASTQs; the natural "raw" is the morphology image plus `transcripts.parquet`, but those read as
processed too — the sequencing-shaped GEO template fits Xenium poorly, so the raw design is a separate
decision (see TODO).

**Reference and size.** A near-identical study is on GEO as **GSE345302** ("...tumor microenvironment
during multiple myeloma evolution", human MM bone-marrow Xenium) — its per-sample deposit uses the same
canonical file set. The morphology images dominate size (~3 GB each), so the total is likely **> 1 TB**;
GEO requires contacting them before submitting above that threshold.

### Editing the tracking sheet with openpyxl — column placement

When appending a column to a tab, anchor on the **last column that has a header value**
(`max(c.column for c in ws[1] if c.value is not None) + 1`), **not** `ws.max_column`. The workbook
carries empty-but-styled cells well past the real table, so `ws.max_column` overshoots (e.g. the
`main` tab once reported 38 when the table ended at 18) and drops new columns into a gap.
`scripts/populate_hise_uuids.py` and `scripts/add_main_columns.py` use the correct anchor.

## Outstanding needs / TODO

- [x] **Define Xenium scope** — done (2026-07-16). Deliverable is to find all Xenium samples,
      organize their locations and metadata, split the multi-tissue outputs, and produce QC
      metrics for everything (not reformatting or re-ingesting data).
- [ ] **Fix swapped TMA labels (TMA00145 / TMA00142)** — these two Xenium outputs have their TMA
      IDs swapped. Two corrections are needed: (1) swap the TIS ID annotations for TMA00145 and
      TMA00142 in the tracking sheet, and (2) run `xeniumranger rename` on both outputs to correct
      the TMA ID embedded in the filename and output metadata.
- [ ] **Annotate multi-tissue outputs** — 8 combined outputs each bundle 2 FH subjects. We are
      **not** separating them into per-tissue outputs; instead, document which subjects are in each
      output and how a future user can partition them. The two subjects occupy distinct regions of
      the slide, so they can be separated by **manual geometric partitioning — splitting the cells
      by a coordinate boundary** (cross-referenced with the FFPE block layout). Outputs (identified
      as `main`-tab rows sharing one `xenium gbucket output name`):
    - EXP-01104 `TMA00142` (slide 0031088): FH1010 + FH1026
    - EXP-01402 `TMAB00157` (slide 0037044): FH1012 + FH1020
    - EXP-01403 `TMA00159` (slide 0037058): FH1006 + FH2004
    - EXP-01403 `TMA00161` (slide 0037058): FH1004 + FH2003
    - EXP-01506 `TMA00238` (slide 0052114): FH1009 + FH1018
    - EXP-01506 `TMA00240` (slide 0052114): FH1002 + FH1016
    - EXP-01507 `TMA00244` (slide 0052117): FH1017 + FH1018
    - EXP-01507 `TMA00246` (slide 0052117): FH1005 + FH1008
- [ ] **Finish per-tissue-sample QC across all runs** — complete the aggregated Xenium
      `metrics_summary.csv` QC table for every run, mirroring the `qc` tab of
      `data/xenium data tracking.xlsx` (currently 88 regions / 17 runs / 45 columns).
      **Open question for the scientists:** since the 8 combined outputs are no longer being
      separated, decide whether we still want individual per-sample QC metrics for the two bundled
      subjects (which would require partitioning the cells by coordinate boundary first — see
      "Annotate multi-tissue outputs"), or whether we settle for QC at the output level (one row per
      combined output, without separating the bundled subjects).
- [x] **Add HISE UUIDs to the tracking sheet** — done (2026-07-17). Populated a `hise uuid`
      column on the `main` tab for all 95 rows as of that date, by matching each row's
      `xenium gbucket output name` against a live
      `hise.list_files_in_project_store('PCD_Spatial')` listing (name → UUID) via
      `scripts/populate_hise_uuids.py`. All 87 unique output names resolved; the 82 already in the
      `hise`-tab snapshot matched live exactly, and the 5 outputs ingested 2026-07-09 were resolved
      from the live query. (Re-checked 2026-09-18: `main` now has 97 rows, every one carrying a
      UUID, and `hise uuid` is column 17 of 22 — five columns were appended after it since.)
- [x] **Ingest the 9 located outputs into HISE** — done (2026-07-09). Zipped the 5 whole outputs
      from `gs://temp-xenium-hise-transfer/` (store mode, `output-…/` at the archive root) and
      uploaded them to the `pcd_spatial` watchfolder bucket
      `gs://wfpcd_spatialaifid-b34b9c66-aa5d-45d0-9569-f78ca8f6813f/`; all 5 auto-ingested and
      verified available in HISE. Filled the AIFI-prefixed `xenium gbucket output name` for the
      9 `main`-tab rows: slides **52114 / 52117** (EXP-01506 / EXP-01507) as 4
      TMA-named outputs (TMA00238/240/244/246), and slide **37039 / TMA00167** (EXP-01408, FH1017)
      as 1 FH-subject-named output (`FH1017_MM_Post_Transplant_1yr`). The 4 TMA outputs still bundle
      2 FH subjects each — see "Annotate multi-tissue outputs". (The project store itself is not a
      `main`-tab column; it is recorded on the `hise` tab as `store_name`, `PCD_Spatial` throughout.)
- [x] **Diff Google bucket vs. HISE** — located all 9 un-ingested samples as **5 whole unzipped
      outputs** in `gs://temp-xenium-hise-transfer/` (top-level folders), matched by slide id +
      section token via `scripts/match_gbucket_to_hise.py`. Folder names recorded in the
      `xenium gbucket output name` column of the tracking sheet's `main` tab.
- [x] **Assemble the GEO `geo`-tab metadata** — done (2026-09-02). Surveyed the five analyzer versions
      (`scripts/xenium_zip_tree.py`), agreed a 10-file curated processed set, and filled the `geo` tab
      for all 91 samples (`scripts/fill_geo_tab.py`): 10 `*processed data file` columns of
      sample-id-prefixed names, `TBD` in the 2 raw columns, `treatment` from `main`. Added
      `xenium analyzer software version` and a `downstream analysis (Fig2 D-F)` flag to the `main` tab
      (`scripts/add_main_columns.py`). See "Data provenance & key decisions" for the convention.
- [ ] **Design the Xenium raw-data deposit** — the `geo` tab's 2 `*raw data file` columns are `TBD`;
      GEO validation will reject that. Decide what "raw" means for Xenium (candidate: morphology image
      + `transcripts.parquet`, or the whole per-sample bundle) and fill the columns.
- [ ] **Deposit the 6-sample downstream analysis** — as Series-level supplementary files: the processed
      AnnData (`commercial/data/objects/combined_clustered.h5ad`), a per-cell annotation CSV, and the
      per-resolution marker-gene tables in `commercial/data/resolution_sweep/` (one CSV per clustering
      resolution, ranking the genes that separate each cluster from all the others). The object carries
      two usable cell-type annotations: **`cell_type_res1_1` is the headline annotation** (18 clusters,
      15 distinct labels), and `cell_type_res06` is a coarser alternative (10 clusters, 10 labels). It
      also carries a third column, `cell_type_res08`, which is **withdrawn — do not deposit it as-is**
      (see the next TODO). The old `deg_leiden_wilcoxon.csv` still exists (now at `commercial/data/archive/`) but is
      **superseded**: it came out of an old exploratory notebook rather than the current pipeline, and
      the `resolution_sweep/` tables are its modern equivalent. It is marked for deletion but kept for
      now, because regenerating it would mean a destructive re-run of the pipeline — so it is not a
      deposit file. The six Fig2 sections are flagged on the `main` tab. Needed because Leiden cluster
      numbering is non-deterministic, so the annotations cannot be regenerated without the object.
- [ ] **Drop the stale `cell_type_res08` column before GEO deposit** — the object
      `commercial/data/objects/combined_clustered.h5ad` carries a third cell-type column,
      `cell_type_res08`, whose 14 labels are wrong. The names themselves were never careless: they were
      hand-curated from marker genes in the old notebook
      `commercial/notebooks/archive/ALL_MM_01July2026_fixed_v2.ipynb` (cell 17), against that notebook's own
      clustering at Leiden resolution 0.8. ("Leiden" is the algorithm that sorts cells into groups;
      "resolution" is a dial — turn it up and you get more, smaller groups.) What went wrong is that the
      notebook clustered a **different set of cells**: it normalized before filtering, which made its
      count filter do nothing and left 251,444 cells, where the current `02b_qc_filtering.py` filters on
      raw counts first and leaves 214,721. Those 36,723 extra cells changed the grouping — the notebook
      found **18** clusters at resolution 0.8, the pipeline finds **14**; `commercial/README.md`
      explains the ordering bug in full. The archived script
      `commercial/scripts/archive/04b_labeling_res08.py` then applied the notebook's 18 names to the
      pipeline's 14 clusters **by number**, so the name curated for the notebook's cluster 3 landed on a
      pipeline cluster 3 built from different cells, and names 14–17 matched nothing and were dropped
      without a word. `02b`'s ordering is the correct one — the pipeline fixed a real bug — so the fix
      is to rebuild the mapping against the current clustering, **never** to restore the old filter
      order. Nothing reads that column today — `06_umap_celltypes.py` and `07_spatial_celltypes.py` use only
      `cell_type_res1_1` — so no figure is affected. But this object is the file that gets uploaded to
      GEO, where a mislabeled column becomes a public
      artifact someone else may trust. **Either delete the column before deposit, or confirm it is
      harmless to ship.**
- [ ] **Pre-submission GEO logistics** — email GEO about the > 1 TB total before submitting; decide
      compression/`.gz` naming and finalize deposited filenames; optionally fill the `description`
      column (left blank).

## Setup & usage

<!-- FILL: Python deps, HISE access requirements, typical run order. -->
