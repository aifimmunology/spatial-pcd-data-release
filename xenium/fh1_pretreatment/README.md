# fh1_pretreatment — FH1 cohort, MM Pre-Treatment samples

The 13 single-patient samples from the non-commercial FH1 cohort taken at the MM
Pre-Treatment timepoint. These are the samples behind paper Figure 3. The tracking sheet's
`downstream analysis (Fig3)` column flags 14; TIS04397-001 (patient FH1010) is left out because
it was run on a 377-gene panel while the other 13 used the 480-gene panel. Keeping it would
shrink the shared gene set to the 258 genes the two panels have in common.

The sample list, with HISE UUIDs and patient IDs, is at the top of `scripts/01_data_ingest.py`.

## The pipeline

Run from the repository root, in order:

```
pixi run python xenium/fh1_pretreatment/scripts/01_data_ingest.py
pixi run python xenium/fh1_pretreatment/scripts/02a_qc_calibration.py
pixi run python xenium/fh1_pretreatment/scripts/02b_qc_filtering.py
pixi run python xenium/fh1_pretreatment/scripts/03_processing.py
```

Each script is self-contained: it sets its own paths and parameters at the top and inlines its
own logic. There is no shared `config.py` or `_pipeline.py`. `02a_qc_calibration` only draws
diagnostic plots used to pick the filter cutoffs — it writes no data, so a plain re-run can
skip it.

| Stage | Reads | Writes |
|---|---|---|
| `01_data_ingest` | `data/samples/{sample}.h5ad` if present, otherwise downloads from HISE into `data/outs/` and saves `data/samples/{sample}.h5ad` | `data/objects/combined_raw.h5ad` |
| `02a_qc_calibration` | `combined_raw.h5ad` | QC calibration plots in `figures/qc_calibration/` (no data) |
| `02b_qc_filtering` | `combined_raw.h5ad` | `data/objects/combined_filtered.h5ad` |
| `03_processing` | `combined_filtered.h5ad` | `data/objects/combined_clustered.h5ad`, `data/resolution_sweep/deg_leiden_res*.csv`, plots in `figures/` |

**`02a_qc_calibration`** draws per-sample count/gene histograms and knee plots (plus a summary
table) from the raw object so the filter cutoffs can be read off the distributions. It writes
only figures.

**`02b_qc_filtering`** drops cells with fewer than 15 transcripts or fewer than 10 detected
genes, then drops the gene AZU1 — before normalizing, so AZU1's counts do not inflate any
cell's library size. AZU1 was removed because it behaved like background signal rather than
real expression (the evidence is in `figures/archive/lognorm_azu1_kept/qc_diagnostics/`).

**`03_processing`** builds the clustering. The steps, in order:
1. **Pearson residuals.** Each count is replaced by how far it sits from what you would
   expect for a cell with that many total transcripts. This is different from the usual
   "scale every cell to the same total, then take the log". The residuals are used only to
   build the embedding; the stored expression values stay log-normalized so marker tests and
   expression plots read normally.
2. **PCA, 30 components.** This compresses 479 genes into 30 summary axes.
3. **Harmony on `cohort_subject`.** Harmony is a batch correction that shifts each patient's
   cells so that patient identity alone does not drive the clustering.
4. **Neighbor graph and UMAP.** The UMAP uses `min_dist=0.5`, matching the archived Pearson run.
5. **Leiden clustering at 11 resolutions, 0.1 to 1.1.** Higher resolution gives more,
   smaller clusters.
6. **Wilcoxon marker-gene tables**, one per resolution.

The seed is 0 everywhere.

`data/objects/` does not exist yet. The pipeline has not been run since it was written, and
every earlier result is under `archive/`.

## archive/ — earlier runs (read-only)

Every earlier run keeps its own scripts, data and figures under the same name:
`scripts/archive/<run>/`, `data/archive/<run>/`, `figures/archive/<run>/`. Never write into
these folders.

**`lognorm_azu1_kept/`**: the first run.
- AZU1 is kept (480 genes).
- Normalization is the standard kind: each cell is scaled to the same total, then log-transformed.
- PCA runs on that log-normalized data.
- Holds the QC diagnostics (`04`) — the three-test screen that identified AZU1 as background —
  which live only here, not in the current pipeline. (The QC calibration step, `02a`, has since
  been restored to the live pipeline above.)

**`azu1_removed/`**: the same input with AZU1 dropped, clustered two ways.
- Standard normalization: `03_processing.py` → `combined_clustered.h5ad`.
- Pearson residuals: `05_normalization_compare.py` → `combined_pearson_clustered.h5ad`. This
  also compares the two clusterings; `06` summarizes that comparison on one page.
- `07` draws the figures for the Pearson clustering (`figures/archive/azu1_removed_pearson/`).
- `08` writes its marker tables (`data/archive/azu1_removed/resolution_sweep_pearson/`).

The current pipeline reproduces the **Pearson** branch of `azu1_removed`. It uses the same
samples, filter, AZU1 drop, normalization, PCA, Harmony, neighbor, UMAP and Leiden settings.
Its output still won't match that branch byte for byte, for three reasons:
- **The starting files.** It starts from the cached `data/samples/` files only if they are
  present. Those are the same files the archived run used.
- **What gets stored.** It keeps fewer extra fields in the `.h5ad`: no stored PCA gene
  weights, no `plasma_score`.
- **File names.** Its marker tables are named `deg_<key>.csv`, without the `_wilcoxon` suffix.

Cluster labels should agree. That is checked by comparing the two runs' labels with the
**adjusted Rand index**, a 0-to-1 score where 1 means every pair of cells is grouped the same
way in both runs and 0 means the agreement is no better than chance.

## notebooks/

`fh1_pretreatment.ipynb` is exploration on the `lognorm_azu1_kept` objects: re-embedding with
different UMAP settings and highlighting single clusters.
