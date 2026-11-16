# commercial — Xenium clustering and cell-type annotation pipeline

The six commercial-biobank sections (Accio, OHSU, TriMetis) flagged `yes` in the tracking
sheet's `downstream analysis (Fig2 D-F)` column. This folder contains a linear pipeline that
takes raw Xenium outputs from HISE through quality control, batch-corrected clustering, and
manual cell-type annotation, and ends by drawing the UMAP and spatial plots used in panels
D, E and F of paper Figure 2. Run the scripts in numbered order; each one reads a
file written by an earlier step.

There are nine scripts. The numbering skips `04b` — that script was withdrawn and now lives
in `archive/`. The "Archived" section at the end of this document explains why.

---

## Data flow

```
HISE — six sample ZIPs, UUIDs hardcoded in 01_data_ingest.py
    │
    └─ 01_data_ingest.py
         ├─ xenium/commercial/data/outs/{sample}/          (extracted Xenium output folders)
         ├─ xenium/commercial/data/samples/TIS05393.h5ad
         ├─ xenium/commercial/data/samples/TIS04788.h5ad
         ├─ xenium/commercial/data/samples/TIS05392.h5ad
         ├─ xenium/commercial/data/samples/TIS05394.h5ad
         ├─ xenium/commercial/data/samples/TIS04789.h5ad
         ├─ xenium/commercial/data/samples/TIS05397_008.h5ad
         └─ xenium/commercial/data/objects/combined_raw.h5ad
              │
              ├─ 02a_qc_calibration.py     (reads samples/*.h5ad; draws plots, writes no h5ad)
              │    ├─ xenium/commercial/figures/qc_calibration/qc_distributions.{png,pdf}
              │    └─ xenium/commercial/figures/qc_calibration/spatial_qc_maps.{png,pdf}
              │
              └─ 02b_qc_filtering.py
                   └─ xenium/commercial/data/objects/combined_filtered.h5ad
                        │
                        └─ 03_processing.py
                             ├─ xenium/commercial/data/objects/combined_clustered.h5ad
                             ├─ xenium/commercial/data/resolution_sweep/
                             │      deg_leiden_res{R}_wilcoxon.csv          (11 files)
                             ├─ xenium/commercial/figures/processing/umap_harmony_diagnostic.{png,pdf}
                             └─ xenium/commercial/figures/processing/umap_resolution_sweep.{png,pdf}
                                  │
                                  ├─ 04a_labeling_res06.py
                                  │    └─ combined_clustered.h5ad  (adds column cell_type_res06;
                                  │                                 nothing downstream reads it)
                                  │
                                  └─ 04c_labeling_res1_1.py
                                       └─ combined_clustered.h5ad  (adds column cell_type_res1_1;
                                                                    this is the headline annotation)
                                            │
                                            ├─ 05_umap_by_sample.py
                                            │    └─ xenium/commercial/figures/umap_by_sample/
                                            │           umap_by_sample.{png,pdf}
                                            │
                                            ├─ 06_umap_celltypes.py
                                            │    ├─ xenium/commercial/figures/umap_celltypes/
                                            │    │      umap_celltypes_res1_1_ondata.{png,pdf}
                                            │    │      umap_celltypes_res1_1_legend.{png,pdf}
                                            │    └─ combined_clustered.h5ad  (rewritten — see below)
                                            │
                                            └─ 07_spatial_celltypes.py
                                                 └─ xenium/commercial/figures/spatial_celltypes/
                                                        spatial_celltypes_res1_1_{sample}.{png,pdf}
                                                        (6 samples)
```

All paths are relative to the repository root (`/home/workspace/PCD_data`), except inside the
scripts themselves, which use absolute paths built from two constants near the top of each
file: `DATA_DIR` (`/home/workspace/PCD_data/xenium/commercial/data`) and `OBJECTS_DIR`
(`DATA_DIR/objects`).

**06 must run before 07.** This is not just a convention — the two scripts share data through
the h5ad file, and running them out of order silently produces mismatched figure colors. The
`06_umap_celltypes.py` section below explains the mechanism.

---

## Terms used in this document

These words appear throughout the script descriptions. Each is defined once here.

**AnnData and `.h5ad`.** AnnData is the Python object that holds a single-cell dataset: one
big table of numbers plus everything you know about its rows and columns. `.h5ad` is the file
format that object is saved to. One `.h5ad` file on disk equals one AnnData object in memory.

**The parts of an AnnData object.** The pipeline writes into four of them, and it matters
which, because different parts survive different operations:

- `.X` — the main matrix. One row per cell, one column per gene. Here that is 214,721 rows
  by 377 columns after filtering.
- `.obs` — a spreadsheet with one row per cell, holding per-cell facts. Columns include
  `sample` (which tissue the cell came from), `total_counts`, and every cluster assignment
  and cell-type label the pipeline produces.
- `.obsm` — extra matrices with one row per cell but many columns. Holds `X_pca`, `X_umap`
  and `spatial` (each cell's physical x, y position on the slide).
- `.layers` — alternative versions of `.X` with the same shape. `layers["counts"]` holds the
  original whole-number transcript counts, kept because `.X` itself gets overwritten by
  normalization.
- `.uns` — a free-form bag of everything else. Relevant here because scanpy stores the
  **color palette** for each categorical `.obs` column in `.uns`, under a key like
  `cell_type_res1_1_colors`.

**Counts.** A count is the number of transcripts of one gene detected in one cell. Counts are
whole numbers: 0, 1, 47 — never 2.7.

**Library size and normalization.** A cell's library size is its total count across all 377
genes. Cells differ enormously in library size for reasons that are technical rather than
biological: in `combined_raw.h5ad` the median cell in TIS04789 has 12 total counts while the
median cell in TIS05397_008 has 31. Normalization rescales every cell so that its total comes
out the same, so that "this cell has more of gene X" means more *relative to that cell's own
output*, not just that the cell was captured better.

**PCA (principal component analysis).** A way of summarizing 377 gene measurements per cell
with a much smaller set of numbers — here 30. Each of the 30 is a weighted blend of the
original genes, chosen so the blends capture as much of the cell-to-cell variation as
possible. The point is speed and noise reduction: the later steps work on 30 numbers per
cell instead of 377.

**Batch effect and Harmony.** A batch effect is a systematic difference between samples that
comes from how they were processed, not from the biology. If uncorrected, cells cluster by
which slide they came from rather than by what kind of cell they are. **Harmony** is an
algorithm that takes the 30 PCA numbers per cell and nudges them so that cells of the same
biological type from different samples land on top of each other, while genuinely different
cell types stay apart. You tell it which `.obs` column names the batches; here that column is
`sample`, so the six tissue IDs are the six batches.

**Neighbor graph.** A structure that records, for each cell, which other cells are most
similar to it. Built here from the 15 nearest neighbors of each cell in the Harmony-corrected
30-dimensional space. Both UMAP and Leiden read this graph rather than the raw numbers.

**UMAP.** A method that places each cell at an (x, y) position on a flat plot so that cells
that are neighbors in the graph land near each other. It is a picture, not a measurement:
distances on a UMAP are not quantitative, and the axes have no units. Every scatter plot in
Figures D and E is a UMAP.

**Leiden clustering.** An algorithm that cuts the neighbor graph into groups, so that cells
within a group are densely connected to each other and sparsely connected to cells outside
it. Each group gets a number — 0, 1, 2, and so on. The numbers are arbitrary labels with no
ordering; cluster 3 is not "between" clusters 2 and 4.

**Resolution.** The single dial that controls how finely Leiden cuts. A low resolution
produces few large clusters; a high resolution produces many small ones. In this dataset
resolution 0.1 gives 4 clusters and resolution 1.1 gives 18. There is no correct value — you
pick one by looking at whether the resulting clusters make biological sense.

**Cluster numbers are not stable across resolutions.** Cluster 5 at resolution 0.8 and
cluster 5 at resolution 1.1 are different sets of cells that happen to share a number. This
is the single most important fact about the annotation scripts, and ignoring it is what
retired `04b` (see "Archived").

**Marker gene.** A gene that is expressed much more strongly in one cluster than in the rest
of the cells. Marker genes are how you decide what a cluster *is* — a cluster whose top
markers are hemoglobin genes is erythroid.

**Wilcoxon rank-sum test.** The statistical test used to find marker genes. For one gene and
one cluster, it ranks every cell in the dataset by how much of that gene it expresses, then
asks whether the cluster's cells sit systematically higher in the ranking than everyone else.
It works on ranks rather than on the count values themselves, which means a handful of cells
with enormous counts cannot single-handedly declare a marker. The output per gene–cluster
pair is a score, a p-value, and a log fold change (how many doublings separate the cluster's
average from everyone else's).

**Knee plot.** A plot used in `02a` to choose a count threshold. Sort all cells from highest
total count to lowest, then plot count against rank position, both axes logarithmic. Real
cells form a high plateau on the left; debris and empty captures form a low tail on the
right. The bend between them — the "knee" — is a candidate cutoff.

---

## Figure specifications

All three figures are drawn from the same file, `combined_clustered.h5ad`, and all three are
saved twice: a `.png` at 150 dots per inch for quick viewing, and a `.pdf` for the manuscript.

In every figure the scatter points are **rasterized**: the points are baked into the PDF as a
pixel image while the axes, ticks and text stay as editable vector objects. Without this, a
PDF holding 214,721 individually drawn circles is slow to open and awkward to edit.

### UMAP colored by sample (paper Figure 2D)

One UMAP of all 214,721 cells, each point colored by which of the six tissues it came from
(the `sample` column in `.obs`). The color legend sits in the right margin, outside the plot
axes. This panel exists to show that Harmony worked: if cells separated into six
non-overlapping islands by color, the batch correction would have failed.

### UMAP colored by cell type (paper Figure 2E)

The same UMAP, colored by the `cell_type_res1_1` annotation (15 labels). Two versions are
saved:

1. **Labels on data** (`..._ondata`) — each cell-type name is drawn directly on the plot at
   the center of its cells, at 8-point type. No separate legend.
2. **Legend in margin** (`..._legend`) — nothing written on the plot; the color key sits in
   the right margin, outside the axes.

Only the resolution-1.1 annotation is drawn. `cell_type_res06` also exists in the file but is
not plotted.

### Per-tissue spatial maps colored by cell type (paper Figure 2F)

Six plots, one per tissue. Each plot shows that tissue's cells at their real physical
positions on the slide — the (x, y) micron coordinates stored in `.obsm['spatial']` — colored
by the same `cell_type_res1_1` annotation used in Figure E. The legend sits in the right
margin.

The point of this panel is that Figures D and E show which cell types exist, while Figure F
shows *where they sit relative to each other in the tissue*.

---

## How to run

```bash
cd /home/workspace/PCD_data
python xenium/commercial/scripts/01_data_ingest.py
python xenium/commercial/scripts/02a_qc_calibration.py   # inspect plots, then set thresholds in 02b
python xenium/commercial/scripts/02b_qc_filtering.py
python xenium/commercial/scripts/03_processing.py
python xenium/commercial/scripts/04a_labeling_res06.py
python xenium/commercial/scripts/04c_labeling_res1_1.py
python xenium/commercial/scripts/05_umap_by_sample.py
python xenium/commercial/scripts/06_umap_celltypes.py
python xenium/commercial/scripts/07_spatial_celltypes.py
```

Three things to know before running the whole chain:

- **`02a` is an inspection step.** It writes no data file. Run it, look at the plots, and
  decide whether the thresholds hardcoded in `02b` are right before continuing.
- **`04a` is optional.** Its output column feeds no figure. Skipping it changes nothing about
  Figures D, E or F.
- **`06` must run before `07`.** See the `06_umap_celltypes.py` section.

`archive/04b_labeling_res08.py` is deliberately absent from this list. Do not run it.

---

## Scripts

### `01_data_ingest.py` — Download, read, and concatenate

Downloads each Xenium sample as a ZIP from HISE, unpacks it, reads the unpacked folder with
`spatialdata_io.xenium`, and saves each sample as its own `.h5ad`. Then stacks all six into
one object and saves that too.

**Inputs:** a `XENIUM_SAMPLES` dictionary hardcoded at the top of the script, mapping sample
name to HISE file UUID:

| Sample | HISE UUID | Tissue |
|---|---|---|
| TIS05393 | `70f7c0d6-8d7e-4b86-acf3-ca80137fd849` | Multiple Myeloma |
| TIS04788 | `38fe47a0-c09f-4373-9a8a-8a957c40aa9f` | Multiple Myeloma |
| TIS05392 | `0c879d26-970f-47fc-a705-88b103b4de80` | Multiple Myeloma |
| TIS05394 | `d7b4b14e-9856-4169-a363-213b4e29a99e` | Multiple Myeloma |
| TIS04789 | `7c3113c3-23bf-4840-8f45-61505048322c` | Multiple Myeloma |
| TIS05397_008 | `16752455-bd23-4f88-8d05-2108e7778496` | Normal bone marrow (2024) |

**Outputs:**

- `xenium/commercial/data/outs/{sample}/` — the unpacked Xenium output folder for each
  sample, straight from the ZIP.
- `xenium/commercial/data/samples/{sample}.h5ad` — six files, one per sample.
- `xenium/commercial/data/objects/combined_raw.h5ad` — all six stacked into one object,
  399,545 cells by 377 genes. A `sample` column is added to `.obs` recording which tissue
  each cell came from.

Per-sample cell counts in the combined raw object:

| Sample | Cells | Median counts per cell |
|---|---:|---:|
| TIS05393 | 114,863 | 15 |
| TIS04788 | 49,534 | 17 |
| TIS05392 | 155,811 | 18 |
| TIS05394 | 17,491 | 17 |
| TIS04789 | 38,748 | 12 |
| TIS05397_008 | 23,098 | 31 |
| **Total** | **399,545** | 16 |

**Parameters to configure:** the UUIDs in `XENIUM_SAMPLES`. Change these only if samples are
added or swapped.

**Notes:**

The samples are stacked with `join='inner'`, meaning only genes present in all six samples
are kept. All six use the same 377-gene panel, so nothing is dropped in practice.

The download step is skipped for any sample whose folder already exists under `outs/`. The
HISE fetch itself (`hisepy.cache_files`) is only skipped if `hisepy` has the ZIP in its own
cache.

**`outs/` is currently empty.** The extracted bundles came to roughly 27 GB and were deleted,
because everything downstream reads from `samples/*.h5ad` instead. Nothing in the pipeline
breaks. If you re-run this script it will simply re-download all six from HISE, which takes a
while — so only re-run `01` if you actually need the raw bundles back or a sample list has
changed.

---

### `02a_qc_calibration.py` — QC inspection (writes no data)

Reads the six per-sample files and draws plots that let you choose filtering thresholds. It
does not modify or write any `.h5ad`.

**Inputs:** `xenium/commercial/data/samples/{sample}.h5ad` (all six).

**Outputs:**

- `xenium/commercial/figures/qc_calibration/qc_distributions.{png,pdf}` — six rows (one per sample) by three columns:
  a histogram of log10(total counts), a histogram of log10(genes detected), and a knee plot.
- `xenium/commercial/figures/qc_calibration/spatial_qc_maps.{png,pdf}` — six rows by two columns, each panel showing
  the sample's cells at their physical slide positions, colored by `total_counts` or by
  `n_genes_by_counts`. The color scale is capped at the 99th percentile (`vmax='p99'`) so
  that a few extreme cells do not flatten the whole color range.
- A summary table printed to the terminal: cell count, median counts, and median genes per
  sample.

**Parameters to configure:** none. This script is deliberately passive — it shows you the
data so you can pick numbers for `02b`.

**Notes:**

`n_genes_by_counts` is recomputed by hand here as the number of genes with at least one
count in the cell, rather than being taken from scanpy's `calculate_qc_metrics` output.

Read the knee plot to choose `MIN_COUNTS`. Read the spatial maps to check *where* the
low-count cells are. Low-count cells forming a rim at the tissue edge are normal in Xenium —
those cells were sliced through. Low-count cells scattered evenly across the whole tissue
would instead point at a problem with the run itself.

---

### `02b_qc_filtering.py` — Cell filtering and normalization

Drops cells that fall below the count and gene thresholds, then rescales the survivors to a
common library size.

**Inputs:** `xenium/commercial/data/objects/combined_raw.h5ad`

**Outputs:** `xenium/commercial/data/objects/combined_filtered.h5ad` — 214,721 cells by 377
genes, with `.X` normalized and the original whole-number counts preserved in
`.layers['counts']`.

**Parameters to configure:**

- `MIN_COUNTS` (currently `15`) — drop any cell with fewer than 15 total transcript counts.
- `MIN_GENES` (currently `10`) — drop any cell with fewer than 10 distinct genes detected.

**Notes — how much this removes, and why it is so much:**

With the current thresholds the filter removes **184,824 of 399,545 cells, or 46.3%**. Per
sample:

| Sample | Before | After | Removed | % removed |
|---|---:|---:|---:|---:|
| TIS05393 | 114,863 | 56,889 | 57,974 | 50.5% |
| TIS04788 | 49,534 | 26,971 | 22,563 | 45.6% |
| TIS05392 | 155,811 | 87,165 | 68,646 | 44.1% |
| TIS05394 | 17,491 | 9,719 | 7,772 | 44.4% |
| TIS04789 | 38,748 | 15,908 | 22,840 | 58.9% |
| TIS05397_008 | 23,098 | 18,069 | 5,029 | 21.8% |
| **Total** | **399,545** | **214,721** | **184,824** | **46.3%** |

Losing nearly half the cells looks alarming until you line the threshold up against the data.
The median cell in the raw combined object has **16** total counts. `MIN_COUNTS` is **15**.
The cutoff sits essentially at the middle of the distribution, so removing about half the
cells is the arithmetic consequence of where the threshold was placed, not a sign that
something went wrong.

The spread across samples follows the same logic. TIS04789's median cell has 12 counts, below
the threshold, so it loses the most (58.9%). TIS05397_008's median cell has 31 counts, twice
the threshold, so it loses the least (21.8%).

The practical consequence: if you want to keep more cells, lower `MIN_COUNTS`, and expect the
gain to be uneven — TIS04789 will recover far more cells than TIS05397_008.

**Notes — normalization and the counts layer:**

The raw matrix is copied into `.layers['counts']` *before* normalization, on line 28. This
matters because `sc.pp.normalize_total` overwrites `.X` in place. After this script,
`.X` holds rescaled fractional values and `.layers['counts']` holds the original whole
numbers. Any method that requires whole-number input — a negative-binomial model, for
instance — must read the layer, not `.X`.

`normalize_total` is called with no `target_sum`, which means scanpy uses its default: every
cell is rescaled so its total equals the median library size across the dataset.

Finally, `.X` is converted to a sparse matrix if it is not one already. A sparse matrix stores
only the non-zero entries. With 377 genes and a median of 16 counts per cell, the great
majority of entries are zero, so this is a large saving in file size.

---

### `03_processing.py` — PCA, Harmony, UMAP, and the resolution sweep

The computational core. Reduces the data with PCA, removes batch effects with Harmony, builds
a neighbor graph, computes a UMAP, and then clusters at eleven different resolutions, writing
marker genes for each.

**Inputs:** `xenium/commercial/data/objects/combined_filtered.h5ad`

**Outputs:**

- `xenium/commercial/data/objects/combined_clustered.h5ad` — the filtered object plus
  `.obsm['X_pca']`, `.obsm['X_pca_harmony']`, `.obsm['X_umap']`, and eleven
  `leiden_res{R}` columns in `.obs`.
- `xenium/commercial/data/resolution_sweep/deg_leiden_res{R}_wilcoxon.csv` — **eleven** files,
  one per resolution. Each row is one gene–cluster pair, carrying the Wilcoxon score, the
  p-value, the adjusted p-value, and the log fold change.
- `xenium/commercial/figures/processing/umap_harmony_diagnostic.{png,pdf}` — three UMAP panels side by side,
  colored by `total_counts`, `n_genes_by_counts`, and `sample`.
- `xenium/commercial/figures/processing/umap_resolution_sweep.{png,pdf}` — a grid three panels wide, one UMAP per
  resolution, with cluster numbers printed on the plot.

**Parameters to configure:**

- `RESOLUTIONS` — currently `[round(0.1 * i, 2) for i in range(1, 12)]`, which is 0.10, 0.20,
  … , 1.10. That is **eleven** values, and is why eleven CSVs are written.
- `n_comps=30` in the PCA call. Also passed `use_highly_variable=False`, meaning all 377
  genes contribute — reasonable here, because a targeted 377-gene panel was already curated
  by hand and there is no long tail of uninformative genes to discard.
- `n_neighbors=15` in the neighbor graph.
- `random_state=0` on PCA, Harmony, neighbors, UMAP and Leiden, so repeated runs on the same
  input give the same answer.

**Notes — the `X_pca` overwrite:**

After Harmony finishes, the corrected coordinates are written to **two** places (lines 42–43):

```python
combined.obsm['X_pca_harmony'] = z
combined.obsm['X_pca'] = z  # downstream steps use X_pca by default
```

The second line is the load-bearing one. Scanpy's `neighbors`, `umap` and `leiden` all reach
for `X_pca` unless told otherwise. By overwriting `X_pca` with the corrected values, every
later step uses the batch-corrected embedding automatically, with no further arguments. The
consequence to be aware of: after this script, `.obsm['X_pca']` does **not** contain raw PCA
output. The uncorrected coordinates are not kept anywhere.

There is also a shape guard just above it. Different versions of `harmonypy` return the
corrected matrix either as cells-by-components or components-by-cells; the script checks
which and transposes if needed.

**Notes — the diagnostic UMAP:**

`umap_harmony_diagnostic` is how you check Harmony actually worked. Look at the `sample`
panel: the six colors should be mixed throughout the plot. Six separate single-colored blobs
would mean the correction did not take, and every cluster downstream would be a sample rather
than a cell type. The `total_counts` and `n_genes_by_counts` panels catch a different failure
— if one region of the UMAP is uniformly low on both, that region is probably poor-quality
cells grouping together rather than a real cell population.

**Notes — how many clusters each resolution gives:**

| Resolution | Clusters |
|---:|---:|
| 0.10 | 4 |
| 0.20 | 6 |
| 0.30 | 5 |
| 0.40 | 8 |
| 0.50 | 8 |
| 0.60 | 10 |
| 0.70 | 13 |
| 0.80 | 14 |
| 0.90 | 15 |
| 1.00 | 16 |
| 1.10 | 18 |

Cluster count generally rises with resolution but not strictly: resolution 0.30 gives **5**
clusters, one *fewer* than resolution 0.20's 6. That is expected behavior, not a bug. Leiden
re-partitions the graph from scratch at each resolution rather than subdividing the previous
answer, so a small change in the dial can reshuffle boundaries rather than just split them.
The practical consequence is the one stated in the Terms section: cluster 3 at one resolution
has nothing to do with cluster 3 at another.

**Notes — choosing a resolution:**

Use the sweep grid and the DEG CSVs together. The grid shows you how the clusters carve up
the UMAP; the CSVs tell you what genes distinguish them. A resolution is too high when
splitting a cluster in two yields two groups with no distinguishing markers between them.
Resolution 1.1 was the one chosen for the figures (see `04c`).

**Notes — what is currently on disk does not match this script:**

`combined_clustered.h5ad` and the `resolution_sweep/` folder were produced by an **earlier
version** of this script whose `RESOLUTIONS` ran to 1.5. As a result the folder holds fifteen
CSVs (`res0.10` through `res1.50`), and the object carries fifteen `leiden_res{R}` columns
rather than eleven. The script as it stands today writes eleven of each.

The practical consequence: re-running `03_processing.py` will overwrite the object and drop
the `leiden_res1.20` through `leiden_res1.50` columns. The four surplus CSVs would be left
behind in the folder, since nothing deletes them. `leiden_res1.10` — the column `04c` needs —
is inside the current range, so the annotation step is unaffected either way.

---

### `04a_labeling_res06.py` — Cell-type annotation at resolution 0.6

Turns the cluster numbers in `leiden_res0.60` into readable cell-type names and saves them as
a new `.obs` column. Writes no figure.

**Inputs:** `xenium/commercial/data/objects/combined_clustered.h5ad`

**Outputs:** the same file, rewritten with a `cell_type_res06` column added to `.obs`.

**Parameters to configure:** `cluster_to_cell_type_res06`, a dictionary from cluster-number
string to label. `KEY` is set to `'leiden_res0.60'`.

| Cluster | Cell type | Cells |
|---|---|---:|
| 0 | Erythroid Precursor 1 | 17,214 |
| 1 | Erythroid Precursor 2 | 15,484 |
| 2 | Neoplastic Plasma Cell 1 (MYC/CXCR4-high) | 8,447 |
| 3 | Neoplastic Plasma Cell 2 (LILRB4-high) | 12,836 |
| 4 | Endothelial | 30,097 |
| 5 | Stromal Fibroblast | 12,599 |
| 6 | Neoplastic Plasma Cell 3 | 61,804 |
| 7 | Cytotoxic T Cell | 6,123 |
| 8 | Macrophage | 11,914 |
| 9 | Monocyte | 38,203 |
| 10 | Proliferating Myeloid Precursor | **never applied — see below** |

**Notes — this annotation is not used by any figure:**

`cell_type_res06` is written into the h5ad and then read by nothing. `06_umap_celltypes.py` and
`07_spatial_celltypes.py` both plot `cell_type_res1_1` only. Treat this script as a live alternative
annotation kept for comparison, not as a step that feeds Figure 2. Skipping it does not change
any figure.

**Notes — the mapping has 11 entries but the clustering has 10 clusters:**

`leiden_res0.60` contains clusters `"0"` through `"9"` — ten of them. The dictionary above has
eleven keys, `"0"` through `"10"`. The surplus key `"10"` ("Proliferating Myeloid Precursor")
matches no cluster and is therefore never used. `cell_type_res06` ends up with ten labels, and
"Proliferating Myeloid Precursor" is not among them.

The script does not warn you about this. Its guard is:

```python
missing = set(combined.obs[KEY].unique()) - set(cluster_to_cell_type_res06.keys())
```

That subtraction is *clusters in the data* minus *keys in the mapping*. It catches the case
where the data has a cluster the mapping forgot. It cannot catch the reverse case — extra keys
in the mapping — because extra keys never appear in that set. With ten clusters and eleven
keys the result is empty, so no error is raised.

**What is affected, concretely.** Exactly one key is affected: `"10"`. It is the last key in
the dictionary, it matches no cluster, and its label "Proliferating Myeloid Precursor" appears
nowhere in `cell_type_res06`. Clusters `"0"` through `"9"` — all ten that exist — are not
affected. Each received the label filed under its own number, which is verified by the cell
counts in the table above: `leiden_res0.60` cluster 0 holds 17,214 cells and
`cell_type_res06` has exactly 17,214 cells labeled Erythroid Precursor 1, and the same one-to-one
match holds for all ten clusters. No cell carries a label belonging to a different cluster.

**So this is a milder problem than the one that retired `04b`, and it is worth being precise
about the difference.** In `04b`, four of eighteen keys went unused and the fourteen that
applied came from a different clustering, so the labels landed on the wrong cells. Here, no
label landed on the wrong cells. One label simply never landed at all.

**What it still tells you.** A mapping written with eleven entries was curated against a
clustering that had eleven groups. This object's resolution 0.6 has ten. So the mapping and
the clustering it is being applied to did not come from the same run, even though the ten
labels that survived are internally consistent. The practical consequence: treat
`cell_type_res06` as unverified until the mapping is rebuilt by reading the marker genes of
the current ten-cluster `leiden_res0.60` run, and in particular do not assume the missing
Proliferating Myeloid Precursor population is absent from the data — at resolution 1.1 it is a
real population of 6,479 cells (cluster 8). At resolution 0.6 it is presumably merged into one
of the ten clusters, unnamed.

Because nothing downstream reads this column, none of this affects Figures D, E or F. It is a
blocker only if you decide to use `cell_type_res06` for something.

**Warning:** this mapping is valid only for the exact clustering stored in the current
`combined_clustered.h5ad`. If `03_processing.py` is re-run, cluster numbers can move, and the
dictionary must be rebuilt by reading the new marker genes from scratch.

---

### `04c_labeling_res1_1.py` — Cell-type annotation at resolution 1.1 (the headline annotation)

**This is the annotation used in Figures E and F.** Same shape as `04a`, but for
`leiden_res1.10`, which has 18 clusters. Writes no figure.

**Inputs:** `xenium/commercial/data/objects/combined_clustered.h5ad`

**Outputs:** the same file, rewritten with a `cell_type_res1_1` column added to `.obs`.

**Parameters to configure:** `cluster_to_cell_type_res1_1`. `KEY` is set to
`'leiden_res1.10'`. The mapping has 18 entries, matching the 18 clusters exactly.

| Cluster | Cell type | Cells in cluster |
|---|---|---:|
| 0 | Erythroid Precursor 1 | 18,598 |
| 1 | Erythroid Precursor 2 | 3,408 |
| 2 | Erythroid Precursor 2 | 11,882 |
| 3 | Neoplastic Plasma Cell 2 (BCL2L11-high, apoptosis-primed) | 6,790 |
| 4 | Neoplastic Plasma Cell 3 | 10,133 |
| 5 | Neoplastic Plasma Cell 6 | 31,841 |
| 6 | Stromal Fibroblast | 9,521 |
| 7 | Neoplastic Plasma Cell 2 (BCL2L11-high, apoptosis-primed) | 32,373 |
| 8 | Proliferating Myeloid Precursor | 6,479 |
| 9 | Macrophage | 7,133 |
| 10 | Monocyte | 12,483 |
| 11 | Macrophage | 9,778 |
| 12 | Cytotoxic T Cell | 3,840 |
| 13 | Neoplastic Plasma Cell 5 (MYC-high, possibly proliferative) | 817 |
| 14 | Neoplastic Plasma Cell 4 | 14,400 |
| 15 | Neoplastic Plasma Cell 1 | 27,681 |
| 16 | Endothelial | 7,370 |
| 17 | TREM2+ Macrophage (Tumor-Supportive, Adjacent) | 194 |

**Notes — 18 clusters, 15 labels:**

Three label names are deliberately reused, so the 18 clusters collapse into 15 distinct cell
types. Clusters 1 and 2 are both called Erythroid Precursor 2; clusters 3 and 7 are both
Neoplastic Plasma Cell 2; clusters 9 and 11 are both Macrophage. This is intentional. At
resolution 1.1 Leiden splits some populations more finely than the marker genes justify, and
merging them at the labeling step is how that is undone without dropping to a coarser
resolution that would under-split something else.

The consequence for the figures: Figures E and F show **15** colors, not 18. The resulting
per-label totals are:

| Cell type | Cells |
|---|---:|
| Neoplastic Plasma Cell 2 (BCL2L11-high, apoptosis-primed) | 39,163 |
| Neoplastic Plasma Cell 6 | 31,841 |
| Neoplastic Plasma Cell 1 | 27,681 |
| Erythroid Precursor 1 | 18,598 |
| Macrophage | 16,911 |
| Erythroid Precursor 2 | 15,290 |
| Neoplastic Plasma Cell 4 | 14,400 |
| Monocyte | 12,483 |
| Neoplastic Plasma Cell 3 | 10,133 |
| Stromal Fibroblast | 9,521 |
| Endothelial | 7,370 |
| Proliferating Myeloid Precursor | 6,479 |
| Cytotoxic T Cell | 3,840 |
| Neoplastic Plasma Cell 5 (MYC-high, possibly proliferative) | 817 |
| TREM2+ Macrophage (Tumor-Supportive, Adjacent) | 194 |

Two labels are very small. TREM2+ Macrophage has 194 cells out of 214,721, which is 0.09% of
the dataset, and Neoplastic Plasma Cell 5 has 817, or 0.38%. Both will be nearly invisible on
the figures unless deliberately highlighted, and both are small enough that any claim resting
on them should be checked against the marker genes directly rather than read off the plot.

**Notes — this mapping was checked and is sound:**

This is the annotation behind the published figures, so its correctness was verified directly
against the object rather than assumed from the fact that 18 equals 18. Matching counts alone
would not be enough: eighteen keys could still miss a cluster and duplicate another. The
actual check was on the key *sets*, in both directions:

- The mapping's keys are exactly `"0"`, `"1"`, … , `"17"` — eighteen of them, no gaps, no
  duplicates, no non-numeric entries.
- `leiden_res1.10`'s clusters are exactly `"0"` through `"17"` — the same eighteen.
- Clusters minus keys is empty, so no cluster went unlabeled. This is the direction the
  script's own guard checks.
- **Keys minus clusters is also empty**, so no label was silently dropped. This is the
  direction the guard cannot see, and it is what went wrong in `04b`.

The two sets are identical. Every cluster received a label and every label was used.

The guard itself still has the same blind spot as `04a`'s — it is the same
clusters-minus-keys subtraction, and it would not warn you about a surplus key. It simply has
nothing to miss here. If you edit this dictionary, re-check the keys-minus-clusters direction
by hand, because nothing in the script will do it for you.

**Warning:** as with `04a`, this mapping is tied to the exact clustering in the current
`combined_clustered.h5ad` and must be rebuilt from marker genes if `03_processing.py` is
re-run.

---

### `05_umap_by_sample.py` — UMAP colored by sample (Figure 2D)

**Inputs:** `xenium/commercial/data/objects/combined_clustered.h5ad`

**Outputs:**

- `xenium/commercial/figures/umap_by_sample/umap_by_sample.png`
- `xenium/commercial/figures/umap_by_sample/umap_by_sample.pdf`

**Parameters to configure:** point size (`size=6`) and output resolution (`dpi=150`) are
hardcoded in the plotting call.

**Notes:** the legend is placed with `legend_loc='right margin'`, so the six sample names sit
outside the axes. After the figure is built, every scatter collection is marked rasterized
before saving, which keeps the PDF light while leaving the axis labels and legend text
editable in a vector editor. This script only reads the h5ad; it does not write it back.

---

### `06_umap_celltypes.py` — UMAP colored by cell type (Figure 2E)

**Inputs:** `xenium/commercial/data/objects/combined_clustered.h5ad`

**Outputs:**

- `xenium/commercial/figures/umap_celltypes/umap_celltypes_res1_1_ondata.{png,pdf}`
- `xenium/commercial/figures/umap_celltypes/umap_celltypes_res1_1_legend.{png,pdf}`
- `xenium/commercial/data/objects/combined_clustered.h5ad` — **rewritten**, see below.

**Parameters to configure:** `ANNOTATION_COLS`, a list of (column name, filename tag) pairs.
It currently holds exactly one entry:

```python
ANNOTATION_COLS = [
    ('cell_type_res1_1', 'res1_1'),
]
```

The list is a loop variable, so adding `('cell_type_res06', 'res06')` to it would produce a
second pair of figures with `res06` in their filenames. Nothing else needs changing.

**Notes — why this script rewrites the h5ad, and why order matters:**

The last three lines of the script save `combined_clustered.h5ad` back to disk. No cell data
changed. What changed is the color palette, and the chain of reasoning runs like this:

1. When scanpy plots a categorical column for the first time, it has to pick a color for each
   category. It assigns them from a default palette and stores the result in `.uns`, under
   `cell_type_res1_1_colors`.
2. That assignment lives only in memory. When `06` exits, it is gone unless the object is
   written back to disk.
3. `07_spatial_celltypes.py` opens the same file and plots the same column. If the colors were never
   saved, scanpy assigns them again from scratch in a fresh process.
4. The colors chosen for a subset of categories need not match those chosen for the full set,
   and `07` plots one sample at a time. So Figure F's Macrophage could come out a different
   color from Figure E's Macrophage.
5. Saving in `06` fixes the palette in the file, so `07` reads the stored colors instead of
   inventing new ones.

**The practical consequence: run `06` before `07`.** If you run `07` on its own — against a
file where `06` has never run — the spatial maps and the UMAPs will use different colors for
the same cell types, and the figure panels will not read as one set. The failure is silent;
nothing errors, the colors are just wrong.

One further consequence of the write-back: because `06` rewrites the whole object, re-running
`04a` or `04c` afterwards will itself rewrite the file again and can reset the stored palette.
If you re-label, re-run `06` before `07`.

---

### `07_spatial_celltypes.py` — Per-tissue spatial maps colored by cell type (Figure 2F)

**Inputs:** `xenium/commercial/data/objects/combined_clustered.h5ad` — after `06` has run.

**Outputs:** `xenium/commercial/figures/spatial_celltypes/spatial_celltypes_res1_1_{sample}.{png,pdf}` — one
pair of files per sample, twelve files in total:

- `spatial_celltypes_res1_1_TIS04788.{png,pdf}`
- `spatial_celltypes_res1_1_TIS04789.{png,pdf}`
- `spatial_celltypes_res1_1_TIS05392.{png,pdf}`
- `spatial_celltypes_res1_1_TIS05393.{png,pdf}`
- `spatial_celltypes_res1_1_TIS05394.{png,pdf}`
- `spatial_celltypes_res1_1_TIS05397_008.{png,pdf}`

**Parameters to configure:** `ANNOTATION_COLS`, identical in form and content to `06`'s — a
single entry, `('cell_type_res1_1', 'res1_1')`. Keep the two lists in step: if you add an
annotation here without adding it to `06`, its colors will not have been saved and the
mismatch described above returns.

**Notes — the `X_spatial` copy:**

Each sample is pulled out of the combined object and given an extra entry in `.obsm`:

```python
sub.obsm['X_spatial'] = sub.obsm['spatial']
```

This is a naming workaround, not a transformation — the same coordinates under a second name.
`spatialdata_io` stores each cell's physical slide position under the key `spatial`, but
`sc.pl.embedding(basis='spatial')` looks for a key named `X_spatial`. Copying the array under
the expected name lets the plot find it. Nothing about the coordinates changes.

**Notes:** samples are iterated in sorted order, so the files come out alphabetically by
sample name, starting with TIS04788. Each sample is copied out with `.copy()` before plotting
so the added `X_spatial` key does not accumulate on the shared object.

---

## Archived

### `archive/04b_labeling_res08.py` — Cell-type annotation at resolution 0.8 (withdrawn)

This script is kept for the record and **should not be run**. It is excluded from the "How to
run" list above.

**What it did.** It set `KEY = 'leiden_res0.80'` and carried a mapping with **18** entries,
keyed `"0"` through `"17"`. It wrote the result into an `.obs` column called
`cell_type_res08`.

**The problem.** `leiden_res0.80` in `combined_clustered.h5ad` has **14** clusters, numbered
`"0"` through `"13"`. The mapping had eighteen keys for fourteen clusters.

**Why nothing warned about it.** The guard was the same subtraction used in `04a` and `04c`:

```python
missing = set(combined.obs[KEY].unique()) - set(mapping.keys())
```

That is *clusters present in the data* minus *keys present in the mapping*. It is built to
catch a cluster the mapping forgot. It cannot see keys the data has no cluster for, because
those keys never enter the left-hand set. With 14 clusters and 18 keys the subtraction gives
an empty set, no error fires, and the four surplus labels — Neoplastic Plasma Cell 6,
Lymphocyte/Dendritic Cell/Vessel, Proliferating Myeloid Precursor, and Lipid-Laden
Macrophage — are dropped without a word.

**Where the mapping actually came from.** It was not invented for this script. The identical
18-entry dictionary appears in `xenium/notebooks/archive/ALL_MM_01July2026_fixed_v2.ipynb`,
cell 17, with a comment stating it was curated by hand from the marker-gene plot drawn in the
cell above it. That preceding cell, cell 15, runs:

```python
sc.tl.leiden(combined, flavor='igraph', n_iterations=2, directed=False,
             resolution=0.8, random_state=0)
```

So the mapping is a genuine resolution-0.8 mapping, honestly curated from the marker genes of
a resolution-0.8 clustering that had eighteen clusters. It was then copied verbatim into
`ALL_MM_01July2026_fixed_v2_resolution_sweep.ipynb` cell 20, re-pointed at `leiden_res0.80`,
and from there into `04b`. Nobody rebuilt it along the way.

**Why the same resolution gives 18 clusters there and 14 here.** The Leiden call in the
notebook and the one in `03_processing.py` use identical settings — same flavor, same
iterations, same `random_state`, same resolution. The clusterings differ because the *cells
going in* differ, and that traces to one ordering difference upstream:

1. The notebook normalizes first and filters second. Cell 7 calls
   `sc.pp.normalize_total`; cell 10 then calls `sc.pp.filter_cells(min_counts=15)`.
2. `normalize_total` rescales every cell so its total equals the dataset median, which here is
   **16.000**. After that call, every cell's total is exactly 16.000 — verified by re-running
   the notebook's ordering on `combined_raw.h5ad`.
3. `filter_cells(min_counts=15)` reads the totals out of `.X`. Since 16.000 is above 15, the
   count filter has essentially nothing left to remove. It drops only the 11,882 cells whose
   raw total was 0, because zero stays zero under rescaling.
4. So in the notebook the count threshold was inert, and only `min_genes=10` did real work.
   Re-running that ordering leaves **251,444** cells.
5. `02b_qc_filtering.py` reverses the order — it filters on raw counts first (line 33) and
   normalizes afterwards (line 50), which is the correct order and makes the threshold mean
   what it says. It leaves **214,721** cells.

The notebook therefore clustered 251,444 cells and the pipeline clusters 214,721, a difference
of 36,723 cells. A different population produces a different partition, so resolution 0.8
lands on 18 clusters in one and 14 in the other.

As a side note, this also explains the old version of this README: it reported 251,444
surviving cells, which is the notebook's number, not `02b`'s.

**Why the fourteen labels that did land are untrustworthy.** Pandas applies a mapping by key,
so cluster `"5"` received whatever label was filed under `"5"`. The labels were applied
consistently and without error. The problem is what those keys *meant*. Key `"5"` was read off
group 5 of the notebook's 18-cluster partition of 251,444 cells. It got stamped onto group 5
of the pipeline's 14-cluster partition of 214,721 cells. Those are different sets of cells
that merely share a number. As stated in the Terms section, a cluster number carries no
meaning outside the exact run that produced it — and that holds across two runs at the *same*
resolution just as much as across two resolutions. So `cell_type_res08` is not reliable.

**What was done instead.** Kevin chose to annotate at resolution 1.1 of the pipeline's own
clustering, which has 18 clusters. That is `04c_labeling_res1_1.py`, and `cell_type_res1_1` is
the annotation behind Figures E and F. `04c` carries its own separately curated 18-entry
mapping — it is **not** this dictionary reused. Of the eighteen keys the two share, only three
agree (`"0"` Erythroid Precursor 1, `"1"` Erythroid Precursor 2, `"9"` Macrophage); the other
fifteen disagree. The `04c` mapping also reuses three label names across pairs of clusters,
which this one never does, and it drops Plasmablast, Lymphocyte/Dendritic Cell/Vessel and
Lipid-Laden Macrophage entirely. The coincidence that both dictionaries have eighteen entries
is exactly that — a coincidence of cluster count, not shared ancestry.

**The withdrawn mapping**, preserved for reference:

| Cluster | Cell type |
|---|---|
| 0 | Erythroid Precursor 1 |
| 1 | Erythroid Precursor 2 |
| 2 | Endothelial |
| 3 | Neoplastic Plasma Cell 1 |
| 4 | Neoplastic Plasma Cell 2 (BCL2L11-high, apoptosis-primed) |
| 5 | Neoplastic Plasma Cell 3 |
| 6 | TREM2+ Macrophage (Tumor-Supportive, Adjacent) |
| 7 | Monocyte |
| 8 | Neoplastic Plasma Cell 4 |
| 9 | Macrophage |
| 10 | Cytotoxic T Cell |
| 11 | Neoplastic Plasma Cell 5 (MYC-high, possibly proliferative) |
| 12 | Plasmablast (B-to-Plasma Cell Transitional) |
| 13 | Stromal Fibroblast |
| 14 | Neoplastic Plasma Cell 6 *(never applied)* |
| 15 | Lymphocyte/Dendritic Cell/Vessel (Intravascular) *(never applied)* |
| 16 | Proliferating Myeloid Precursor *(never applied)* |
| 17 | Lipid-Laden Macrophage (Tumor-Adjacent) *(never applied)* |

**The column is still in the file.** `cell_type_res08` physically exists inside
`combined_clustered.h5ad` today, with 14 categories, alongside its stored palette
`cell_type_res08_colors` in `.uns`. Nothing reads it. It is inert, but it is there, so do not
be surprised to find it when inspecting the object — and do not use it.

---

## Stale outputs currently on disk

Earlier versions of these scripts wrote to different paths. Those older files are kept under
`xenium/commercial/figures/archive/`. They are **not** produced by the pipeline as it stands and
will not be refreshed by re-running it:

- `archive/Fig_D/FigD_umap_by_sample.{png,pdf}`, `archive/Fig_E/FigE_*.{png,pdf}`,
  `archive/Fig_F/FigF_*.{png,pdf}` — the old figure folders. The current scripts write to
  `umap_by_sample/`, `umap_celltypes/` and `spatial_celltypes/` with no figure-letter prefix.
- `archive/old_notebooks/` — plots drawn by the notebooks in `xenium/commercial/notebooks/archive/`
  before this pipeline existed, including `UMAP_celltypes_res06.{png,pdf}` and
  `UMAP_celltypes_res08.{png,pdf}` from a time when the labeling scripts also drew a figure.
  Neither `04a` nor `04c` saves any figure today.

The resolution-sweep mismatch described in the `03_processing.py` section belongs to the same
category: fifteen CSVs on disk, eleven written by the current script.
