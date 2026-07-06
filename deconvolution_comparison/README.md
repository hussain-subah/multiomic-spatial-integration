# Deconvolution method comparison
Benchmarks alternative cell-type deconvolution tools against the pipeline's
existing Cell2location / SpaceJam result, on the same GeoMx WTA ROI data.

The GeoMx platform used here is **ROI/segment-level (bulk-like), with no
per-spot spatial coordinates** anywhere in the pipeline. That constraint
determines which tools are applicable (see the feasibility table below) and
means the comparison is **agreement-based, not accuracy-based**: there is no
ground-truth composition for the real ROIs, so we measure how well each
method's inferred proportions concord with the Cell2location baseline and
with each other, not which method is "correct." (For a known-ground-truth
accuracy check on synthetic data, see the separate pseudobulk validation:
`python/pseudobulk_utils.py` / `notebooks/06_pseudobulk_validation.py`.)

## Feasibility

| Tool | Type | Status here | Notes |
|------|------|-------------|-------|
| **Cell2location** | Bayesian, scRNA-seq ref | **Baseline** | Existing pipeline result; comparison target, no new adapter. |
| **SpatialDecon** | Log-normal regression, GeoMx-native | Applicable (best fit) | NanoString's own DSP deconvolution — the most platform-appropriate method here. Uses negative-probe background. `adapter_spatialdecon.R` |
| **MuSiC** | Bulk, scRNA-seq ref | Applicable | Purpose-built for bulk; uses cross-subject variance. `adapter_music.R` |
| **Bisque** | Bulk, scRNA-seq ref | Applicable | Reference-based decomposition for bulk. `adapter_bisque.R` |
| **DWLS** | Bulk, signature matrix | Applicable | Uses the signature matrix (independent raw-reference signature by default). `adapter_dwls.R` |
| **SPOTlight** | NMF+NNLS, sc ref | Applicable | Runs on ROI mixtures; needs marker genes. `adapter_spotlight.R` |
| **RCTD / spacexr** | Spatial spot ref | Applicable* | *Runs with **dummy coordinates** — spatial benefit is lost; acts as a robust bulk-ish reference method. `adapter_rctd.R` |
| **BayesPrism** | Bayesian, scRNA-seq ref | Applicable | Works on bulk-like mixtures; heavier to run. `adapter_bayesprism.R` |
| **STdeconvolve** | Reference-free (LDA) | Applicable (indirect) | Produces topics, not cell types; topics are matched to reference signatures post-hoc, so comparison is approximate. `adapter_stdeconvolve.R` |
| **CIBERSORTx** | Bulk, signature matrix | **Input-prep only** | Requires a license token + web/Docker submission; cannot run as a plain function call. `scripts/prepare_cibersortx_inputs.R` prepares the input files. |
| **CARD** | Spatially-informed | **Excluded** | Requires per-spot (x,y) coordinates, which this data does not have. `CARDfree` would discard its spatial premise, so it is not included. |

## Inputs (set paths in each script's config header)

- **Mixture (ROI x gene)**: `data/CAA-AD_expression_wide.csv`
  (exported by `python.geomx_anndata_utils.export_expression_csv`, called in
  `notebooks/01_wta_to_anndata.py`).
- **snRNA-seq reference (counts + labels)**: the MatrixMarket triple exported
  by `scripts/run_sn_reference_export.R`
  (`AD_CAA_counts.mtx`, `AD_CAA_metadata.csv`, `AD_CAA_genes.csv`), or the
  full Seurat RDS checkpoint.
- **Signature matrix (gene x celltype)** for the signature-based tools (DWLS,
  SpatialDecon, CIBERSORTx). Two options, set via `signature_source` in
  `run_all_deconvolution.R`:
  - **`"independent"` (default)**: a per-cell-type mean aggregated from the raw
    snRNA-seq reference counts (`build_mean_signature()`), written to
    `data/independent_reference_signature.csv`. Not produced by cell2location,
    so it does **not** share a basis with the Cell2location/SpaceJam baseline —
    this is the fair comparison.
  - **`"inferred"`**: the cell2location NB-regression signatures
    (`AD+CAA_inferred_signatures.csv`, from `notebooks/02_regression_signatures.py`).
    Shares a basis with the baseline, so DWLS/SpatialDecon agreement with it is
    partly built-in.
  Note a signature can only come from the labeled reference, never from the
  GeoMx target count matrix (that is the unlabeled mixture being deconvolved).
- **Per-ROI negative-probe background**: `data/CAA-AD_negprobe_background.csv`
  (exported by `python.geomx_anndata_utils.export_negprobe_background_from_targets`,
  called in `notebooks/01_wta_to_anndata.py`) — the mean of each ROI's negative
  control targets, read from the negative rows of the raw TargetCountMatrix
  (where the GeoMx negatives live). Used by SpatialDecon; if absent, the adapter
  falls back to a proxy background with a warning.
- **Cell2location baseline (comparison target)**:
  `results/cell_proportions/spatial_celltype_proportions_for_R.csv`.

## Output contract

Every adapter returns, and every run script writes, the same standard
long-format table so the comparison step is method-agnostic:

```
results/deconvolution_comparison/<method>_proportions.csv
    columns: method, ROI_ID, celltype, proportion   (proportions sum to ~1 per ROI)
```

## How to run

1. Ensure the inputs above exist (run notebooks 01-02 and the reference export first).
2. `Rscript deconvolution_comparison/scripts/run_all_deconvolution.R`
   — runs every R adapter whose package is installed; skips (with a message)
   any whose package is missing, so one absent dependency doesn't block the rest.
3. `Rscript deconvolution_comparison/scripts/prepare_cibersortx_inputs.R`
   — writes CIBERSORTx-format signature + mixture TSVs; submit those separately,
   then drop the returned proportions into the standard output location.
4. `Rscript deconvolution_comparison/scripts/run_deconvolution_comparison.R`
   — loads all available `*_proportions.csv` plus the Cell2location baseline and
   writes concordance summaries.

## Status

**None of this has been executed.**
