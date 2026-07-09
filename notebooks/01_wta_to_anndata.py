
import scanpy as sc
import os

data_path = "/N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/"

results_dir = "results/geomx_exports"
os.makedirs(results_dir, exist_ok=True)

from python.geomx_anndata_utils import (
    build_geomx_wta_anndata,
    split_anndata_by_obs,
    export_expression_csv,
    export_negprobe_background_from_targets,
)

data_path = "/N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/"

adata_wta = build_geomx_wta_anndata(
    data_path=data_path,
    target_counts_file="CAA-AD_Raw_TargetCountMatrix.csv",
    probe_counts_file="CAA-AD_Raw_BioProbeCountMatrix.csv",
    feature_annotations_file="CAA-AD_Feature_Annotations.csv",
    sample_annotations_file="CAA-AD_Sample_Annotations.csv",
    output_file="multiomic-spatial-integration/data/CAA-AD_AnnData.h5ad",
)

adata_wta

adata_wta.X.shape
adata_wta.obsm["negProbes"].shape
adata_wta.obs.head()
adata_wta.var.head()

# Plain CSV export of raw counts (ROI x gene), consumed directly by
# R/pathway_proportion_utils.R for the cell-type-proportion / pathway
# linkage analysis -- avoids re-deriving GeoMx ETL logic in R.
export_expression_csv(
    adata_wta,
    output_file=f"{results_dir}/CAA-AD_expression_wide.csv",
)

# Per-ROI negative-probe background (read from the negative targets in the
# raw TargetCountMatrix), consumed by the SpatialDecon adapter in the
# deconvolution method comparison.
export_negprobe_background_from_targets(
    data_path=data_path,
    target_counts_file="CAA-AD_Raw_TargetCountMatrix.csv",
    feature_annotations_file="CAA-AD_Feature_Annotations.csv",
    output_file=f"{results_dir}/CAA-AD_negprobe_background.csv",
)

adata_by_disease = split_anndata_by_obs(
    adata_wta,
    obs_col="disease_status",
    values=["AD-CAA", "Control"]
)

adata_wta_ADCAA = adata_by_disease["AD-CAA"]
adata_wta_Control = adata_by_disease["Control"]

