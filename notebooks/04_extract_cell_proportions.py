# 04 — Extract and Annotate Spatial Cell-Type Proportions

#This notebook extracts posterior cell-type abundance estimates from trained SpaceJam/Pyro models and converts them into analysis-ready tables.

#Outputs include:

#- absolute cell-type abundance per ROI
#- relative cell-type proportions per ROI
#- long-format ROI × cell-type tables
#- metadata-annotated abundance tables

import scanpy as sc
import pandas as pd
import torch

from python.abundance_extraction_utils import (
    load_variational_means,
    extract_spot_factors,
    normalize_spot_factors,
    make_abundance_dataframe,
    make_long_abundance_table,
    save_abundance_outputs
)
## Load GeoMx AnnData metadata

adata_wta = sc.read_h5ad("data/CAA-AD_AnnData.h5ad")
adata_wta.obs.head()

## Load posterior means from trained spatial models
## Extract spot factors
#`spot_factors` represent inferred absolute cell-type abundance per ROI.
## Normalize to relative proportions
spacejam_results = "/N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/SpaceJam_results"

ad_spot_abs = torch.load(f"{spacejam_results}/ADCAA_spot_factors_abs.pt", map_location="cpu").detach().cpu().numpy()
ad_spot_rel = torch.load(f"{spacejam_results}/ADCAA_spot_factors_rel.pt", map_location="cpu").detach().cpu().numpy()

ctrl_spot_abs = torch.load(f"{spacejam_results}/CTRL_spot_factors_abs.pt", map_location="cpu").detach().cpu().numpy()
ctrl_spot_rel = torch.load(f"{spacejam_results}/CTRL_spot_factors_rel.pt", map_location="cpu").detach().cpu().numpy()

## Add cell-type labels

#These labels should match the regression-derived signature columns.
signature_path = "/N/u/echimal/Quartz/Desktop/CLR_MRI/Human_GeoMx_Sep2025/Regression-model/AD+CAA_inferred_signatures.csv"

signature_df = pd.read_csv(signature_path, index_col=0)
celltype_labels = pd.DataFrame({"celltype": signature_df.columns})
## Build abundance data frames

ad_obs = adata_wta.obs[adata_wta.obs["disease_status"] == "AD-CAA"].copy()
ctrl_obs = adata_wta.obs[adata_wta.obs["disease_status"] == "Control"].copy()

ad_abs_df = make_abundance_dataframe(
    spot_factors=ad_spot_abs,
    obs=ad_obs,
    celltype_labels=celltype_labels["celltype"].tolist()
)

ad_rel_df = make_abundance_dataframe(
    spot_factors=ad_spot_rel,
    obs=ad_obs,
    celltype_labels=celltype_labels["celltype"].tolist()
)

ctrl_abs_df = make_abundance_dataframe(
    spot_factors=ctrl_spot_abs,
    obs=ctrl_obs,
    celltype_labels=celltype_labels["celltype"].tolist()
)

ctrl_rel_df = make_abundance_dataframe(
    spot_factors=ctrl_spot_rel,
    obs=ctrl_obs,
    celltype_labels=celltype_labels["celltype"].tolist()
)

## Create long-format table for statistics

ad_long = make_long_abundance_table(
    abundance_df=ad_rel_df,
    celltype_labels=celltype_labels["celltype"].tolist(),
    abundance_col="rel_abundance"
)

ctrl_long = make_long_abundance_table(
    abundance_df=ctrl_rel_df,
    celltype_labels=celltype_labels["celltype"].tolist(),
    abundance_col="rel_abundance"
)

abundance_long = pd.concat([ad_long, ctrl_long], axis=0)

abundance_long.head()

## Save outputs

save_abundance_outputs(
    output_dir="results/cell_proportions",
    ad_abs_df=ad_abs_df,
    ad_rel_df=ad_rel_df,
    ctrl_abs_df=ctrl_abs_df,
    ctrl_rel_df=ctrl_rel_df,
    long_df=abundance_long
)




