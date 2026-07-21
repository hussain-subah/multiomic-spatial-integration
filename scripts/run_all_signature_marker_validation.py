#!/usr/bin/env python3

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr, spearmanr


PROJECT_DIR = Path(
    "/N/u/echimal/Quartz/Desktop/CLR_MRI/"
    "Human_GeoMx_Sep2025/multiomic-spatial-integration"
)

SPATIAL_FILE = PROJECT_DIR / "data" / "CAA-AD_AnnData.h5ad"

PROPORTION_FILE = (
    PROJECT_DIR
    / "results"
    / "cell_proportions"
    / "spatial_celltype_proportions_for_R.csv"
)

MARKER_FILE = (
    PROJECT_DIR
    / "results"
    / "all_marker_validation"
    / "top45_reference_markers.csv"
)

SIGNATURE_FILE = (
    PROJECT_DIR
    / "results"
    / "regression_model"
    / "AD+CAA_inferred_signatures.csv"
)

OUTPUT_DIR = (
    PROJECT_DIR
    / "results"
    / "all_marker_validation"
)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def safe_correlation(x, y, method):
    """Return correlation and p-value while handling constant vectors."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    keep = np.isfinite(x) & np.isfinite(y)
    x = x[keep]
    y = y[keep]

    if len(x) < 3:
        return np.nan, np.nan, len(x)

    if np.nanstd(x) == 0 or np.nanstd(y) == 0:
        return np.nan, np.nan, len(x)

    if method == "spearman":
        statistic, p_value = spearmanr(x, y)
    elif method == "pearson":
        statistic, p_value = pearsonr(x, y)
    else:
        raise ValueError(f"Unsupported correlation method: {method}")

    return statistic, p_value, len(x)


print("Loading spatial AnnData...")

adata = sc.read_h5ad(SPATIAL_FILE)

if "counts" not in adata.layers:
    raise KeyError("Spatial AnnData does not contain layers['counts'].")

print(adata)


print("\nLoading inferred proportions...")

props = pd.read_csv(PROPORTION_FILE)

required_prop_columns = {
    "ROI_ID",
    "celltype",
    "rel_abundance",
    "disease_status",
    "region",
    "Scan_ID",
}

missing_prop_columns = required_prop_columns - set(props.columns)

if missing_prop_columns:
    raise KeyError(
        "Proportion table is missing columns: "
        f"{sorted(missing_prop_columns)}"
    )

print("Proportion table shape:", props.shape)
print("Unique ROIs:", props["ROI_ID"].nunique())
print("Unique cell types:", props["celltype"].nunique())


print("\nLoading reference markers...")

markers = pd.read_csv(MARKER_FILE)

if "cluster" not in markers.columns:
    raise KeyError(
        "Marker table does not contain a 'cluster' column."
    )

gene_column_candidates = [
    "gene",
    "features",
    "feature",
]

gene_column = next(
    (
        column
        for column in gene_column_candidates
        if column in markers.columns
    ),
    None,
)

if gene_column is None:
    first_column = markers.columns[0]

    if first_column not in {
        "p_val",
        "avg_log2FC",
        "pct.1",
        "pct.2",
        "p_val_adj",
        "cluster",
    }:
        gene_column = first_column
    else:
        raise KeyError(
            "Could not identify the marker gene-name column. "
            f"Columns: {markers.columns.tolist()}"
        )

markers[gene_column] = markers[gene_column].astype(str)
markers["cluster"] = markers["cluster"].astype(str)

print("Marker gene column:", gene_column)
print("Cell types with markers:", markers["cluster"].nunique())


print("\nLoading regression signatures...")

signatures = pd.read_csv(
    SIGNATURE_FILE,
    index_col=0,
)

signatures.index = signatures.index.astype(str)

print("Signature shape:", signatures.shape)


print("\nPreparing normalized GeoMx expression...")

spatial_genes = adata.var_names.astype(str)

counts = adata.layers["counts"]

if not isinstance(counts, np.ndarray):
    counts = counts.toarray()

counts = np.asarray(counts, dtype=np.float64)

library_size = counts.sum(axis=1, keepdims=True)

cpm = (
    counts
    / np.maximum(library_size, 1.0)
) * 1_000_000.0

log_cpm = np.log1p(cpm)

expression_df = pd.DataFrame(
    log_cpm,
    index=adata.obs_names.astype(str),
    columns=spatial_genes,
)


print("\nValidating abundance against independent marker scores...")

results = []
roi_scores = []

celltypes = sorted(
    set(props["celltype"].astype(str))
    & set(markers["cluster"])
)

for celltype in celltypes:

    marker_genes = (
        markers.loc[
            markers["cluster"] == celltype,
            gene_column,
        ]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    present_genes = [
        gene
        for gene in marker_genes
        if gene in expression_df.columns
    ]

    abundance = (
        props.loc[
            props["celltype"].astype(str) == celltype,
            [
                "ROI_ID",
                "rel_abundance",
                "disease_status",
                "region",
                "Scan_ID",
            ],
        ]
        .drop_duplicates(subset="ROI_ID")
        .copy()
    )

    abundance["ROI_ID"] = abundance["ROI_ID"].astype(str)

    if len(present_genes) == 0:

        results.append(
            {
                "celltype": celltype,
                "n_markers_reference": len(marker_genes),
                "n_markers_spatial": 0,
                "spearman_rho": np.nan,
                "spearman_p": np.nan,
                "pearson_r": np.nan,
                "pearson_p": np.nan,
                "n_rois": 0,
                "marker_score_sd": np.nan,
                "mean_rel_abundance": abundance[
                    "rel_abundance"
                ].mean(),
            }
        )

        continue

    marker_score = (
        expression_df[present_genes]
        .mean(axis=1)
        .rename("marker_score")
        .reset_index()
        .rename(columns={"index": "ROI_ID"})
    )

    merged = abundance.merge(
        marker_score,
        on="ROI_ID",
        how="inner",
    )

    rho, p_rho, n_rho = safe_correlation(
        merged["rel_abundance"],
        merged["marker_score"],
        method="spearman",
    )

    pearson_r, pearson_p, _ = safe_correlation(
        merged["rel_abundance"],
        merged["marker_score"],
        method="pearson",
    )

    results.append(
        {
            "celltype": celltype,
            "n_markers_reference": len(marker_genes),
            "n_markers_spatial": len(present_genes),
            "spearman_rho": rho,
            "spearman_p": p_rho,
            "pearson_r": pearson_r,
            "pearson_p": pearson_p,
            "n_rois": n_rho,
            "marker_score_sd": merged["marker_score"].std(),
            "mean_rel_abundance": merged[
                "rel_abundance"
            ].mean(),
        }
    )

    score_output = merged.copy()
    score_output.insert(0, "validated_celltype", celltype)
    roi_scores.append(score_output)

results_df = pd.DataFrame(results)

results_df["spearman_p_adj"] = (
    results_df["spearman_p"]
    .rank(method="first")
)

valid_p = results_df["spearman_p"].notna()

if valid_p.any():
    from statsmodels.stats.multitest import multipletests

    results_df.loc[
        valid_p,
        "spearman_p_adj",
    ] = multipletests(
        results_df.loc[
            valid_p,
            "spearman_p",
        ],
        method="fdr_bh",
    )[1]

results_df["validation_status"] = np.select(
    [
        results_df["n_markers_spatial"] < 3,
        results_df["spearman_rho"].isna(),
        (
            results_df["spearman_rho"] > 0
        )
        & (
            results_df["spearman_p_adj"] < 0.05
        ),
        results_df["spearman_rho"] > 0,
    ],
    [
        "insufficient_markers",
        "unevaluable",
        "positive_significant",
        "positive_nonsignificant",
    ],
    default="nonpositive",
)

results_df = results_df.sort_values(
    "spearman_rho",
    ascending=False,
    na_position="last",
)

results_df.to_csv(
    OUTPUT_DIR / "all_celltype_marker_validation.csv",
    index=False,
)

if roi_scores:
    roi_scores_df = pd.concat(
        roi_scores,
        ignore_index=True,
    )

    roi_scores_df.to_csv(
        OUTPUT_DIR / "all_celltype_roi_marker_scores.csv",
        index=False,
    )


print("\nChecking marker identity within regression signatures...")

signature_gene_mean = signatures.mean(axis=1)
signature_gene_sd = signatures.std(axis=1).replace(0, np.nan)

signature_z = (
    signatures
    .sub(signature_gene_mean, axis=0)
    .div(signature_gene_sd, axis=0)
)

identity_results = []

for celltype in celltypes:

    marker_genes = (
        markers.loc[
            markers["cluster"] == celltype,
            gene_column,
        ]
        .dropna()
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    present_signature_genes = [
        gene
        for gene in marker_genes
        if gene in signature_z.index
    ]

    if len(present_signature_genes) == 0:
        continue

    scores = (
        signature_z.loc[
            present_signature_genes
        ]
        .mean(axis=0)
        .sort_values(ascending=False)
    )

    best_signature = str(scores.index[0])
    self_score = (
        float(scores[celltype])
        if celltype in scores.index
        else np.nan
    )

    self_rank = (
        int(
            scores.index
            .tolist()
            .index(celltype)
        ) + 1
        if celltype in scores.index
        else np.nan
    )

    identity_results.append(
        {
            "celltype": celltype,
            "n_markers_signature": len(
                present_signature_genes
            ),
            "best_signature": best_signature,
            "best_score": float(scores.iloc[0]),
            "self_score": self_score,
            "self_rank": self_rank,
            "self_is_best": best_signature == celltype,
        }
    )

identity_df = pd.DataFrame(identity_results)

identity_df.to_csv(
    OUTPUT_DIR / "regression_signature_marker_identity.csv",
    index=False,
)

# ============================================================
# Combine spatial marker concordance with signature identity
# ============================================================

combined = results_df.merge(
    identity_df,
    on="celltype",
    how="left",
)

# ============================================================
# Nuanced QC classification
# ============================================================

combined["qc_category"] = np.select(
    [
        # Too few spatially detected markers to evaluate reliably
        combined["n_markers_spatial"] < 3,

        # Correlation could not be calculated
        combined["spearman_rho"].isna(),

        # Abundance does not positively track its own marker score
        combined["spearman_rho"] <= 0,

        # Positive relationship, but it does not survive FDR correction
        (
            combined["spearman_rho"] > 0
        )
        & (
            combined["spearman_p_adj"] >= 0.05
        ),

        # Spatial signal is supported, but subtype identity is ambiguous
        (
            combined["spearman_rho"] > 0
        )
        & (
            combined["spearman_p_adj"] < 0.05
        )
        & (
            combined["self_rank"] > 3
        ),

        # Spatial marker concordance is significant and the intended
        # signature ranks within the top three matches
        (
            combined["spearman_rho"] > 0
        )
        & (
            combined["spearman_p_adj"] < 0.05
        )
        & (
            combined["self_rank"] <= 3
        ),
    ],
    [
        "insufficient_markers",
        "unevaluable",
        "marker_concordance_failure",
        "weak_marker_concordance",
        "lineage_valid_subtype_ambiguous",
        "validated",
    ],
    default="unevaluable",
)

# Only the most serious categories are treated as potential problems.
# Subtype ambiguity is retained as a warning, not a model failure.
combined["potential_problem"] = combined["qc_category"].isin(
    [
        "insufficient_markers",
        "unevaluable",
        "marker_concordance_failure",
    ]
)

# ============================================================
# Human-readable QC interpretation
# ============================================================

combined["qc_interpretation"] = np.select(
    [
        combined["qc_category"] == "validated",

        combined["qc_category"]
        == "lineage_valid_subtype_ambiguous",

        combined["qc_category"]
        == "weak_marker_concordance",

        combined["qc_category"]
        == "marker_concordance_failure",

        combined["qc_category"]
        == "insufficient_markers",

        combined["qc_category"]
        == "unevaluable",
    ],
    [
        (
            "Positive, FDR-significant marker concordance and "
            "acceptable regression-signature identity."
        ),
        (
            "Positive, FDR-significant spatial marker concordance, "
            "but subtype markers overlap more strongly with another "
            "related regression factor."
        ),
        (
            "Positive marker concordance, but the association does "
            "not survive FDR correction."
        ),
        (
            "Inferred abundance is non-positive or negatively "
            "associated with its independent spatial marker score."
        ),
        (
            "Too few reference markers were detected in the GeoMx "
            "dataset for reliable evaluation."
        ),
        (
            "Marker-abundance correlation could not be evaluated."
        ),
    ],
    default="Not classified.",
)

# ============================================================
# Add subtype identity details
# ============================================================

combined["self_is_top3"] = (
    combined["self_rank"].notna()
    & (
        combined["self_rank"] <= 3
    )
)

combined["identity_status"] = np.select(
    [
        combined["self_rank"].isna(),

        combined["self_rank"] == 1,

        (
            combined["self_rank"] > 1
        )
        & (
            combined["self_rank"] <= 3
        ),

        combined["self_rank"] > 3,
    ],
    [
        "identity_unevaluable",
        "self_best_match",
        "self_within_top3",
        "subtype_identity_ambiguous",
    ],
    default="identity_unevaluable",
)

# ============================================================
# Sort by QC severity
# ============================================================

qc_order = {
    "marker_concordance_failure": 1,
    "insufficient_markers": 2,
    "unevaluable": 3,
    "weak_marker_concordance": 4,
    "lineage_valid_subtype_ambiguous": 5,
    "validated": 6,
}

combined["qc_order"] = combined["qc_category"].map(
    qc_order
)

combined = combined.sort_values(
    [
        "qc_order",
        "spearman_rho",
        "self_rank",
    ],
    ascending=[
        True,
        True,
        False,
    ],
    na_position="last",
)

# ============================================================
# Save complete validation summary
# ============================================================

combined.to_csv(
    OUTPUT_DIR / "combined_factor_validation_summary.csv",
    index=False,
)

# Save smaller category-specific tables for easier review
combined.loc[
    combined["potential_problem"]
].to_csv(
    OUTPUT_DIR / "factors_requiring_attention.csv",
    index=False,
)

combined.loc[
    combined["qc_category"]
    == "lineage_valid_subtype_ambiguous"
].to_csv(
    OUTPUT_DIR / "spatially_valid_subtype_ambiguous_factors.csv",
    index=False,
)

combined.loc[
    combined["qc_category"] == "validated"
].to_csv(
    OUTPUT_DIR / "validated_factors.csv",
    index=False,
)

# ============================================================
# Console summaries
# ============================================================

print("\nValidation completed.")

print("\nSpatial marker-concordance status counts:")
print(
    results_df["validation_status"]
    .value_counts(dropna=False)
    .to_string()
)

print("\nNuanced QC category counts:")
print(
    combined["qc_category"]
    .value_counts(dropna=False)
    .to_string()
)

print("\nSignature identity status counts:")
print(
    combined["identity_status"]
    .value_counts(dropna=False)
    .to_string()
)

print("\nLowest marker-abundance correlations:")
print(
    combined[
        [
            "celltype",
            "n_markers_spatial",
            "spearman_rho",
            "spearman_p_adj",
            "self_rank",
            "best_signature",
            "mean_rel_abundance",
            "qc_category",
        ]
    ]
    .sort_values(
        "spearman_rho",
        ascending=True,
        na_position="last",
    )
    .head(15)
    .to_string(index=False)
)

print("\nFactors requiring attention:")
attention_columns = [
    "celltype",
    "n_markers_reference",
    "n_markers_spatial",
    "spearman_rho",
    "spearman_p_adj",
    "self_rank",
    "best_signature",
    "mean_rel_abundance",
    "qc_category",
    "qc_interpretation",
]

attention_df = combined.loc[
    combined["potential_problem"],
    attention_columns,
]

if attention_df.empty:
    print("None.")
else:
    print(
        attention_df.to_string(index=False)
    )

print(
    "\nSpatially supported factors with ambiguous "
    "subtype identity:"
)

ambiguous_columns = [
    "celltype",
    "n_markers_spatial",
    "spearman_rho",
    "spearman_p_adj",
    "self_rank",
    "best_signature",
    "best_score",
    "self_score",
    "mean_rel_abundance",
]

ambiguous_df = combined.loc[
    combined["qc_category"]
    == "lineage_valid_subtype_ambiguous",
    ambiguous_columns,
]

if ambiguous_df.empty:
    print("None.")
else:
    print(
        ambiguous_df.to_string(index=False)
    )

print("\nValidated factors:")
validated_columns = [
    "celltype",
    "n_markers_spatial",
    "spearman_rho",
    "spearman_p_adj",
    "self_rank",
    "best_signature",
    "mean_rel_abundance",
]

validated_df = combined.loc[
    combined["qc_category"] == "validated",
    validated_columns,
].sort_values(
    "spearman_rho",
    ascending=False,
)

if validated_df.empty:
    print("None.")
else:
    print(
        validated_df.to_string(index=False)
    )

print("\nOutputs saved under:")
print(OUTPUT_DIR)

print("\nPrimary outputs:")
print(
    OUTPUT_DIR
    / "combined_factor_validation_summary.csv"
)
print(
    OUTPUT_DIR
    / "factors_requiring_attention.csv"
)
print(
    OUTPUT_DIR
    / "spatially_valid_subtype_ambiguous_factors.csv"
)
print(
    OUTPUT_DIR
    / "validated_factors.csv"
)

