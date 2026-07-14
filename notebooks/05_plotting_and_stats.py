# =============================================================================
# 05 — Exploratory Plotting and Statistical Handoff
#
# This notebook performs exploratory visualization of spatial cell-type
# abundance estimates derived from the SpaceJam / Cell2Location workflow.
#
# Goals:
#   - confirm ROI metadata structure
#   - visualize inferred cell-type proportions
#   - inspect disease-, pathology-, and region-associated trends
#   - generate quick non-parametric screening results
#   - validate the table used for formal R mixed-effects modeling
#
# Formal inference is performed in R using beta mixed-effects models.
# =============================================================================

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from python.plotting_utils import (
    load_abundance_tables,
    summarize_roi_counts,
    plot_celltype_by_group,
    plot_celltype_by_pathology_and_disease,
    plot_celltype_by_region_and_disease,
    plot_mean_composition_stacked,
    run_mannwhitney_screen,
)


# =============================================================================
# Configuration
# =============================================================================

PROJECT_DIR = Path(
    "/N/u/echimal/Quartz/Desktop/CLR_MRI/"
    "Human_GeoMx_Sep2025/multiomic-spatial-integration"
)

CELL_PROPORTION_DIR = (
    PROJECT_DIR / "results" / "cell_proportions"
)

FIGURE_DIR = (
    PROJECT_DIR / "figures" / "exploratory"
)

FIGURE_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

wide_file = (
    CELL_PROPORTION_DIR /
    "cell2location_abundance_wide_with_meta.csv"
)

long_file = (
    CELL_PROPORTION_DIR /
    "cell2location_abundance_long_rel.csv"
)

r_handoff_file = (
    CELL_PROPORTION_DIR /
    "spatial_celltype_proportions_for_R.csv"
)

stats_output_file = (
    CELL_PROPORTION_DIR /
    "celltype_disease_stats_mannwhitney.csv"
)


# =============================================================================
# Plot settings
# =============================================================================

sns.set(
    style="whitegrid",
    context="talk",
)

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42


# =============================================================================
# Validate required files
# =============================================================================

required_files = [
    wide_file,
    long_file,
]

missing_files = [
    path
    for path in required_files
    if not path.exists()
]

if missing_files:
    raise FileNotFoundError(
        "Required Notebook 4 output files were not found:\n" +
        "\n".join(str(path) for path in missing_files)
    )

print("Wide input:", wide_file)
print("Long input:", long_file)


# =============================================================================
# Load abundance tables
# =============================================================================

df_wide, df_long = load_abundance_tables(
    wide_file=str(wide_file),
    long_file=str(long_file),
)

print("\nWide table shape:", df_wide.shape)
print("Long table shape:", df_long.shape)

print("\nLong-table preview:")
print(df_long.head())


# =============================================================================
# Validate the long-format table
# =============================================================================

required_columns = {
    "ROI_ID",
    "Scan_ID",
    "disease_status",
    "pathology",
    "region",
    "celltype",
    "rel_abundance",
}

missing_columns = (
    required_columns -
    set(df_long.columns)
)

if missing_columns:
    raise KeyError(
        "The long-format abundance table is missing required columns: "
        f"{sorted(missing_columns)}"
    )

df_long["rel_abundance"] = pd.to_numeric(
    df_long["rel_abundance"],
    errors="raise",
)

if df_long["rel_abundance"].isna().any():
    raise ValueError(
        "Missing relative-abundance values were found."
    )

if (df_long["rel_abundance"] < 0).any():
    raise ValueError(
        "Negative relative-abundance values were found."
    )

if (df_long["rel_abundance"] > 1).any():
    raise ValueError(
        "Relative-abundance values greater than one were found."
    )

n_rois = df_long["ROI_ID"].nunique()
n_celltypes = df_long["celltype"].nunique()

expected_rows = n_rois * n_celltypes

print("\nDataset validation:")
print("ROIs:", n_rois)
print("Cell types:", n_celltypes)
print("Rows:", len(df_long))
print("Expected ROI × cell-type rows:", expected_rows)
print(
    "Relative abundance range:",
    df_long["rel_abundance"].min(),
    "to",
    df_long["rel_abundance"].max(),
)

if len(df_long) != expected_rows:
    raise ValueError(
        "The long table is not a complete ROI × cell-type grid: "
        f"{len(df_long)} rows observed versus "
        f"{expected_rows} expected."
    )

roi_sums = (
    df_long
    .groupby("ROI_ID", observed=True)["rel_abundance"]
    .sum()
)

print(
    "ROI composition sum range:",
    roi_sums.min(),
    "to",
    roi_sums.max(),
)

if not (
    (roi_sums - 1.0).abs() < 1e-5
).all():
    raise ValueError(
        "One or more ROI-level compositions do not sum to one."
    )

print("Long-format abundance validation passed.")


# =============================================================================
# ROI count sanity check
# =============================================================================

roi_summary = summarize_roi_counts(
    df_long,
    roi_col="ROI_ID",
    group_col="disease_status",
)

if roi_summary is not None:
    print("\nROI counts by disease status:")
    print(roi_summary)


# =============================================================================
# Plot selected cell types by disease status
# =============================================================================

available_celltypes = set(
    df_long["celltype"].astype(str)
)

disease_plot_celltypes = [
    "Microglia1",
    "Astrocytes1",
    "Endothelial",
]

for celltype in disease_plot_celltypes:

    if celltype not in available_celltypes:
        print(
            f"Skipping {celltype}: not present in the abundance table."
        )
        continue

    plot_celltype_by_group(
        df_long,
        celltype=celltype,
        group="disease_status",
        save_file=str(
            FIGURE_DIR /
            f"{celltype}_by_disease_status.pdf"
        ),
    )


# =============================================================================
# Plot selected cell types by pathology and disease status
# =============================================================================

pathology_plot_celltypes = [
    "Oligodendrocytes1",
    "SMC",
]

for celltype in pathology_plot_celltypes:

    if celltype not in available_celltypes:
        print(
            f"Skipping {celltype}: not present in the abundance table."
        )
        continue

    plot_celltype_by_pathology_and_disease(
        df_long,
        celltype=celltype,
        save_file=str(
            FIGURE_DIR /
            f"{celltype}_by_pathology_disease.pdf"
        ),
    )


# =============================================================================
# Plot selected cell type by region and disease status
# =============================================================================

region_plot_celltype = "Microglia3"

if region_plot_celltype in available_celltypes:

    plot_celltype_by_region_and_disease(
        df_long,
        celltype=region_plot_celltype,
        save_file=str(
            FIGURE_DIR /
            f"{region_plot_celltype}_by_region_disease.pdf"
        ),
    )

else:
    print(
        f"Skipping {region_plot_celltype}: "
        "not present in the abundance table."
    )


# =============================================================================
# Mean composition by disease status
# =============================================================================

mean_composition = plot_mean_composition_stacked(
    df_long,
    group_col="disease_status",
    save_file=str(
        FIGURE_DIR /
        "mean_celltype_composition_by_disease.pdf"
    ),
)

print("\nMean composition:")
print(mean_composition)


# =============================================================================
# Quick exploratory non-parametric screen
#
# This is not the final statistical model. It is used only to identify
# exploratory trends before formal mixed-effects modeling.
# =============================================================================

stats_df = run_mannwhitney_screen(
    df_long,
    group_col="disease_status",
    group_1="AD-CAA",
    group_2="Control",
    abundance_col="rel_abundance",
    min_n=5,
)

stats_df.to_csv(
    stats_output_file,
    index=False,
)

print("\nTop exploratory Mann–Whitney results:")
print(stats_df.head(10))

print(
    "\nSaved exploratory statistics:",
    stats_output_file,
)


# =============================================================================
# Validate the formal R handoff file
#
# Formal R model:
#
# rel_abundance ~ contrast_variable + (1 | Scan_ID)
#
# with:
#
# glmmTMB::beta_family(link = "logit")
# =============================================================================

# Notebook 4 already produced this file. Confirm that it is identical to
# the validated long table before leaving it in place.

if r_handoff_file.exists():

    existing_r_df = pd.read_csv(
        r_handoff_file,
        low_memory=False,
    )

    comparison_columns = [
        "ROI_ID",
        "Scan_ID",
        "disease_status",
        "pathology",
        "region",
        "celltype",
        "rel_abundance",
    ]

    left = (
        df_long[comparison_columns]
        .sort_values(
            ["ROI_ID", "celltype"]
        )
        .reset_index(drop=True)
    )

    right = (
        existing_r_df[comparison_columns]
        .sort_values(
            ["ROI_ID", "celltype"]
        )
        .reset_index(drop=True)
    )

    if left.shape != right.shape:
        raise ValueError(
            "The existing R handoff table has a different shape "
            "from the validated long table."
        )

    metadata_columns = [
        "ROI_ID",
        "Scan_ID",
        "disease_status",
        "pathology",
        "region",
        "celltype",
    ]

    if not left[metadata_columns].equals(
        right[metadata_columns]
    ):
        raise ValueError(
            "The existing R handoff table has different metadata "
            "or ordering from the validated long table."
        )

    abundance_difference = (
        left["rel_abundance"] -
        right["rel_abundance"]
    ).abs().max()

    if abundance_difference > 1e-10:
        raise ValueError(
            "The existing R handoff table differs numerically from "
            "the validated long table. Maximum difference: "
            f"{abundance_difference}"
        )

    print(
        "\nExisting R handoff file matches the validated long table:"
    )
    print(r_handoff_file)

else:
    # Recreate it only if it is missing.
    df_long.to_csv(
        r_handoff_file,
        index=False,
    )

    print(
        "\nR handoff file was missing and has been created:"
    )
    print(r_handoff_file)


print("\nNotebook 5 completed successfully.")
