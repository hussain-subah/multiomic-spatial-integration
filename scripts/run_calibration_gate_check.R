# ============================================================
# Calibration validation gate
#
# AD-CAA and Control were deconvolved by two independently trained
# SpaceJam/Cell2location models with no shared calibration anchor. Before
# any cross-condition absolute-abundance contrast (Disease/Overall/
# MaxPathology) is trusted, this checks whether the per-ROI calibration
# ratio (total abs_abundance / total_counts) shifts by disease_status more
# than plausible within-condition noise -- a proxy for an arbitrary,
# condition-specific scale on top of any real biology.
#
# The Amyloid-effect contrast (entirely within AD-CAA, a single posterior)
# is gauge-safe by construction and is always "proceed" regardless of this
# check.
#
# Must be run before scripts/run_absolute_abundance_contrasts.R, which
# looks for this script's decision file.
# ============================================================

library(glmmTMB)
library(dplyr)
library(tidyr)
library(readr)

source("R/absolute_calibration_utils.R")

setwd("D:/repo/multiomic-spatial-integration-branchv2")

abs_long_file <- "results/cell_proportions/roi_celltype_abundance_long_abs.csv"
total_counts_file <- "results/cell_proportions/roi_total_counts.csv"
output_dir <- "results/spatial_stats_absolute"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

abs_df <- read.csv(abs_long_file, stringsAsFactors = FALSE)
total_counts_df <- read.csv(total_counts_file, stringsAsFactors = FALSE)

# ============================================================
# Calibration ratio table + gate decision
# ============================================================

calib_df <- compute_calibration_ratio_table(
  abs_df = abs_df,
  total_counts_df = total_counts_df
)

write.csv(
  calib_df,
  file.path(output_dir, "calibration_ratio_table.csv"),
  row.names = FALSE
)

gate_result <- run_calibration_gate_check(
  calib_df = calib_df,
  biological_threshold_log2 = 0.5
)

write.csv(
  gate_result$within_condition_cv,
  file.path(output_dir, "calibration_gate_within_condition_cv.csv"),
  row.names = FALSE
)

# Machine-readable summary (used by plot_calibration_gate_diagnostic() so it
# doesn't need to re-parse the human-readable decision .txt or refit the model)
gate_summary_df <- data.frame(
  decision = gate_result$decision,
  disease_coefficient_log2 = gate_result$disease_coefficient_log2,
  disease_se_log2 = gate_result$disease_se_log2,
  biological_threshold_log2 = gate_result$biological_threshold_log2,
  region_pattern_correlation = gate_result$region_pattern_correlation
)

write.csv(
  gate_summary_df,
  file.path(output_dir, "calibration_gate_summary.csv"),
  row.names = FALSE
)

# ============================================================
# Per-contrast decision: Amyloid effect is exempt (single posterior,
# gauge-safe by construction); Disease/Overall/MaxPathology all share the
# same cross-condition gate decision, since all three compare AD-CAA vs
# Control ROIs drawn from the two independently trained models.
# ============================================================

decision_lines <- c(
  paste0("Amyloid: proceed  # exempt, entirely within AD-CAA, single posterior"),
  paste0("Disease: ", gate_result$decision),
  paste0("Overall: ", gate_result$decision),
  paste0("MaxPathology: ", gate_result$decision),
  "",
  paste0("disease_coefficient_log2: ", round(gate_result$disease_coefficient_log2, 4)),
  paste0("disease_se_log2: ", round(gate_result$disease_se_log2, 4)),
  paste0("biological_threshold_log2: ", gate_result$biological_threshold_log2),
  paste0("region_pattern_correlation: ", round(gate_result$region_pattern_correlation, 4))
)

writeLines(decision_lines, file.path(output_dir, "calibration_gate_decision.txt"))

message(
  "Calibration gate decision (cross-condition contrasts): ", gate_result$decision, "\n",
  "disease_status coefficient: ", round(gate_result$disease_coefficient_log2, 3),
  " log2 (threshold ", gate_result$biological_threshold_log2, ")\n",
  "Written to ", file.path(output_dir, "calibration_gate_decision.txt")
)

if (gate_result$decision == "do_not_proceed") {
  message(
    "\ndo_not_proceed: only the Amyloid effect should be treated as reliable ",
    "from the absolute-abundance pipeline. Disease/Overall/MaxPathology are ",
    "exploratory-only -- scripts/run_absolute_abundance_contrasts.R will still ",
    "run them (stamped gate_status = 'do_not_proceed'), but they should not ",
    "be reported as confirmatory."
  )
}

message("Calibration gate check completed successfully.")
