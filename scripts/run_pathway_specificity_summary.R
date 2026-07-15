# ============================================================
# Pathway-proportion linkage: specificity summary
#
# The per-celltype GSEA results in results/pathway_proportion_link/ include
# a large fraction of generic translation/ribosome/OXPHOS/olfactory-receptor
# pathways that recur as "significant" across most cell types and strata --
# the signature of a shared technical or compositional confound rather than
# cell-type-specific biology (proportions sum to ~1, so a real signal in one
# dominant cell type mechanically anti-correlates with many others).
#
# This script counts, per pathway, how many (celltype, stratum) combinations
# it's significant in, flags generic/housekeeping pathways, and surfaces the
# low-frequency, non-housekeeping hits per stratum as the more credible,
# cell-type-restricted signal.
# ============================================================

setwd("D:/repo/multiomic-spatial-integration-branchv2")

source("R/viz_theme.R")
source("R/pathway_proportion_utils.R")
source("R/viz_pathway.R")

pathway_dir <- "results/pathway_proportion_link"
fig_dir <- "results/figures/pathway_specificity"

PADJ_CUTOFF <- 0.05
MAX_FREQUENCY <- 10   # keep pathways significant in at most this many of the 414 (celltype, stratum) combinations
TOP_N_PER_STRATUM <- 15

combined <- load_pathway_enrichment_dir(pathway_dir)
message(sprintf(
  "Loaded %d pathway x celltype x stratum rows from %d files.",
  nrow(combined), length(unique(paste(combined$stratum, combined$celltype)))
))

specificity <- compute_pathway_specificity(combined, padj_cutoff = PADJ_CUTOFF)

message(sprintf(
  "%d / %d (%.1f%%) significant hits are flagged as generic/housekeeping.",
  sum(specificity$significant$is_housekeeping), nrow(specificity$significant),
  100 * mean(specificity$significant$is_housekeeping)
))

top_specific <- top_specific_pathways_by_stratum(
  specificity, max_frequency = MAX_FREQUENCY, top_n = TOP_N_PER_STRATUM
)

message(sprintf(
  "%d cell-type-restricted, non-housekeeping hits kept across %d strata (max_frequency <= %d).",
  nrow(top_specific), length(unique(top_specific$stratum)), MAX_FREQUENCY
))

utils::write.csv(
  specificity$frequency,
  file.path(pathway_dir, "pathway_specificity_summary.csv"),
  row.names = FALSE
)

utils::write.csv(
  top_specific,
  file.path(pathway_dir, "top_specific_pathways_by_stratum.csv"),
  row.names = FALSE
)

message("\nTop specific hit per stratum:")
for (stratum in unique(top_specific$stratum)) {
  row <- top_specific[top_specific$stratum == stratum, ][1, ]
  message(sprintf(
    "  %-28s %-20s %-60s NES=%.2f padj=%.2e (global n_sig=%d)",
    stratum, row$celltype, substr(row$pathway, 1, 58),
    row$NES, row$padj, row$n_significant
  ))
}

plot_specificity_dotplot(top_specific, fig_dir)

message("\nWrote ", file.path(pathway_dir, "pathway_specificity_summary.csv"))
message("Wrote ", file.path(pathway_dir, "top_specific_pathways_by_stratum.csv"))
message("Pathway specificity summary completed.")
