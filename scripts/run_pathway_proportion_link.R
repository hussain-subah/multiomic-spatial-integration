# ============================================================
# Pathway analysis linked to cell-type proportions
#
# For each cell type, ranks genes by (Spearman) correlation between
# per-ROI gene expression and that cell type's per-ROI relative
# abundance, then optionally runs rank-based gene set enrichment (fgsea)
# on the ranked list. Cell-type proportions are compositional and can't
# be fed into pathway tools directly -- this is the proportion-driven
# analysis that stands in for that.
# ============================================================

library(fgsea)

source("R/pathway_proportion_utils.R")
source("R/enrichment_utils.R")

# ============================================================
# Config -- adjust paths for your environment
# ============================================================

expression_csv <- "results/geomx_exports/CAA-AD_expression_wide.csv"
proportions_csv <- "results/cell_proportions/roi_celltype_abundance_long.csv"
#roi_id_col <- "ROI_ID"
output_dir <- "results/pathway_proportion_link"

# Set to a local GMT path (e.g. MSigDB Hallmark/KEGG/Reactome/GO) to enable
# pathway enrichment; leave NULL to only produce ranked gene lists.
gmt_file <- "resources/gmt/Reactome_2022.gmt" #KEGG_2021_Human.gmt # or #WikiPathway_2023_Human.gmt

# NULL = all cell types found in proportions_csv; or a character vector to
# restrict to specific cell types.
celltypes <- NULL

# NULL = pooled across all ROIs; or e.g. "pathology" / "region" /
# "disease_pathology_group" / "disease_pathology_region_group" (see below)
# to run the analysis separately within each level of that column (avoids
# conflating a cell-type/gene association with a between-condition
# difference). Note: per-stratum n gets small once region is folded in too
# (as low as ~13 ROIs for AD-CAA/Amyloid/Capillaries) -- expect noisier
# rankings/GSEA there than in the coarser disease_pathology_group split.
stratify_by <- "disease_pathology_region_group"

# Optional pre-filters applied to proportions_df before anything else runs,
# for restricting to an intersection of columns that stratify_by alone can't
# express (e.g. pathology == "Amyloid" AND region == "Capillaries").
# NULL = no filtering on that column.
pathology_filter <- NULL
region_filter <- NULL

roi_id_col <- "ROI_ID"
proportion_col <- "rel_abundance"
celltype_col <- "celltype"
method <- "spearman"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# Load data
# ============================================================

expr_raw <- load_roi_expression(expression_csv)
expr_norm <- normalize_expression_cpm(expr_raw)

proportions_df <- read.csv(proportions_csv, stringsAsFactors = FALSE)

# Derived Control / AD-AmyloidFree / AD-Amyloid grouping. There's no
# Control x Amyloid combination in this dataset (Control ROIs are always
# AmyloidFree), so this collapses the theoretical 2x2 of disease_status x
# pathology down to the 3 groups that actually exist.
proportions_df$disease_pathology_group <- ifelse(
  proportions_df$disease_status == "Control",
  "Control",
  paste0("AD_", proportions_df$pathology)
)

# Same 3 groups, further split by vascular region (all 3 x 3 combinations
# exist in this dataset, with per-stratum ROI counts from ~13 to ~33).
proportions_df$disease_pathology_region_group <- paste(
  proportions_df$disease_pathology_group,
  proportions_df$region,
  sep = "_"
)

if (!is.null(pathology_filter)) {
  proportions_df <- proportions_df[proportions_df$pathology == pathology_filter, ]
}
if (!is.null(region_filter)) {
  proportions_df <- proportions_df[proportions_df$region == region_filter, ]
}
if (nrow(proportions_df) == 0) {
  stop("pathology_filter/region_filter left 0 rows in proportions_df -- check the values match the CSV.")
}

if (is.null(celltypes)) {
  celltypes <- unique(proportions_df[[celltype_col]])
}

# ============================================================
# Per cell type (x stratum, if set)
# ============================================================

strata <- if (is.null(stratify_by)) {
  list(all = proportions_df)
} else {
  split(proportions_df, proportions_df[[stratify_by]])
}

skipped <- data.frame(
  celltype = character(0),
  stratum = character(0),
  reason = character(0),
  stringsAsFactors = FALSE
)

for (stratum_name in names(strata)) {

  stratum_df <- strata[[stratum_name]]

  for (ct in celltypes) {

    message(sprintf("Ranking genes for celltype = %s, stratum = %s", ct, stratum_name))

    ranked <- tryCatch(
      rank_genes_by_celltype_association(
        expr_mat = expr_norm,
        proportions_df = stratum_df,
        celltype = ct,
        roi_id_col = roi_id_col,
        proportion_col = proportion_col,
        celltype_col = celltype_col,
        method = method
      ),
      error = function(e) {
        message("  Skipped: ", conditionMessage(e))
        skipped <<- rbind(skipped, data.frame(
          celltype = ct,
          stratum = stratum_name,
          reason = conditionMessage(e),
          stringsAsFactors = FALSE
        ))
        NULL
      }
    )

    if (is.null(ranked)) next

    ct_dir <- file.path(output_dir, stratum_name, ct)
    dir.create(ct_dir, recursive = TRUE, showWarnings = FALSE)

    write.csv(
      ranked,
      file = file.path(ct_dir, "ranked_genes.csv"),
      row.names = FALSE
    )

    sig_genes <- extract_top_genes(
      data.frame(gene = ranked$gene, p_val_adj = ranked$p_adj),
      padj_cutoff = 0.05
    )
    sig_genes <- format_gene_list(sig_genes)

    write.table(
      sig_genes,
      file = file.path(ct_dir, "significant_gene_list.txt"),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )

    enrich_res <- run_pathway_enrichment(ranked, gmt_file = gmt_file)

    if (!is.null(enrich_res)) {
      enrich_res$leadingEdge <- vapply(
        enrich_res$leadingEdge,
        paste,
        collapse = ";",
        FUN.VALUE = character(1)
      )

      write.csv(
        enrich_res,
        file = file.path(ct_dir, "pathway_enrichment.csv"),
        row.names = FALSE
      )
    }
  }
}

if (nrow(skipped) > 0) {
  message(sprintf("\n%d celltype/stratum combination(s) skipped:", nrow(skipped)))
  for (i in seq_len(nrow(skipped))) {
    message(sprintf(
      "  - celltype = %s, stratum = %s: %s",
      skipped$celltype[i], skipped$stratum[i], skipped$reason[i]
    ))
  }

  write.csv(
    skipped,
    file = file.path(output_dir, "skipped_celltypes.csv"),
    row.names = FALSE
  )
} else {
  message("\nNo celltype/stratum combinations were skipped.")
}

message("Pathway-proportion linkage analysis completed.")
