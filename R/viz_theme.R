#' Shared theme and colorblind-safe palettes for analysis figures
#'
#' Color is assigned by the job it does, following three
#' validated-by-construction, colorblind-safe choices:
#' - sequential magnitude  -> viridis (perceptually uniform, CVD-safe)
#' - signed / diverging     -> blue-grey-red centered at 0 (avoids red-green;
#'   grey neutral midpoint)
#' - categorical identity    -> Okabe-Ito (the canonical CVD-safe qualitative
#'   palette; used in fixed order, never cycled, max 8 entries)
#'
#' Cell types (46 of them) are NEVER encoded by color -- no palette survives
#' that many categories -- they are always put on an axis (heatmaps, facets,
#' sorted bars).
#'
#' @keywords internal
NULL


#' Okabe-Ito categorical palette (fixed order, colorblind-safe)
#'
#' @param n Number of colors needed (max 8).
#' @return Character vector of hex colors.
#' @export
okabe_ito_palette <- function(n = 8) {
  pal <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999"
  )
  if (n > length(pal)) {
    stop("okabe_ito_palette: only ", length(pal), " colorblind-safe categorical ",
         "colors are available; ", n, " requested. Encode the extra categories on ",
         "an axis instead of by color.", call. = FALSE)
  }
  pal[seq_len(n)]
}


#' Minimal recessive theme for analysis figures
#'
#' @param base_size Base font size.
#' @return A ggplot2 theme object.
#' @export
theme_analysis <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.3),
      axis.title = ggplot2::element_text(color = "grey20"),
      axis.text = ggplot2::element_text(color = "grey30"),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(color = "grey40"),
      legend.title = ggplot2::element_text(color = "grey20"),
      strip.text = ggplot2::element_text(face = "bold", color = "grey20")
    )
}


#' Diverging fill scale centered at zero (blue-grey-red)
#'
#' For signed quantities (log2 odds ratios, correlations): blue = negative,
#' red = positive, grey neutral at the midpoint.
#'
#' @param name Legend title.
#' @param limits Optional symmetric limits; if NULL, ggplot picks them.
#' @return A ggplot2 fill scale.
#' @export
scale_fill_diverging <- function(name = NULL, limits = NULL) {
  ggplot2::scale_fill_gradient2(
    name = name,
    low = "#2166AC", mid = "grey92", high = "#B2182B",
    midpoint = 0, limits = limits
  )
}


#' Diverging color scale centered at zero (blue-grey-red)
#'
#' @param name Legend title.
#' @param limits Optional symmetric limits.
#' @return A ggplot2 color scale.
#' @export
scale_color_diverging <- function(name = NULL, limits = NULL) {
  ggplot2::scale_color_gradient2(
    name = name,
    low = "#2166AC", mid = "grey70", high = "#B2182B",
    midpoint = 0, limits = limits
  )
}


#' Sequential viridis fill scale (magnitude)
#'
#' @param name Legend title.
#' @param ... Passed to `ggplot2::scale_fill_viridis_c`.
#' @return A ggplot2 fill scale.
#' @export
scale_fill_magnitude <- function(name = NULL, ...) {
  ggplot2::scale_fill_viridis_c(name = name, option = "viridis", ...)
}


#' Save a ggplot to both PDF and PNG under a figures directory
#'
#' @param plot A ggplot object.
#' @param filename Base filename (no extension).
#' @param output_dir Directory to write into.
#' @param width,height Dimensions in inches.
#' @return Invisibly, the PDF path.
#' @export
save_figure <- function(plot, filename, output_dir, width = 8, height = 6) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(output_dir, paste0(filename, ".pdf"))
  png_path <- file.path(output_dir, paste0(filename, ".png"))
  ggplot2::ggsave(pdf_path, plot, width = width, height = height)
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = 200)
  message("Wrote ", pdf_path, " and .png")
  invisible(pdf_path)
}
