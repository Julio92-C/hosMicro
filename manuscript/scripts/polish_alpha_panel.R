#!/usr/bin/env Rscript
# polish_alpha_panel.R — project-specific polish of the hospital paper's custom
# richness + shannon per-sample bar panel. Standalone: reads only the pipeline's
# saved alpha_diversity.csv and this project's config.yaml. Sources NO pipeline
# module (R/NN_*.R / utils) — keeps the polish decoupled from the pipeline.
#
# Layout is intentionally unchanged from the original composite (A = richness,
# B = shannon, each faceted Air | Surface, grey/pink bars, red overall-mean +
# black per-group-mean dashed lines, per-bar value labels). Only the styling is
# tidied: KW/BH stats moved from the bold title into a light subtitle, smaller
# value labels + rotated x-labels, consistent typography, bold A/B tags.
#
# Run:  Rscript projects/hospital_microbiome/polish_alpha_panel.R
# Spec: docs/superpowers/specs/2026-09-23-hospital-alpha-panel-polish-design.md

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ---- Resolve this script's own dir so CWD doesn't matter --------------------
.args      <- commandArgs(trailingOnly = FALSE)
.file_arg  <- sub("^--file=", "", .args[grep("^--file=", .args)])
script_dir <- dirname(normalizePath(.file_arg))
cfg        <- yaml::read_yaml(file.path(script_dir, "config.yaml"))

root    <- cfg$project_root
ds_dir  <- file.path(root, cfg$outputs$datasets_dir)
fig_dir <- file.path(root, cfg$outputs$figures_dir, "alpha_diversity")
csv     <- file.path(ds_dir, "alpha_diversity.csv")
out_png <- file.path(fig_dir, "alpha_richness_shannon_panel.png")   # overwrite

stopifnot(file.exists(csv))

# ---- Data + colours (single source of truth: cfg$colors$type) ---------------
grp  <- cfg$metadata$group_cols[[1]]                 # "type"
cols <- unlist(cfg$colors[[grp]])                    # Air/Surface -> hex
# The cfg$colors recolor lives on main; a branch reset before that commit
# won't have it. Fall back to the paper's Air-grey / Surface-pink scheme so
# this polish renders the intended colours regardless of branch state.
if (length(cols) == 0) cols <- c(Air = "#828385", Surface = "#ed8cb6")
div  <- readr::read_csv(csv, show_col_types = FALSE)

metrics <- c(richness = "Richness", shannon = "Shannon")

# Export stays 13x7in; Word shrinks a full-width image to ~50% of the text
# column, halving apparent font size. `word_scale` multiplies an explicit,
# BALANCED type hierarchy (below) so every element stays legible on the Word
# page while the title stays only slightly larger than the body — no runaway
# heading (the bug from blanket-scaling base_size). Bump word_scale to go
# larger/smaller uniformly; tweak individual `sz` entries to rebalance.
word_scale <- 2
sz <- list(
  title    = 8   * word_scale,   # metric heading (Richness / Shannon)
  subtitle = 6   * word_scale,   # KW / p_adj line
  ax_title = 7   * word_scale,   # "Sample" / "Richness"
  ax_text  = 6   * word_scale,   # y-axis ticks
  x_text   = 4.5 * word_scale,   # sample IDs
  strip    = 6.5 * word_scale,   # Air / Surface facet labels
  tag      = 9   * word_scale    # bold A / B
)

# ---- Family-adjust KW across the alpha indices (BH), as the pipeline (R/06) --
# R/06 adjusts over every index it reports (richness, shannon, simpson, chao1),
# so use the same family to match the text and Figure 1 (Shannon p_adj 0.47).
fam <- intersect(c("richness", "shannon", "simpson", "chao1"), colnames(div))
raw_p <- vapply(fam, function(m) {
  sub <- div[!is.na(div[[m]]), , drop = FALSE]
  suppressWarnings(kruskal.test(reformulate(grp, m), data = sub)$p.value)
}, numeric(1))
adj_p <- p.adjust(raw_p, method = "BH")

# Natural sample order within each panel (A2 before A10, S2 before S10).
div$sample <- factor(div$sample, levels = unique(div$sample[order(
  substr(div$sample, 1, 1), as.numeric(gsub("[^0-9]", "", div$sample)))]))

# ---- One tidied faceted bar panel per metric --------------------------------
make_panel <- function(m) {
  sub        <- div[!is.na(div[[m]]), , drop = FALSE]
  pretty     <- metrics[[m]]
  overall    <- mean(sub[[m]], na.rm = TRUE)
  grp_means  <- sub |>
    dplyr::group_by(.data[[grp]]) |>
    dplyr::summarise(grp_mean = mean(.data[[m]], na.rm = TRUE), .groups = "drop")
  subtitle   <- sprintf("KW p = %s   |   p_adj (BH) = %s",
                        format(raw_p[[m]], digits = 3),
                        format(adj_p[[m]], digits = 3))

  ggplot2::ggplot(sub, ggplot2::aes(x = sample, y = .data[[m]],
                                    fill = .data[[grp]])) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_hline(yintercept = overall, linetype = "dashed",
                        colour = "red", linewidth = 0.5) +
    ggplot2::geom_hline(data = grp_means,
                        ggplot2::aes(yintercept = .data$grp_mean),
                        linetype = "dashed", colour = "black", linewidth = 0.6) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::facet_wrap(stats::as.formula(paste("~", grp)),
                        scales = "free_x", nrow = 1) +
    ggplot2::labs(title = pretty, subtitle = subtitle, x = "Sample", y = pretty) +
    ggplot2::theme_classic(base_size = sz$ax_text) +
    ggplot2::theme(
      legend.position = "none",
      plot.title      = ggplot2::element_text(face = "bold", hjust = 0.5,
                                              size = sz$title),
      plot.subtitle   = ggplot2::element_text(hjust = 0.5, colour = "grey30",
                                              size = sz$subtitle),
      axis.title      = ggplot2::element_text(size = sz$ax_title),
      axis.text.y     = ggplot2::element_text(size = sz$ax_text),
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1,
                                              size = sz$x_text),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = NA),
      strip.text      = ggplot2::element_text(face = "bold", size = sz$strip),
      plot.margin     = ggplot2::margin(10, 8, 4, 6)
    )
}

panel <- make_panel("richness") / make_panel("shannon") +
  patchwork::plot_annotation(tag_levels = "A") &
  ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold",
                                                  size = sz$tag))

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(out_png, panel, width = 13, height = 7, dpi = 300)
cat(sprintf("Wrote %s\n", out_png))
