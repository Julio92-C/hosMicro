#!/usr/bin/env Rscript
# polish_relative_abundance.R — project-specific relative-abundance figures
# for the hospital paper, replacing the fig05 composite (whose 50-species
# legend is crushed by the project-wide composite font/legend knobs).
#
# Emits two standalone per-sample stacked-bar figures, faceted Air | Surface,
# laid out after the collaborator's reference genus figure:
#   relative_abundance_genus_top10.{png,tiff}     top 10 genera PER GROUP
#                                                 (Air U Surface) + Other
#   relative_abundance_species_top30.{png,tiff}   top 30 species + Other
# plus the plotted tables under <datasets_dir>/relative_abundance/.
#
# Standalone and read-only with respect to the pipeline: reads the
# pipeline's noncontaminants_list.csv + this project's config.yaml and only
# SOURCES R/utils_taxa.R (a single pure function, clean_taxa_names) so the
# species filter is byte-identical to R/05. No pipeline file is modified and
# no pipeline output is overwritten.
#
# Method notes
#   * Genus figure: Kraken2 genus-rank (G) clade counts; genera ranked by
#     mean relative abundance within each group, top 10 of Air and top 10 of
#     Surface combined (reference-figure approach). Y = proportion (0-1).
#   * Species rows, count floor (count > cfg$filters$min_count_per_sample)
#     and per-sample % follow R/05_relative_abundance.R.
#   * Species-rank rows only (Kraken2 rank code "S"): strain / subgenus /
#     group rows carry clade counts that would double-count their species.
#   * Genus = the rank-G ancestor in the Kraken2 lineage (the pipeline's
#     `genus` column is empty for this project); genus % = sum of its
#     species' %, so both figures share the same denominator.
#   * Top-N ranked by mean relative abundance across samples (not total
#     read count, which would favour the deepest samples).
#   * Microbial taxa only: each taxid's lineage is rebuilt from the original
#     tree-ordered Kraken2 reports (metaOmics/Datasets/*_kraken2_db1_combined_
#     reports.txt). Kept: Bacteria, Archaea, Viruses, and Eukaryota outside
#     Embryophyta (land plants) and Metazoa (animals, incl. human) — i.e.
#     fungi, protists and microalgae stay. Excluded taxa are listed in
#     <datasets_dir>/relative_abundance/excluded_non_microbial_taxa.csv.
#   * Colours: 10 colour-blind-aware genus hues (Paul Tol "muted" + orange);
#     in the species figure every species takes a shade of its genus' hue, so
#     the two figures read as one system. "Other" is light grey.
#
# Run:  Rscript projects/hospital_microbiome/polish_relative_abundance.R

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
  library(ggplot2)
})
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- Resolve paths (CWD-independent) -----------------------------------------
.args      <- commandArgs(trailingOnly = FALSE)
.file_arg  <- sub("^--file=", "", .args[grep("^--file=", .args)])
script_dir <- dirname(normalizePath(.file_arg))
repo_root  <- Sys.getenv("MERIDIAN_HOME", "../meridian")   # meridian checkout (R/utils_taxa.R)
cfg        <- yaml::read_yaml(file.path(script_dir, "config.yaml"))
source(file.path(repo_root, "R", "utils_taxa.R"))   # clean_taxa_names()

root    <- cfg$project_root
ds_dir  <- file.path(root, cfg$outputs$datasets_dir)
fig_dir <- file.path(root, cfg$outputs$figures_dir, "relative_abundance")
out_ds  <- file.path(ds_dir, "relative_abundance")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_ds,  recursive = TRUE, showWarnings = FALSE)

grp  <- cfg$metadata$group_cols[[1]]                 # "type"
sid  <- cfg$metadata$sample_id_col
grp_levels <- names(unlist(cfg$colors[[grp]])) %||% c("Air", "Surface")

# ---- Knobs ---------------------------------------------------------------------
top_species  <- 30
top_genus    <- 10
other_label  <- "Other"
width_in     <- 16.5     # same canvas width as the fig02-04 composites, so
height_in    <- 7.5      #   Word shrinks every figure by the same factor
base_pt      <- 16       # axis ticks + axis titles + facet strips
title_pt     <- 18       # "Sample" / "Relative abundance"
legend_pt    <- 14       # legend entries + title
sample_pt    <- 11       # 59 sample IDs (largest that still fits the bars)
genus_pal <- c("#332288", "#88CCEE", "#44AA99", "#117733", "#999933",
               "#DDCC77", "#CC6677", "#882255", "#AA4499", "#E69F00")
extra_pal <- c("#0072B2", "#D55E00", "#009E73", "#56B4E9", "#6A3D9A",
               "#B15928", "#1F78B4", "#33A02C", "#FB9A99", "#BC80BD",
               "#8DD3C7", "#F0E442")
other_col <- "#D9D9D9"
# Genera dropped as likely Kraken2 false positives (malaria parasites in
# hospital air/surface samples).
exclude_genera <- c("Plasmodium")

# ---- Species table (same filter as R/05, but full binomials kept) -------------
raw <- readr::read_csv(file.path(ds_dir, "noncontaminants_list.csv"),
                       show_col_types = FALSE,
                       col_types = readr::cols(.default = "c", count = "d"))
# ---- Microbial filter from the Kraken2 taxonomy tree ---------------------------
# Combined kreports list taxa in tree (pre-order) order with the name indented
# by 2 spaces per depth level, so a depth-indexed stack yields the lineage.
kraken_lineage_class <- function(paths) {
  out <- list()
  for (p in paths) {
    ln <- readLines(p, warn = FALSE)
    ln <- ln[!startsWith(ln, "#")]
    parts <- strsplit(ln, "\t")
    rank  <- vapply(parts, function(x) x[length(x) - 2L], character(1))
    taxid <- vapply(parts, function(x) x[length(x) - 1L], character(1))
    nm    <- vapply(parts, function(x) x[length(x)], character(1))
    depth <- (nchar(nm) - nchar(sub("^ +", "", nm))) %/% 2L
    nm    <- sub("^ +", "", nm)
    stack <- character(); rstack <- character()
    cls   <- character(length(nm))
    genus <- rep(NA_character_, length(nm))
    for (i in seq_along(nm)) {
      stack  <- c(stack[seq_len(depth[i])], nm[i])
      rstack <- c(rstack[seq_len(depth[i])], rank[i])
      g <- which(rstack == "G")
      if (length(g)) genus[i] <- stack[max(g)]
      cls[i] <- if (any(stack %in% c("Bacteria", "Archaea"))) "prokaryote"
        else if ("Viruses" %in% stack) "virus"
        else if ("Embryophyta" %in% stack) "plant"
        else if ("Metazoa" %in% stack) "animal"
        else if ("Eukaryota" %in% stack) "microbial_eukaryote"
        else "other"
    }
    out[[p]] <- data.frame(taxid = taxid, rank = rank, genus = genus,
                           class = cls, stringsAsFactors = FALSE)
  }
  d <- do.call(rbind, out)
  d[!duplicated(d$taxid), ]
}
kreports <- list.files(file.path(root, "Datasets"),
                       pattern = "_kraken2_db1_combined_reports\\.txt$",
                       full.names = TRUE)
stopifnot(length(kreports) > 0)
lineage   <- kraken_lineage_class(kreports)
microbial <- c("prokaryote", "virus", "microbial_eukaryote")
raw$tax_class <- lineage$class[match(raw$taxid, lineage$taxid)]
raw$tax_class[is.na(raw$tax_class)] <- "unmatched"
raw$tax_rank  <- lineage$rank[match(raw$taxid, lineage$taxid)]
raw$tax_genus <- lineage$genus[match(raw$taxid, lineage$taxid)]

min_count <- cfg$relative_abundance$min_count_per_sample %||%
             cfg$filters$min_count_per_sample %||% 5
meta <- readr::read_csv(file.path(root, cfg$metadata$file),
                        show_col_types = FALSE, col_types = readr::cols(.default = "c"))

# Shared filters for both levels: count floor, microbial lineage, excluded
# genera (likely Kraken2 false positives), study samples only.
excluded_log <- list()
prep <- function(df, level) {
  df <- dplyr::filter(df, !is.na(.data$name), .data$count > min_count)
  bad <- !(df$tax_class %in% microbial) | (df$genus %in% exclude_genera)
  if (any(bad)) {
    excluded_log[[level]] <<- df[bad, ] |>
      dplyr::mutate(reason = ifelse(.data$tax_class %in% microbial,
                                    "excluded genus (likely false positive)",
                                    paste0("non-microbial (", .data$tax_class, ")"))) |>
      dplyr::group_by(level = level, .data$name, .data$genus, .data$reason) |>
      dplyr::summarise(total_count = sum(.data$count), n_samples = dplyr::n(),
                       .groups = "drop")
  }
  cat(sprintf("[%s] kept %.1f%% of reads after microbial/genus filters\n", level,
              100 * sum(df$count[!bad]) / sum(df$count)))
  df[!bad, ] |>
    dplyr::inner_join(meta[, c(sid, grp)], by = c("sample" = sid)) |>
    dplyr::filter(.data[[grp]] %in% grp_levels)
}

# Species: rank "S" only. Kraken2 `_all` counts are clade totals, so keeping
# strain (S1, S2...) or subgenus/group rows too would double-count reads.
td <- cfg$taxa_display %||% list()
td$abbreviate_genus <- FALSE
td$max_name_len     <- 1000   # no truncation: keep distinct species distinct
sp <- clean_taxa_names(raw[!is.na(raw$tax_rank) & raw$tax_rank == "S", ],
                       modifyList(cfg, list(taxa_display = td)))
first_word <- sub("\\s.*$", "", sp$name)
sp$genus <- ifelse(!is.na(sp$tax_genus), sp$tax_genus, first_word)
sp <- prep(sp, "species")

# Genus: rank "G" rows (genus clade counts, so reads classified only to genus
# level are included — same basis as the collaborator's genus figure).
gn <- raw[!is.na(raw$tax_rank) & raw$tax_rank == "G", ]
gn$genus <- gn$name
gn <- prep(gn, "genus")

# Shallow samples: same rule as the pipeline's taxa stages
# (taxa_filter.min_sample_reads, R/utils_taxa.R::community_filter) - summed
# species-rank reads over taxa above the count floor. Applied to both figures
# so genus and species panels show the same samples.
min_sample_reads <- cfg$taxa_filter$min_sample_reads %||% 0
if (min_sample_reads > 0) {
  depth   <- tapply(sp$count, sp$sample, sum)
  shallow <- names(depth)[depth < min_sample_reads]
  shallow <- union(shallow, setdiff(unique(gn$sample), names(depth)))
  if (length(shallow)) {
    cat(sprintf("[samples] dropped %d with < %d species reads: %s\n",
                length(shallow), min_sample_reads,
                paste(sort(shallow), collapse = ", ")))
  }
  sp <- sp[!sp$sample %in% shallow, ]
  gn <- gn[!gn$sample %in% shallow, ]
}

# ---- Aggregate + top-N -----------------------------------------------------------
# Relative abundance as a proportion of each sample's total at that level.
# Ranking uses the mean relative abundance across samples — per group when
# per_group = TRUE (top N of Air  U  top N of Surface, as in the reference
# figure), otherwise across all samples.
rel_table <- function(df, level_col, top_n, per_group) {
  agg <- df |>
    dplyr::group_by(.data$sample, .data[[grp]], taxon = .data[[level_col]]) |>
    dplyr::summarise(count = sum(.data$count), .groups = "drop") |>
    dplyr::group_by(.data$sample) |>
    dplyr::mutate(prop = .data$count / sum(.data$count)) |>
    dplyr::ungroup()
  n_by_grp <- agg |> dplyr::distinct(.data$sample, .data[[grp]]) |>
    dplyr::count(.data[[grp]], name = "n_samp")
  grp_mean <- agg |>
    dplyr::group_by(.data[[grp]], .data$taxon) |>
    dplyr::summarise(s = sum(.data$prop), .groups = "drop") |>
    dplyr::left_join(n_by_grp, by = grp) |>
    dplyr::mutate(mean_prop = .data$s / .data$n_samp)
  overall <- agg |>
    dplyr::group_by(.data$taxon) |>
    dplyr::summarise(mean_prop = sum(.data$prop) / dplyr::n_distinct(agg$sample),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$mean_prop))
  if (per_group) {
    top <- character()
    for (g in grp_levels) {
      r <- grp_mean[grp_mean[[grp]] == g, ]
      top <- c(top, utils::head(r$taxon[order(-r$mean_prop)], top_n))
    }
    top <- unique(top)
  } else {
    top <- utils::head(overall$taxon, top_n)
  }
  plot_df <- agg |>
    dplyr::mutate(taxon = ifelse(.data$taxon %in% top, .data$taxon, other_label)) |>
    dplyr::group_by(.data$sample, .data[[grp]], .data$taxon) |>
    dplyr::summarise(count = sum(.data$count), prop = sum(.data$prop),
                     .groups = "drop")
  list(plot_df = plot_df, top = top, rank = overall, grp_mean = grp_mean)
}

gen <- rel_table(gn, "genus", top_genus,   per_group = TRUE)
spc <- rel_table(sp, "name",  top_species, per_group = FALSE)

# ---- Colours -------------------------------------------------------------------
mix <- function(a, b, w) grDevices::colorRampPalette(c(a, b))(101)[round(w * 100) + 1]
all_pal <- c(genus_pal, extra_pal)

genus_of  <- setNames(sp$genus, sp$name)[!duplicated(sp$name)]
sp_genera <- unique(genus_of[spc$top])
extra_genera <- setdiff(sp_genera, gen$top)
genus_order  <- c(gen$top, extra_genera)
genus_cols   <- setNames(rep_len(all_pal, length(genus_order)), genus_order)

# Species legend/stack order: grouped by genus (genus-figure order first),
# most abundant species first within each genus; each species is a shade of
# its genus' colour (wide dark->light ramp so sibling species stay distinct).
sp_top <- tibble::tibble(taxon = spc$top, genus = unname(genus_of[spc$top]),
                         mean_prop = spc$rank$mean_prop[match(spc$top, spc$rank$taxon)]) |>
  dplyr::mutate(genus = factor(.data$genus, levels = genus_order)) |>
  dplyr::arrange(.data$genus, dplyr::desc(.data$mean_prop))
sp_split <- split(sp_top$taxon, sp_top$genus, drop = TRUE)
species_cols <- unlist(lapply(sp_split, function(tx) {
  base <- genus_cols[[as.character(genus_of[[tx[1]]])]]
  k <- length(tx)
  if (k == 1) return(base)
  grDevices::colorRampPalette(c(mix(base, "black", 0.5), base,
                                mix(base, "white", 0.7)))(k)
}), use.names = FALSE)
names(species_cols) <- unlist(sp_split, use.names = FALSE)

# ---- Labels: italic taxon names (full genus), "sp." kept upright ----------------
pm_label <- function(shown) {
  esc <- function(s) gsub('(["\\\\])', "\\\\\\1", s)
  out <- ifelse(
    shown == other_label, sprintf('"%s"', other_label),
    ifelse(grepl("\\ssp\\.", shown),
           sprintf('italic("%s")~"%s"', esc(sub("\\s.*$", "", shown)),
                   esc(sub("^\\S+\\s+", "", shown))),
           sprintf('italic("%s")', esc(shown))))
  parse(text = out)
}

# ---- Plot (layout after the collaborator's reference figure) --------------------
sample_levels <- function(s) {
  s <- unique(s)
  s[order(gsub("[0-9]+", "", s), suppressWarnings(as.numeric(gsub("\\D", "", s))))]
}

make_plot <- function(plot_df, levels_ord, cols, legend_title, ncol) {
  plot_df$taxon  <- factor(plot_df$taxon, levels = levels_ord)
  plot_df$sample <- factor(plot_df$sample, levels = sample_levels(plot_df$sample))
  plot_df[[grp]]  <- factor(plot_df[[grp]], levels = grp_levels)
  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$sample, y = .data$prop,
                                        fill = .data$taxon)) +
    ggplot2::geom_col(width = 0.9) +
    ggplot2::facet_grid(cols = ggplot2::vars(.data[[grp]]),
                        scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = cols, breaks = levels_ord,
                               labels = pm_label(levels_ord), name = legend_title) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.01)),
                                breaks = seq(0, 1, 0.25),
                                labels = function(x) sprintf("%.2f", x)) +
    ggplot2::labs(x = "Sample", y = "Relative abundance") +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = ncol)) +
    ggplot2::theme_bw(base_size = base_pt) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      axis.text.x       = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                hjust = 1, size = sample_pt),
      axis.text.y       = ggplot2::element_text(size = base_pt),
      axis.title        = ggplot2::element_text(size = title_pt, face = "bold"),
      strip.background  = ggplot2::element_rect(fill = "grey85", colour = "grey20"),
      strip.text        = ggplot2::element_text(face = "bold", size = title_pt),
      panel.spacing.x   = grid::unit(8, "pt"),
      legend.title      = ggplot2::element_text(size = legend_pt, face = "bold"),
      legend.text       = ggplot2::element_text(size = legend_pt),
      legend.key.height = grid::unit(0.65, "cm"),
      legend.key.width  = grid::unit(0.6, "cm"),
      legend.key.spacing.y = grid::unit(2, "pt"),
      legend.box.spacing = grid::unit(6, "pt"),
      plot.margin       = ggplot2::margin(6, 10, 6, 6)
    )
}

save_both <- function(p, stem) {
  png  <- file.path(fig_dir, paste0(stem, ".png"))
  tiff <- file.path(fig_dir, paste0(stem, ".tiff"))
  ggplot2::ggsave(png, p, width = width_in, height = height_in, dpi = 300,
                  bg = "white")
  ggplot2::ggsave(tiff, p, width = width_in, height = height_in, dpi = 300,
                  bg = "white", device = grDevices::tiff, compression = "lzw")
  cat(sprintf("Wrote %s (+ .tiff)\n", png))
}

# Genus: top 10 per group (union)
g_levels <- c(gen$top, if (other_label %in% gen$plot_df$taxon) other_label)
g_cols   <- c(genus_cols[gen$top], setNames(other_col, other_label))
p_gen <- make_plot(gen$plot_df, g_levels, g_cols, "Genus", ncol = 1)
save_both(p_gen, "relative_abundance_genus_top10")
cat(sprintf("Genus legend: %d genera (top %d per group, union) + %s\n",
            length(gen$top), top_genus, other_label))

# Species: top 30 overall
s_levels <- c(sp_top$taxon, if (other_label %in% spc$plot_df$taxon) other_label)
s_cols   <- c(species_cols, setNames(other_col, other_label))
p_spc <- make_plot(spc$plot_df, s_levels, s_cols, "Species", ncol = 2)
save_both(p_spc, "relative_abundance_species_top30")

# ---- Plotted tables + keys ------------------------------------------------------
readr::write_csv(dplyr::bind_rows(excluded_log) |>
                   dplyr::arrange(.data$level, dplyr::desc(.data$total_count)),
                 file.path(out_ds, "excluded_taxa.csv"))
readr::write_csv(gen$plot_df, file.path(out_ds, "genus_top10_per_group_long.csv"))
readr::write_csv(gen$grp_mean |> dplyr::select(-"s") |>
                   dplyr::arrange(.data[[grp]], dplyr::desc(.data$mean_prop)),
                 file.path(out_ds, "genus_mean_relative_abundance_by_group.csv"))
readr::write_csv(spc$plot_df, file.path(out_ds, "species_top30_long.csv"))
readr::write_csv(
  dplyr::bind_rows(
    tibble::tibble(level = "genus", taxon = names(g_cols), colour = unname(g_cols)),
    tibble::tibble(level = "species", taxon = s_levels,
                   genus = c(as.character(sp_top$genus), NA)[seq_along(s_levels)],
                   colour = unname(s_cols[s_levels]))),
  file.path(out_ds, "relative_abundance_colour_key.csv"))
cat("Done.\n")
