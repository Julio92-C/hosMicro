#!/usr/bin/env Rscript
# polish_ge_heatmap.R — hospital_microbiome genetic-element hotspot heatmap
# (former manuscript Figure 4), regenerated on the corrected data.
#
# Same layout as the original figure: one column per sample (air | surface),
# three row blocks — ARG drug classes, MGE replicon families, VF categories —
# coloured by ln(reads per million bacterial reads + 1). Rows whose abundance
# differs between air and surface (Kruskal-Wallis, BH-adjusted across all rows)
# are starred (* < 0.05, ** < 0.01, *** < 0.001).
#
# Inputs:  <datasets_dir>/abri_kraken2_cleaned.csv (NX-filtered, microbial only)
#          Reports/Bracken/bracken_arranged.csv     (unfiltered Kraken2 counts;
#            taxid 2 = reads classified as Bacteria, the per-million scaling
#            factor — the NX-filtered copy zeroes that row in some batch samples)
# Outputs: <figures_dir>/heatmap/ge_class_heatmap.{png,tiff}
#          <datasets_dir>/heatmap/ge_class_heatmap_values.csv (+ _stats.csv)
# Run from the repo root:
#   Rscript projects/hospital_microbiome/polish_ge_heatmap.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})
`%||%` <- function(a, b) if (is.null(a)) b else a

script_dir <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a))) else "." })
cfg    <- yaml::read_yaml(file.path(script_dir, "config.yaml"))
root   <- cfg$project_root
ds_dir <- file.path(root, cfg$outputs$datasets_dir)
fig_dir <- file.path(root, cfg$outputs$figures_dir, "heatmap")
out_dir <- file.path(ds_dir, "heatmap")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
grp_cols <- unlist(cfg$colors$type)

meta <- read_csv(file.path(root, cfg$metadata$file), show_col_types = FALSE) |>
  filter(type %in% c("Air", "Surface")) |>
  select(sample, type)

# ---- Per-sample bacterial reads (scaling factor) ----------------------------
brk <- read_csv(file.path(root, "Reports", "Bracken", "bracken_arranged.csv"),
                show_col_types = FALSE)
bact <- brk |> filter(taxid == 2) |>
  select(any_of(meta$sample)) |>
  pivot_longer(everything(), names_to = "sample", values_to = "bact_reads") |>
  filter(bact_reads > 0)
smp <- inner_join(meta, bact, by = "sample") |> arrange(type, as.integer(gsub("\\D", "", sample)))

# ---- Gene hits -> row labels -------------------------------------------------
hits <- read_csv(file.path(ds_dir, "abri_kraken2_cleaned.csv"),
                 show_col_types = FALSE, guess_max = 1e5) |>
  filter(sample %in% smp$sample)

beta_lactam <- c("penam", "cephalosporin", "cephamycin", "carbapenem", "penem", "monobactam")
drug_class <- function(res) {
  vapply(strsplit(tolower(res), ";"), function(tok) {
    tok <- trimws(tok)
    tok[tok %in% beta_lactam] <- "beta-lactam"
    tok <- sub("^streptogramin_[ab]$", "streptogramin", tok)
    tok[tok == "glycylcycline"] <- "tetracycline"
    fam <- unique(tok[nzchar(tok)])
    if (length(fam) == 0) return(NA_character_)
    if (length(fam) > 1) return("MDR")
    lab <- gsub("_", " ", fam)
    paste0(toupper(substr(lab, 1, 1)), substr(lab, 2, nchar(lab)))
  }, character(1))
}
replicon_family <- function(gene) sub("\\(.*$", "", sub("_.*$", "", gene))
vf_category <- function(product) {
  # VFDB PRODUCT: "... [<VF name> (VF0571) - <Category> (VFC0301)] [...]"
  pat <- ".* - ([^\\[\\]()]+?) \\(VFC\\d+\\).*"
  out <- ifelse(grepl(pat, product, perl = TRUE),
                trimws(sub(pat, "\\1", product, perl = TRUE)), NA_character_)
  out
}

arg <- hits |> filter(DATABASE == "card") |>
  mutate(row = drug_class(RESISTANCE), block = "ARGs")
mge <- hits |> filter(DATABASE == "plasmidfinder") |>
  mutate(row = replicon_family(GENE), block = "MGEs")
vf  <- hits |> filter(DATABASE == "vfdb") |>
  mutate(row = vf_category(PRODUCT), block = "VFs")

long <- bind_rows(arg, mge, vf) |>
  filter(!is.na(row)) |>
  group_by(block, row, sample) |>
  summarise(reads = n_distinct(sequence), .groups = "drop") |>
  inner_join(smp, by = "sample") |>
  mutate(rpm = reads / bact_reads * 1e6, value = log1p(rpm))
write_csv(long, file.path(out_dir, "ge_class_heatmap_values.csv"))

# ---- Matrix (rows within block ordered by total abundance) -------------------
rows <- long |> group_by(block, row) |> summarise(tot = sum(rpm), .groups = "drop") |>
  mutate(block = factor(block, levels = c("ARGs", "MGEs", "VFs"))) |>
  arrange(block, desc(tot))
mat <- long |> select(block, row, sample, value) |>
  pivot_wider(names_from = sample, values_from = value, values_fill = 0)
mat <- mat[match(paste(rows$block, rows$row), paste(mat$block, mat$row)), ]
m <- as.matrix(mat[, intersect(smp$sample, colnames(mat)), drop = FALSE])
missing_cols <- setdiff(smp$sample, colnames(m))
if (length(missing_cols)) {
  m <- cbind(m, matrix(0, nrow(m), length(missing_cols), dimnames = list(NULL, missing_cols)))
}
m <- m[, smp$sample, drop = FALSE]

# ---- Air vs surface tests per row --------------------------------------------
grp <- smp$type[match(colnames(m), smp$sample)]
p <- apply(m, 1, function(v) tryCatch(kruskal.test(v, factor(grp))$p.value, error = function(e) NA))
padj <- p.adjust(p, method = "BH")
stars <- ifelse(is.na(padj), "", ifelse(padj < 0.001, "***", ifelse(padj < 0.01, "**", ifelse(padj < 0.05, "*", ""))))
labels <- paste0(rows$row, stars)
rownames(m) <- make.unique(labels)
write_csv(data.frame(block = rows$block, row = rows$row, p = p, p_adj = padj, stars = stars),
          file.path(out_dir, "ge_class_heatmap_stats.csv"))

# ---- Draw ---------------------------------------------------------------------
ann_col <- data.frame(Group = grp, row.names = colnames(m))
ann_row <- data.frame(Class = as.character(rows$block), row.names = rownames(m))
ann_colors <- list(Group = grp_cols[c("Air", "Surface")],
                   Class = c(ARGs = "#E9B400", MGEs = "#1FBF4A", VFs = "#8E24AA"))
gaps_col <- sum(grp == "Air")
gaps_row <- cumsum(table(rows$block))[1:2]
pal <- colorRampPalette(c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D", "#67000D"))(100)

draw <- function(file, dev) {
  dev(file)
  pheatmap(m, color = pal, cluster_rows = FALSE, cluster_cols = FALSE,
           gaps_col = gaps_col, gaps_row = gaps_row,
           annotation_col = ann_col, annotation_row = ann_row,
           annotation_colors = ann_colors, annotation_names_row = FALSE,
           border_color = "grey70", cellwidth = 11, cellheight = 11,
           fontsize = 10, fontsize_row = 10, fontsize_col = 9, angle_col = 90,
           legend_breaks = pretty(c(0, max(m))), legend_labels = pretty(c(0, max(m))))
  invisible(dev.off())
}
w <- (ncol(m) * 11 + 420) / 72; h <- (nrow(m) * 11 + 200) / 72
draw(file.path(fig_dir, "ge_class_heatmap.png"),
     function(f) png(f, width = w, height = h, units = "in", res = 300))
draw(file.path(fig_dir, "ge_class_heatmap.tiff"),
     function(f) tiff(f, width = w, height = h, units = "in", res = 300, compression = "lzw"))

cat(sprintf("Heatmap: %d rows (ARGs %d, MGEs %d, VFs %d) x %d samples (air %d, surface %d)\n",
            nrow(m), sum(rows$block == "ARGs"), sum(rows$block == "MGEs"), sum(rows$block == "VFs"),
            ncol(m), sum(grp == "Air"), sum(grp == "Surface")))
cat("Significant rows (BH < 0.05):", paste(labels[nzchar(stars)], collapse = "; "), "\n")
cat("Wrote", file.path(fig_dir, "ge_class_heatmap.png"), "(+ .tiff)\n")
