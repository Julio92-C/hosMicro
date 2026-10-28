#!/usr/bin/env Rscript
# polish_table1.R — hospital_microbiome Table 1: ARG / VF / MGE summary across
# air and surface samples, built only from pipeline outputs.
#
# Replaces the legacy Table 1 (GEs_statsAnova_gTable.R), whose "Samples
# (air/surface)" column counted gene x taxon combinations above per-class read
# thresholds and whose ANOVA used genes, not samples, as the unit. Here every
# count is per sample and the tests are the pipeline's per-sample ones:
#   - Unique genes, reads (unique Abricate-hit reads), samples with >= 1 hit
#   - Genes unique to air / unique to surface / shared
#   - Richness Kruskal-Wallis p_adj (BH, pipeline stats/ge_alpha_kw_padj_summary)
#   - PERMANOVA R2 and p (BH across the three classes) + PERMDISP p, from the
#     resistome / virulome / mobilome stage outputs
#
# Outputs: <datasets_dir>/panels/table1_genetic_elements.{csv,md}
# Run from the repo root: Rscript projects/hospital_microbiome/polish_table1.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

script_dir <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a))) else "." })
cfg    <- yaml::read_yaml(file.path(script_dir, "config.yaml"))
root   <- cfg$project_root
ds_dir <- file.path(root, cfg$outputs$datasets_dir)
out_dir <- file.path(ds_dir, "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta <- read_csv(file.path(root, cfg$metadata$file), show_col_types = FALSE) |>
  filter(type %in% c("Air", "Surface"))

# Denominators: study samples that carry sequencing data after cleaning.
nc_samples <- unique(read_csv(file.path(ds_dir, "noncontaminants_list.csv"),
                              col_select = "sample", show_col_types = FALSE)$sample)
seq_meta <- filter(meta, sample %in% nc_samples)
n_air  <- sum(seq_meta$type == "Air")
n_surf <- sum(seq_meta$type == "Surface")

hits <- read_csv(file.path(ds_dir, "abri_kraken2_cleaned.csv"),
                 show_col_types = FALSE, guess_max = 1e5) |>
  inner_join(select(seq_meta, sample, type), by = "sample")

classes <- tibble::tribble(
  ~class, ~db,             ~family,     ~stage,
  "ARGs", "card",          "Resistome", "resistome",
  "VFs",  "vfdb",          "Virulome",  "virulome",
  "MGEs", "plasmidfinder", "Mobilome",  "mobilome"
)

kw <- read_csv(file.path(ds_dir, "stats", "ge_alpha_kw_padj_summary.csv"),
               show_col_types = FALSE)

parse_permanova <- function(stage) {
  txt <- readLines(file.path(ds_dir, stage, "beta_permanova.txt"))
  m <- regmatches(txt, regexec("^Model\\s+\\d+\\s+\\S+\\s+(\\S+)\\s+\\S+\\s+(\\S+)", txt))
  m <- m[lengths(m) > 0][[1]]
  d <- readLines(file.path(ds_dir, stage, "beta_permdisp.txt"))
  g <- regmatches(d, regexec("^Groups\\s+\\d+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+(\\S+)", d))
  g <- g[lengths(g) > 0][[1]]
  c(r2 = as.numeric(m[2]), p = as.numeric(m[3]), disp_p = as.numeric(g[2]))
}

fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3))

rows <- lapply(seq_len(nrow(classes)), function(i) {
  cl <- classes[i, ]
  h  <- filter(hits, DATABASE == cl$db)
  air  <- unique(h$GENE[h$type == "Air"])
  surf <- unique(h$GENE[h$type == "Surface"])
  pm <- parse_permanova(cl$stage)
  tibble(
    class           = cl$class,
    unique_genes    = n_distinct(h$GENE),
    reads           = n_distinct(h$sequence),
    samples_air     = n_distinct(h$sample[h$type == "Air"]),
    samples_surface = n_distinct(h$sample[h$type == "Surface"]),
    genes_air_only  = length(setdiff(air, surf)),
    genes_surf_only = length(setdiff(surf, air)),
    genes_shared    = length(intersect(air, surf)),
    richness_kw_padj = kw$p_adj[kw$key == cl$family & kw$family == "richness"],
    shannon_kw_padj  = kw$p_adj[kw$key == cl$family & kw$family == "shannon"],
    permanova_r2    = pm[["r2"]],
    permanova_p     = pm[["p"]],
    permdisp_p      = pm[["disp_p"]]
  )
})
tab <- bind_rows(rows) |>
  mutate(permanova_padj = p.adjust(permanova_p, method = "BH"))

write_csv(tab, file.path(out_dir, "table1_genetic_elements.csv"))

md <- c(
  sprintf("| Element class | Unique genes | Reads | Samples with hits: air (n = %d) / surface (n = %d) | Genes air-only / surface-only / shared | Richness KW p_adj | Shannon KW p_adj | PERMANOVA R² (p_adj) | PERMDISP p |",
          n_air, n_surf),
  "|---|---|---|---|---|---|---|---|---|",
  sprintf("| %s | %d | %s | %d / %d | %d / %d / %d | %s | %s | %.3f (%s) | %s |",
          tab$class, tab$unique_genes, format(tab$reads, big.mark = ","),
          tab$samples_air, tab$samples_surface,
          tab$genes_air_only, tab$genes_surf_only, tab$genes_shared,
          fmt_p(tab$richness_kw_padj), fmt_p(tab$shannon_kw_padj),
          tab$permanova_r2, fmt_p(tab$permanova_padj), fmt_p(tab$permdisp_p))
)
writeLines(md, file.path(out_dir, "table1_genetic_elements.md"))
cat(md, sep = "\n")
cat("\nWrote", file.path(out_dir, "table1_genetic_elements.{csv,md}"), "\n")
