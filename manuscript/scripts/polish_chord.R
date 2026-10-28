#!/usr/bin/env Rscript
# polish_chord.R — hospital_microbiome Figure 3 (sample -> pathogen -> ARG ->
# drug class chord diagram), regenerated from the pipeline outputs.
#
# Port of the legacy `R scripts/airSurf_chordDiagram_v1.2.1.R`. The figure
# logic is unchanged (CARD hits only, "G. species" names, vanA-cluster gene
# relabels, >20-char resistance strings -> "Multi-drug", beta-lactam merges,
# the same pathogen list, igv palette, circos layout). Only the input moves:
# the legacy script read Datasets/abri_kraken2_merged1.csv (old global
# decontamination); this reads the pipeline's NX-filtered, microbial-only
# gene table <datasets_dir>/abri_kraken2_cleaned.csv, restricted to the
# metadata's Air/Surface samples.
#
# Outputs: <figures_dir>/network/chord_pathogens_legacy.{png,tiff} and
#          <datasets_dir>/network/chord_pathogens_legacy_links.csv
# Run from the repo root:
#   Rscript projects/hospital_microbiome/polish_chord.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(circlize)
  library(paletteer)
})

script_dir <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a))) else "." })
cfg  <- yaml::read_yaml(file.path(script_dir, "config.yaml"))
root <- cfg$project_root
ds_dir  <- file.path(root, cfg$outputs$datasets_dir)
fig_dir <- file.path(root, cfg$outputs$figures_dir, "network")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ds_dir, "network"), recursive = TRUE, showWarnings = FALSE)

meta <- read_csv(file.path(root, cfg$metadata$file), show_col_types = FALSE)
study_samples <- meta$sample[meta$type %in% c("Air", "Surface")]

hits <- read_csv(file.path(ds_dir, "abri_kraken2_cleaned.csv"),
                 show_col_types = FALSE, guess_max = 1e5) |>
  filter(sample %in% study_samples)

# ---- Legacy transformations (verbatim logic) ---------------------------------
df_all <- hits |>
  filter(grepl("card", DATABASE)) |>
  mutate(
    name = str_replace(name, "^(\\w)\\w+\\s(\\w+)", "\\1. \\2"),
    name = if_else(str_length(name) > 30,
                   str_replace(name, "^((?:\\S+\\s+){2}).*", "\\1"),
                   name),
    name = str_trim(name)) |>
  rename(SAMPLE = sample, NAME = name) |>
  mutate(GENE = case_when(
           grepl("vanR_gene_in_vanA_cluster", GENE) ~ "vanRA",
           grepl("vanH_gene_in_vanA_cluster", GENE) ~ "vanHA",
           grepl("vanS_gene_in_vanA_cluster", GENE) ~ "vanSA",
           TRUE ~ GENE),
         NAME = ifelse(grepl("Rahnella aquatilis CIP 78.65 = ATCC 33071", NAME),
                       "Rahnella aquatilis", NAME)) |>
  arrange(NAME) |>
  mutate(RESISTANCE = ifelse(nchar(RESISTANCE) > 20, "Multi-drug", RESISTANCE))

for (t in c("cephalosporin;penam", "carbapenem;penam", "Cephalosporin;Penem")) {
  df_all$RESISTANCE <- str_replace_all(df_all$RESISTANCE, t, "Beta-lactam")
}
df_all <- mutate(df_all, RESISTANCE = str_to_title(str_replace_all(RESISTANCE, "_", " ")))

pathogens <- "P. aeruginosa|E. faecium|S. aureus|A. baumannii|K. pneumoniae|K. quasipneumoniae|A. lwoffii|E. cloacae complex|E. hormaechei|E. coli|S. aureus|R. aquatilis"
df <- filter(df_all, grepl(pathogens, NAME))

links <- distinct(df, SAMPLE, NAME, GENE, RESISTANCE) |>
  left_join(distinct(meta, sample, type), by = c("SAMPLE" = "sample"))
write_csv(links, file.path(ds_dir, "network", "chord_pathogens_legacy_links.csv"))

# ---- Connection matrix + chord (legacy layout) --------------------------------
categories <- unique(c(sort(df$SAMPLE), sort(df$NAME), sort(df$GENE), sort(df$RESISTANCE)))
mat <- matrix(0, nrow = length(categories), ncol = length(categories),
              dimnames = list(categories, categories))
for (i in seq_len(nrow(df))) {
  mat[df$SAMPLE[i], df$NAME[i]]     <- 1
  mat[df$NAME[i], df$GENE[i]]       <- 1
  mat[df$GENE[i], df$RESISTANCE[i]] <- 1
}

color_sectors <- c(
  paletteer_d("ggsci::default_igv", n = length(unique(df$SAMPLE))),
  paletteer_d("ggsci::default_igv", n = length(unique(df$NAME))),
  paletteer_d("ggsci::default_igv", n = length(unique(df$GENE))),
  paletteer_d("ggsci::default_igv", n = length(unique(df$RESISTANCE)))
)

draw_chord <- function() {
  circos.clear()
  circos.par(track.height = 0.1, start.degree = 130, gap.degree = 2,
             canvas.xlim = c(-1.25, 1.25), canvas.ylim = c(-1.25, 1.25),  # legacy: +-1 clipped long labels
             circle.margin = c(1, 1), unit.circle.segments = 500)
  chordDiagram(mat, transparency = 0.5, annotationTrack = "grid", scale = FALSE,
               directional = 1, diffHeight = mm_h(3), grid.col = color_sectors,
               preAllocateTracks = list(track.height = 0.1,
                                        unit.circle.segments = 100,
                                        start.degree = 90, scale = TRUE))
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.8), cex = 1)
  }, bg.border = NA)
  circos.clear()
}

png(file.path(fig_dir, "chord_pathogens_legacy.png"), width = 12, height = 12,
    units = "in", res = 300)
draw_chord(); invisible(dev.off())
tiff(file.path(fig_dir, "chord_pathogens_legacy.tiff"), width = 12, height = 12,
     units = "in", res = 300, compression = "lzw")
draw_chord(); invisible(dev.off())

smp <- distinct(links, SAMPLE, type)
cat(sprintf("Chord: %d samples (%d air, %d surface), %d pathogens, %d ARGs, %d drug classes\n",
            nrow(smp), sum(smp$type == "Air"), sum(smp$type == "Surface"),
            n_distinct(df$NAME), n_distinct(df$GENE), n_distinct(df$RESISTANCE)))
cat("Pathogens:", paste(sort(unique(df$NAME)), collapse = ", "), "\n")
cat("Drug classes:", paste(sort(unique(df$RESISTANCE)), collapse = ", "), "\n")
cat("Wrote", file.path(fig_dir, "chord_pathogens_legacy.png"), "(+ .tiff)\n")
