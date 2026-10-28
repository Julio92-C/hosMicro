#!/usr/bin/env Rscript
# polish_beta_venn.R — hospital_microbiome community overview figure
# (preprint Figure S2 layout), regenerated from the pipeline outputs.
#
# Reproduces the collaborator's four-panel figure on the corrected data
# (batch-specific NX rule, microbial taxa only):
#   A  Venn (VennDiagram) of species detected in air vs surface samples
#   B  NMDS of Bray-Curtis dissimilarities on species relative abundances,
#      PERMANOVA + PERMDISP (9,999 permutations)
#   C  Venn of taxa (all ranks) carrying ARG / VF / MGE reads, air vs surface
#   D  PCoA of binary Jaccard dissimilarities on the presence/absence of
#      genetic elements per taxon per sample, PERMANOVA + PERMDISP
#
# Species panels (A, B) follow the paper's community rules: Kraken2 rank S,
# > filters.min_count_per_sample reads per sample, samples below
# taxa_filter.min_sample_reads excluded. Genetic-element panels (C, D) use
# every study sample with at least one hit, as in the legacy figure.
#
# Outputs: <figures_dir>/beta_diversity/beta_venn_overview.{png,tiff}, beta_venn_ge.{png,tiff}
#          <datasets_dir>/beta_diversity/beta_venn_overview_stats.csv
# Run from the repo root:
#   Rscript projects/hospital_microbiome/polish_beta_venn.R

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(vegan)
  library(VennDiagram)
  library(patchwork)
})
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
`%||%` <- function(a, b) if (is.null(a)) b else a

script_dir <- local({ a <- grep("^--file=", commandArgs(FALSE), value = TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a))) else "." })
cfg    <- yaml::read_yaml(file.path(script_dir, "config.yaml"))
root   <- cfg$project_root
ds_dir <- file.path(root, cfg$outputs$datasets_dir)
fig_dir <- file.path(root, cfg$outputs$figures_dir, "beta_diversity")
out_dir <- file.path(ds_dir, "beta_diversity")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cols   <- unlist(cfg$colors$type)          # Air grey, Surface pink
floor_n   <- cfg$filters$min_count_per_sample %||% 0
min_reads <- cfg$taxa_filter$min_sample_reads %||% 0
nperm     <- cfg$stats$permanova_permutations %||% 9999
group_lab <- cfg$panels$group_label %||% "Group"

meta <- read_csv(file.path(root, cfg$metadata$file), show_col_types = FALSE) |>
  filter(type %in% c("Air", "Surface")) |>
  select(sample, type)

# ---- Species table (panels A, B) ---------------------------------------------
sp <- read_csv(file.path(ds_dir, "noncontaminants_list.csv"),
               col_select = c("sample", "name", "tax_rank", "count"),
               show_col_types = FALSE) |>
  filter(tax_rank == "S", count > floor_n) |>
  inner_join(meta, by = "sample")
depth <- sp |> group_by(sample) |> summarise(n = sum(count), .groups = "drop")
sp <- filter(sp, sample %in% depth$sample[depth$n >= min_reads])

# ---- Genetic-element table (panels C, D) -------------------------------------
ge <- read_csv(file.path(ds_dir, "abri_kraken2_cleaned.csv"),
               col_select = c("sample", "name", "GENE", "DATABASE"),
               show_col_types = FALSE, guess_max = 1e5) |>
  filter(!is.na(name), !is.na(GENE)) |>
  inner_join(meta, by = "sample")

# ---- Helpers -------------------------------------------------------------------
venn_panel <- function(sets, title) {
  # draw.pairwise.venn mirrors the circles (and misplaces the labels) when
  # area1 < area2, so always pass the larger set first, then mirror back with
  # a 180-degree rotation when needed so Air is always drawn on the left
  # (`inverted` mirrors the circles but not the counts).
  ord <- names(sort(lengths(sets[c("Air", "Surface")]), decreasing = TRUE))
  g <- draw.pairwise.venn(
    area1 = length(sets[[ord[1]]]), area2 = length(sets[[ord[2]]]),
    cross.area = length(intersect(sets$Air, sets$Surface)),
    category = c("", ""), fill = unname(cols[ord]),
    rotation.degree = if (ord[1] == "Surface") 180 else 0,
    print.mode = c("raw", "percent"), sigdigs = 3,
    alpha = 0.6, lty = "blank", cex = 1.3, cat.cex = 1.4,
    fontfamily = "sans", cat.fontfamily = "sans",
    cat.pos = c(-20, 20), cat.dist = 0.05, scaled = FALSE, ind = FALSE)
  # Category labels drawn here: the rotation leaves VennDiagram's own labels
  # stacked at the top centre. Air is always the left circle.
  labs <- grid::gList(
    grid::textGrob("Air", x = 0.34, y = 0.90, gp = grid::gpar(fontsize = 17)),
    grid::textGrob("Surface", x = 0.66, y = 0.90, gp = grid::gpar(fontsize = 17)))
  wrap_elements(full = grid::gTree(children = grid::gList(g, labs))) + ggtitle(title)
}

ordinate <- function(mat, groups, method, binary, title, type = c("pcoa", "nmds")) {
  type <- match.arg(type)
  d <- vegdist(mat, method = method, binary = binary)
  set.seed(12345)
  if (type == "nmds") {
    nm   <- metaMDS(d, k = 2, trymax = 200, autotransform = FALSE, trace = FALSE)
    pts  <- scores(nm, display = "sites")
    axl  <- c("NMDS1", "NMDS2")
    note <- sprintf("\nNMDS stress = %.3f", nm$stress)
  } else {
    pc   <- cmdscale(d, k = 2, eig = TRUE)
    pts  <- pc$points
    expl <- 100 * pc$eig[1:2] / sum(pc$eig[pc$eig > 0])
    axl  <- sprintf(c("PCoA1 (%.1f%%)", "PCoA2 (%.1f%%)"), expl)
    note <- ""
  }
  set.seed(12345)
  pm <- adonis2(d ~ groups, permutations = nperm)
  set.seed(12345)
  pd <- permutest(betadisper(d, groups), permutations = nperm)$tab$`Pr(>F)`[1]
  df <- data.frame(X = pts[, 1], Y = pts[, 2], type = groups)
  fp  <- function(x) if (x < 0.001) "p < 0.001" else sprintf("p = %.3f", x)
  lab <- sprintf("PERMANOVA R² = %.3f, %s\nPERMDISP %s%s",
                 pm$R2[1], fp(pm$`Pr(>F)`[1]), fp(pd), note)
  # Normal-theory 95% ellipses (stat_ellipse's default multivariate-t hugs
  # the points more tightly); stats inset top-left.
  p <- ggplot(df, aes(X, Y, colour = type, fill = type)) +
    stat_ellipse(geom = "polygon", type = "norm", level = 0.95, alpha = 0.15,
                 show.legend = FALSE) +
    geom_point(size = 2.6, alpha = 0.9) +
    annotate("text", x = -Inf, y = Inf, label = lab, hjust = -0.04, vjust = 1.15,
             size = 3.8, lineheight = 1.05) +
    scale_colour_manual(values = cols, name = group_lab) +
    scale_fill_manual(values = cols, name = group_lab) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    labs(x = axl[1], y = axl[2], title = title) +
    theme_bw(base_size = 12) +
    theme(legend.position = "top", panel.grid.minor = element_blank(),
          legend.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15)) +
    guides(colour = guide_legend(override.aes = list(size = 4)))
  list(plot = p,
       stats = data.frame(r2 = pm$R2[1], p = pm$`Pr(>F)`[1], permdisp_p = pd,
                          stress = if (type == "nmds") nm$stress else NA_real_,
                          n_air = sum(groups == "Air"), n_surface = sum(groups == "Surface")))
}

to_matrix <- function(long, feature, value) {
  w <- long |> group_by(sample, feature = .data[[feature]]) |>
    summarise(v = sum(.data[[value]]), .groups = "drop") |>
    pivot_wider(names_from = feature, values_from = v, values_fill = 0)
  m <- as.matrix(w[, -1]); rownames(m) <- w$sample; m
}

# ---- A: species Venn -----------------------------------------------------------
sets_sp <- split(sp$name, sp$type) |> lapply(unique)
pA <- venn_panel(sets_sp, "A")

# ---- B: species Bray-Curtis NMDS ----------------------------------------------
m_sp <- to_matrix(sp, "name", "count")
m_sp <- m_sp / rowSums(m_sp)
grp_sp <- factor(meta$type[match(rownames(m_sp), meta$sample)], levels = c("Air", "Surface"))
oB <- ordinate(m_sp, grp_sp, "bray", FALSE, "B", type = "nmds")

# ---- C: taxa carrying genetic elements Venn -----------------------------------
sets_ge <- split(ge$name, ge$type) |> lapply(unique)
pC <- venn_panel(sets_ge, "C")

# ---- D: binary Jaccard on taxon x element presence/absence ---------------------
ge_pa <- ge |> distinct(sample, type, feature = paste(name, GENE, sep = " | ")) |>
  mutate(present = 1)
m_ge <- to_matrix(ge_pa, "feature", "present")
grp_ge <- factor(meta$type[match(rownames(m_ge), meta$sample)], levels = c("Air", "Surface"))
oD <- ordinate(m_ge, grp_ge, "jaccard", TRUE, "D", type = "pcoa")

# ---- Compose + save -------------------------------------------------------------
fig <- (pA | oB$plot) / (pC | oD$plot)
ggsave(file.path(fig_dir, "beta_venn_overview.png"), fig, width = 13, height = 11, dpi = 300)
ggsave(file.path(fig_dir, "beta_venn_overview.tiff"), fig, width = 13, height = 11, dpi = 300,
       compression = "lzw")

# Genetic-element panels alone (hospital Fig S3 once the species panels moved
# into Figure 1, polish_taxa_overview.R), relabelled A-B.
fig_ge <- venn_panel(sets_ge, "A") | (oD$plot + labs(title = "B"))
ggsave(file.path(fig_dir, "beta_venn_ge.png"), fig_ge, width = 13, height = 5.8, dpi = 300)
ggsave(file.path(fig_dir, "beta_venn_ge.tiff"), fig_ge, width = 13, height = 5.8, dpi = 300,
       compression = "lzw")

venn_counts <- function(s) c(air_only = length(setdiff(s$Air, s$Surface)),
                             surface_only = length(setdiff(s$Surface, s$Air)),
                             shared = length(intersect(s$Air, s$Surface)))
vA <- venn_counts(sets_sp); vC <- venn_counts(sets_ge)
stats <- bind_rows(
  data.frame(panel = "A", what = "species Venn", t(vA), total = sum(vA)),
  data.frame(panel = "C", what = "GE-carrying taxa Venn", t(vC), total = sum(vC)),
  cbind(data.frame(panel = "B", what = "species Bray-Curtis NMDS"), oB$stats),
  cbind(data.frame(panel = "D", what = "GE presence/absence binary Jaccard PCoA"), oD$stats))
write_csv(stats, file.path(out_dir, "beta_venn_overview_stats.csv"))

cat(sprintf("A species: air-only %d, surface-only %d, shared %d (%.1f%% of %d)\n",
            vA[1], vA[2], vA[3], 100 * vA[3] / sum(vA), sum(vA)))
cat(sprintf("B Bray-Curtis NMDS: n = %d air / %d surface, R2 = %.3f, p = %.4f, PERMDISP p = %.4f, stress = %.3f\n",
            oB$stats$n_air, oB$stats$n_surface, oB$stats$r2, oB$stats$p, oB$stats$permdisp_p, oB$stats$stress))
cat(sprintf("C GE taxa: air-only %d, surface-only %d, shared %d (%.1f%% of %d)\n",
            vC[1], vC[2], vC[3], 100 * vC[3] / sum(vC), sum(vC)))
cat(sprintf("D Jaccard: n = %d air / %d surface, R2 = %.3f, p = %.4f, PERMDISP p = %.4f\n",
            oD$stats$n_air, oD$stats$n_surface, oD$stats$r2, oD$stats$p, oD$stats$permdisp_p))
cat("Wrote", file.path(fig_dir, "beta_venn_overview.png"), "(+ .tiff)\n")
