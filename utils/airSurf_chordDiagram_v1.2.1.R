# Load the libraries
library(readr)
library(circlize)
library(plotly)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(grid)
library(paletteer)

# Load dataset
# abri_kraken2_merged <- read_csv("hospital_Microbiome/data/abricate_kraken2_filtered.csv")

## Load clean dataset
abri_kraken2_merged <- read_csv("OneDrive - University of West London/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/Datasets/abri_kraken2_cleaned.csv")


# Filter and mutate the dataset
abri_kraken2_filtered <- abri_kraken2_merged %>%
filter(grepl("card", DATABASE)) %>%
  mutate(GENE = case_when(
    grepl("vanR_gene_in_vanA_cluster", GENE) ~ "vanRA",
    grepl("vanH_gene_in_vanA_cluster", GENE) ~ "vanHA",
    grepl("vanS_gene_in_vanA_cluster", GENE) ~ "vanSA",
    TRUE ~ GENE
  ),
  name = ifelse(grepl("Rahnella aquatilis CIP 78.65 = ATCC 33071", name), "Rahnella aquatilis", name)) %>%
  arrange(name) %>%
  mutate(RESISTANCE = ifelse(nchar(RESISTANCE) > 20, "Multi-drug", RESISTANCE))

# Define modify_strings function
modify_strings <- function(df, column_name, target_strings, replacement_strings) {
  if (length(target_strings) != length(replacement_strings)) {
    stop("target_strings and replacement_strings must be of the same length")
  }
  for (i in seq_along(target_strings)) {
    df <- df %>%
      mutate(!!sym(column_name) := str_replace_all(!!sym(column_name), target_strings[i], replacement_strings[i]))
  }
  return(df)
}

# List of target and replacement strings
target_strings <- c("cephalosporin;penam", "carbapenem;penam", "Cephalosporin;Penem")
replacement_strings <- rep("Beta-lactam", length(target_strings))

# Call Modify String function
abri_kraken2_filtered <- modify_strings(abri_kraken2_filtered, "RESISTANCE", target_strings, replacement_strings)

# Clean the RESISTANCE strings
abri_kraken2_filtered <- abri_kraken2_filtered %>%
  mutate(RESISTANCE = str_replace_all(RESISTANCE, "_", " ") %>% str_to_title())

# Filter by species
df <- abri_kraken2_filtered %>%
  filter(grepl("Pseudomonas aeruginosa|Enterococcus faecium|Staphylococcus aureus|Acinetobacter baumannii|Klebsiella pneumoniae|Klebsiella quasipneumoniae|Acinetobacter lwoffii|Enterobacter cloacae complex|Enterobacter hormaechei|Escherichia coli|Staphylococcus aureus|Rahnella aquatilis", name))

# Create a connection matrix
categories <- unique(c(sort(df$sample), sort(df$name), sort(df$GENE), sort(df$RESISTANCE)))
mat <- matrix(0, nrow = length(categories), ncol = length(categories), dimnames = list(categories, categories))

for (i in 1:nrow(df)) {
  mat[df$sample[i], df$name[i]] <- 1
  mat[df$name[i], df$GENE[i]] <- 1
  mat[df$GENE[i], df$RESISTANCE[i]] <- 1
}

# Set colors for categories
color_sectors <- c(
  paletteer_d("ggsci::default_igv", n = length(unique(df$sample))),
  paletteer_d("ggsci::default_igv", n = length(unique(df$name))),
  paletteer_d("ggsci::default_igv", n = length(unique(df$GENE))),
  paletteer_d("ggsci::default_igv", n = length(unique(df$RESISTANCE)))
)



# Graph parameters
circos.par(track.height = 0.1, 
           start.degree = 115, 
           gap.degree = 2, 
           canvas.xlim = c(-1, 1), 
           canvas.ylim = c(-1, 1), 
           circle.margin = c(1, 1), 
           unit.circle.segments = 500)

# Create the chord diagram
chordDiagram(mat, 
             transparency = 0.5, 
             annotationTrack = "grid", 
             scale = FALSE, directional = 1, 
             diffHeight = mm_h(3), 
             grid.col = color_sectors, 
             preAllocateTracks = list(track.height = 0.1, 
                                      unit.circle.segments = 100, 
                                      start.degree = 75, 
                                      scale = TRUE)
             )

# Add labels to the sectors
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], sector.index, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.8), cex = 1)
}, bg.border = NA)



## reset the graphic parameters and internal variables
circos.clear()

# Clear environment
rm(list = ls())

# Clear packages
detach("package:datasets", unload = TRUE)

# Clear plots
dev.off()  # But only if there IS a plot

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)