
# Install and load required packages
install.packages("vegan")
install.packages("ggplot2")


# Load libraries
library(vegan)
library(ggplot2)
library(readr)
library(dplyr)
library(tidyverse)
library(plotly)
library(paletteer)

# Load dataset
# abri_kraken2_merged <- read_csv("OneDrive - University of West London/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/Datasets/abricate_kraken2_merged.csv")

## Load clean dataset
abri_kraken2_filtered <- read_csv("OneDrive - University of West London/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/Datasets/abri_kraken2_cleaned.csv")


# Import metadata
MetadataLocations <- read_csv("OneDrive - University of West London/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/Datasets/MetadataLocations2.csv")
colnames(MetadataLocations)


abri_kraken2_filtered <- abri_kraken2_filtered %>%
  arrange(GENE)
  

# Count unique number of genes
unique_gene <- unique(abri_kraken2_filtered$GENE)

# Annotation rows
# Pivot the Genes values into columns
df_wide_ann_rows <- abri_kraken2_filtered %>%
  arrange(GENE) %>%
  pivot_wider(names_from = GENE, values_from = c(name), values_fn = length) %>%
  group_by(sample) %>%
  summarise(across(16:293, ~ paste(., collapse = ", "))) %>%
  mutate(across(2:279, ~ str_replace_all(., "(NA,|,NA|,NA,|NA| )", "")))

colnames(df_wide_ann_rows)

# Assign the values of Sample as row names and empty as rows as NA value
df_wide_ann_rows <- df_wide_ann_rows %>%
  na.omit()  %>%
  remove_rownames %>%
  column_to_rownames(var="sample") %>%
  mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))


# Convert from a dataframe to a matrix
sample2arg_data_matrix <- data.matrix(df_wide_ann_rows, rownames.force = NA)

# Count NA in data sets
sum(is.na(sample2arg_data_matrix))

# Replace all NA values with 0 values
sample2arg_data_matrix <- sample2arg_data_matrix %>% replace(is.na(.), 0)

# Calculate Bray-Curtis distances
bray_curtis_dist <- vegdist(sample2arg_data_matrix, method = "bray")

# Perform PCoA
pcoa_result <- cmdscale(bray_curtis_dist, eig = TRUE, k = 2)

# Calculate percentage of variance explained
eig_values <- pcoa_result$eig
var_explained <- eig_values / sum(eig_values) * 100

# Create a data frame for plotting
pcoa_df <- data.frame(sample = rownames(sample2arg_data_matrix),
                      PC1 = pcoa_result$points[, 1],
                      PC2 = pcoa_result$points[, 2])


# Merge the metadata and sample dataset
ann_pcoa <- inner_join(pcoa_df, MetadataLocations, by="sample")

# Annotation col colors for Location
# Description variable
location_col <- as.factor(ann_pcoa$location)

# Generate a color palette based on the number of levels in gene_family
l <- length(levels(location_col))
l_palette <- paletteer_d("ggsci::default_igv", n = l)
l_palette

# Plot the results
pcaoa_plot <- ggplot(ann_pcoa, aes(x = PC1, y = PC2, label = sample, colour = location)) +
  geom_point(size = 3) +
  #geom_text(vjust = 1.5, hjust = 1.5) +
  #stat_ellipse(geom = "polygon", aes(fill = location), alpha = 0.2) +
  stat_ellipse(lwd = 0.8) + # Change line width
  scale_color_manual(values = l_palette) +
  #scale_fill_manual(values = l_palette) +
  labs(title = "",
       color = "Location",
       x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 2), "%)")
       ) +
  theme_classic() +
  theme(text = element_text(size = 16))

# ggplot graph
pcaoa_plot

plotly::ggplotly(pcaoa_plot)


# CLEAN UP #################################################

# Clear environment
rm(list = ls())

# Clear packages
detach("package:datasets", unload = TRUE)

# Clear plots
dev.off()  # But only if there IS a plot

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)
