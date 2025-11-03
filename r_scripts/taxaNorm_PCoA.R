
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
abri_kraken2_merged <- read_csv("Hospital_microbiome/Datasets/abri_kraken2_merged1.csv")


# Unique genetic elements taxa
uniqueTaxa <- length(unique(abri_kraken2_merged$name))


# Create a sample × taxa matrix from your microbiome_df
abundance_matrix <- abri_kraken2_merged %>%
  select(sample, name, count) %>%
  # mutate(
  #   count = (count / sum(count)) * 100 # Convert the count to relative abundance towards to normal distributions
  # ) %>%
  tidyr::pivot_wider(
    names_from = name,
    values_from = count,
    values_fill = 0,
    values_fn = sum # Use `sum` to aggregate duplicate values
  ) %>%
  # distinct(sample, .keep_all = TRUE) %>%
  column_to_rownames("sample")

# Convert the data frame to a matrix
abundance_matrix <- as.matrix(abundance_matrix)

# transform an abundance matrix into a relative abundance matrix
relative_abundance_matrix <- abundance_matrix %>%
  make_relative()

# Check how many missing values you have
sum(is.na(relative_abundance_matrix))
sum(is.nan(as.matrix(relative_abundance_matrix))) # Check for NaNs as well

# Remove rows with any missing values
clean_matrix <- na.omit(relative_abundance_matrix)

# Compute Bray-Curtis distance
bray_curtis_dist <- vegdist(clean_matrix, method = "bray")

# Perform PCoA
pcoa_result <- cmdscale(bray_curtis_dist, eig = TRUE, k = 2)

# Calculate percentage of variance explained
eig_values <- pcoa_result$eig
var_explained <- eig_values / sum(eig_values) * 100

# Create a data frame for plotting
pcoa_df <- data.frame(sample = rownames(clean_matrix),
                      PC1 = pcoa_result$points[, 1],
                      PC2 = pcoa_result$points[, 2])


# Metadata subset
MetadataLocations <- abri_kraken2_merged %>%
  select(sample, type, location, ward)

# # Merge the metadata and sample dataset
ann_pcoa <- inner_join(pcoa_df, MetadataLocations, by="sample")

# Re-arrange ann_pcoa
ann_pcoa_sorted <- ann_pcoa %>%
  distinct() %>%
  remove_rownames %>%
  column_to_rownames(var="sample")

# Check that the numeber of row match
nrow(clean_matrix) # Should be the number of samples you have
nrow(ann_pcoa)     # Should be the number of samples you have


# Run the PERMANOVA test
permanova_results <- adonis2(
  clean_matrix ~ type,          # Formula: distance matrix explained by 'Group'
  data = ann_pcoa_sorted,           # Data frame where 'Group' variable is found
  permutations = 999,           # Number of permutations (999 or 9999 for publication)
  method = "bray"               # Specify distance method if not already specified in the matrix
)

# Print the results
print(permanova_results)

# A Permutational Multivariate Analysis of Variance (PERMANOVA) using Bray-Curtis 
# distances indicated a significant difference in community composition across the 
# different levels of the 'type (Air vs Surface)' variable (F = 7.847, R² = 0.117, p = 0.001). The 
# variable explained approximately 11.7% of the total variation observed in 
# the community data."

# Extract R2 and P-value for the 'type' variable
# The 'type' variable is on the first row of results (after the 'Df' column, before 'Residual')
# The results are stored within the data frame structure of the permanova_results object
r2_value <- permanova_results$R2[1]
p_value <- permanova_results$`Pr(>F)`[1]

# Format the results into a single label string for the plot
# Use paste0 to combine text and round the values
permanova_label <- paste0(
  "PERMANOVA R\u00B2: ", round(r2_value, 3), 
  "\nP-value: ", round(p_value, 4)
)

# You might want to use scientific notation for very small P-values
# permanova_label <- paste0("PERMANOVA R\u00B2: ", round(r2_value, 3), "\nP-value: ", format.pval(p_value, digits = 3))


# # Annotation colors for Location
# # Description variable
# location_col <- as.factor(ann_pcoa$type)
# 
# # Generate a color palette based on the number of levels in gene_family
# l <- length(levels(location_col))
# l_palette <- paletteer_d("ggsci::default_igv", n = l)

# # Annotation colors for sample type
l_palette <- c(Air = "#828385", Surface = "#ed8cb6")

# Plot the results
pcaoa_plot <- ggplot(ann_pcoa, aes(x = PC1, y = PC2, label = sample, colour = type)) +
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
  theme(text = element_text(size = 16)) +
  # --- 3. Add the PERMANOVA results as text annotation ---
  annotate(
    "text", 
    x = Inf, y = Inf,             # Position the text at the top right corner
    label = permanova_label,      # The string we created earlier
    hjust = 1.1, vjust = 1.1,     # Adjust justification slightly inwards
    size = 4,                     # Control font size
    color = "black"               # Control font color
  )

# ggplot graph
pcaoa_plot

# Plot PCoA analysis with the plotly package
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
