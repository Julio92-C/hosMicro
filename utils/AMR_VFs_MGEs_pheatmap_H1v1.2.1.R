# Load the libraries
library(readr)
library(dplyr)
library(tidyverse)
library(paletteer)
library(reshape2)
library(pheatmap)
library(stringr)


# Load dataset
abri_kraken2_merged <- read_csv("OneDrive - University of West London/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/Datasets/abricate_kraken2_filtered.csv")

# Import metadata
MetadataLocations <- read_csv("OneDrive - University of West London/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/Datasets/MetadataLocations2.csv")
colnames(MetadataLocations)

# Pivot the Sample values into columns
df_wide <- abri_kraken2_merged %>%
  arrange(sample) %>%
  pivot_wider(names_from = sample, values_from = taxid) %>%
  group_by(name) %>%
  summarise(across(16:42, ~ paste(., collapse = ", "))) %>%
  mutate(across(2:28, ~ str_replace_all(., "(NA,|,NA|,NA,|NA|,| )", "")))

# Assign the values of col1 as row names and handle empty strings
processed_mgs2arg_data <- df_wide %>%
  remove_rownames() %>% 
  column_to_rownames(var = "name") %>%
  mutate(across(everything(), ~ na_if(., "")))

# Remove rows where all values are NA
sample2arg_data <- processed_mgs2arg_data %>%
  filter(rowSums(is.na(.)) < ncol(.))

# Convert from a dataframe to a matrix and replace NA with 0
sample2arg_data_matrix <- data.matrix(sample2arg_data, rownames.force = NA) %>%
  replace(is.na(.), 0)

# Count NA in data sets
sum(is.na(sample2arg_data_matrix))

# Convert to a binary matrix
sample2arg_data_Bmatrix <- as.matrix((sample2arg_data_matrix > 0) + 0)

# Order the matrix by row names
sample2arg_data_Bmatrix <- sample2arg_data_Bmatrix[order(rownames(sample2arg_data_Bmatrix)), ]

# Rearranging a Binary Matrix by Row Count
rearrange_matrix <- function(mat) {
  # Calculate the row sums (count of '1's)
  ordered_indices <- order(rowSums(mat), decreasing = TRUE)
  
  # Rearrange the matrix
  return(mat[ordered_indices, ])
}

# Rearranging the binary matrix
sorted_matrix <- rearrange_matrix(sample2arg_data_Bmatrix)

# Modify ordering of the clusters using clustering callback option
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

# Annotations col names
# Transpose sample2arg_data
sample2arg_data_transposed <- as.data.frame(t(sample2arg_data))

# Convert rownames into column and select sample
sample2arg_data_transposed <- sample2arg_data_transposed %>%
  mutate(sample = rownames(.)) %>%
  select(sample)

# Covert to a dataframe and trim whitespace
MetadataLocations <- data.frame(MetadataLocations)
MetadataLocations$sample <- str_trim(MetadataLocations$sample)

# Merge the metadata and sample dataset
ann_col <- inner_join(sample2arg_data_transposed, MetadataLocations, by = "sample") %>%
  na.omit() %>%
  remove_rownames() %>%
  column_to_rownames(var = "sample") %>%
  rename(LOCATION = location)

# Generate a color palette based on the number of levels in location
d_palette <- paletteer_d("ggsci::category10_d3", n = nlevels(as.factor(ann_col$LOCATION)))

# Create a named dataframe for the col colors
df4 <- ann_col %>%
  distinct(LOCATION) %>%
  mutate(color = d_palette)

# Create a named vector for rows
location_colors <- setNames(as.character(df4$color), df4$LOCATION)

# Annotation rows
# Pivot the Sample values into columns
df_wide_ann_rows <- abri_kraken2_merged %>%
  arrange(sample) %>%
  pivot_wider(names_from = DATABASE, values_from = c(GENE)) %>%
  group_by(name) %>%
  summarise(across(15:17, ~ paste(., collapse = ", "))) %>%
  mutate(across(2:4, ~ str_replace_all(., "(NA,|,NA|,NA,|NA| )", ""))) %>%
  na.omit() %>%
  remove_rownames() %>%
  column_to_rownames(var = "name") %>%
  mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))

# Classify functions
classify <- function(input_string, type) {
  if (is.na(input_string)) return(NA)
  elements <- unlist(strsplit(input_string, ","))
  unique_elements <- length(unique(elements))
  
  if (type == "virulence") {
    return(ifelse(unique_elements == 1, "Virulent", "Hypervirulent"))
  } else if (type == "amr") {
    return(ifelse(unique_elements == 1, "Drug-resistant", "Multidrug-resistant"))
  } else if (type == "MGE") {
    return(ifelse(unique_elements == 1, "Single MGE", "Multi-MGEs"))
  }
}

# Apply classification functions
ann_rows <- df_wide_ann_rows %>%
  mutate(VFs = sapply(vfdb, classify, type = "virulence"),
         AMR = sapply(card, classify, type = "amr"),
         MGEs = sapply(plasmidfinder, classify, type = "MGE")) %>%
  select(VFs, AMR, MGEs)

# Resistance category variable
Resistance_category <- abri_kraken2_merged %>%
  distinct(name, GENE, RESISTANCE) %>%
  mutate(RESISTANCE = ifelse(nchar(RESISTANCE) > 20, "Multi-drug", RESISTANCE))

# Modify strings function
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

# Call Modify String function
Resistance_category <- modify_strings(Resistance_category, "RESISTANCE", 
                                      c("cephalosporin;penam", "carbapenem;penam", "Cephalosporin;Penem"), 
                                      c("Beta-lactam", "Beta-lactam", "Beta-lactam"))

# Clean string function
clean_string <- function(x) {
  str_to_title(str_replace_all(x, "_", " "))
}

# Apply the function to the RESISTANCE column
Resistance_category <- Resistance_category %>%
  mutate(RESISTANCE = sapply(RESISTANCE, clean_string)) %>%
  distinct(name, .keep_all = TRUE) %>%
  column_to_rownames(var = "name") %>%
  select(RESISTANCE)

# Merge all the annotation rows
ann_row_all <- merge(ann_rows, Resistance_category, by = 0) %>%
  rename(DRUG = RESISTANCE) %>%
  column_to_rownames(var = "Row.names") %>%
  select(DRUG, AMR, VFs, MGEs)

# Safe annotation row to a csv file
# df <- write.csv(ann_row_all, "OneDrive - University of West London/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/Datasets/ann_row_all.csv", row.names = FALSE)

# Generate a color palette for drug categories
r_palette <- paletteer_d("ggsci::default_igv", n = length(unique(ann_row_all$DRUG))-1)

# Create a named dataframe for the Resistance category
ann_row_drug_col <- ann_row_all %>%
  select(DRUG) %>%
  na.omit() %>%
  distinct(DRUG) %>%
  mutate(color = r_palette)

# Create a named vector for rows
drug_colors <- setNames(as.character(ann_row_drug_col$color), ann_row_drug_col$DRUG)

# Annotation colors customization
ann_colors <- list(
  LOCATION = location_colors,
  DRUG = drug_colors,
  AMR = c(`Drug-resistant` = "#a02d02", `Multidrug-resistant` = "#522501"),
  MGEs = c(`Multi-MGEs` = "#5c0404", `Single MGE` = "#fcd2d2"),
  VFs = c(`Virulent` = "#fce5d2", `Hypervirulent` = "#c47a02")
)


# Create heatmap using pheatmap package ##
heatmap_plot <- pheatmap(sorted_matrix[1:70,], display_numbers = FALSE, cluster_cols = TRUE, cluster_rows = FALSE,
                         scale = "none", 
                         clustering_callback = callback,  
                         border_color = "NA", color = c("#CCCCCCFF", "#666666FF"),
                         legend_breaks = c(0, 1),
                         legend_labels = c("Absent", "Present"),
                         annotation_row = ann_row_all,
                         annotation_col = ann_col[4],
                         show_rownames = TRUE,
                         # cutree_cols = 15,
                         # cutree_rows = 20,
                         annotation_colors = ann_colors,
                         fontsize_row = 14,  # Adjust this value for row names
                         fontsize_col = 14,   # Adjust this value for column names
                         fontsize = 14  # Adjust this value for the legend text
                         
                         
)



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
