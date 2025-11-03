# Load the libraries
library(readr)
library(dplyr)
library(tidyverse)
library(paletteer)
library(reshape2)
library(stringr)
library(VennDiagram)




# Load dataset
abri_kraken2_merged <- read_csv("Hospital_microbiome/Datasets/abri_kraken2_merged1.csv")
# colnames(abri_kraken2_merged)


# Annotation rows
# Pivot the Sample values into columns
df_wide_ann_rows <- abri_kraken2_merged %>%
  arrange(sample) %>%
  pivot_wider(names_from = DATABASE, values_from = c(GENE), values_fn = list) %>%
  group_by(name) %>%
  summarise(across(22:24, ~ paste(., collapse = ", "))) %>%
  mutate(across(2:4, ~ str_replace_all(., "(NULL,|,NULL|,NULL,|NULL| )", "")))

# Assign the values of Sample as row names
df_wide_ann_rows <- df_wide_ann_rows %>%
  na.omit()  %>%
  remove_rownames %>%
  column_to_rownames(var="name") %>%
  mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))


# Convert from a dataframe to a matrix
sample2arg_data_matrix <- data.matrix(df_wide_ann_rows, rownames.force = NA)

# Count NA in data sets
sum(is.na(sample2arg_data_matrix))

# Replace all NA values with 0 values
sample2arg_data_matrix <- sample2arg_data_matrix %>% replace(is.na(.), 0)

# Convert to a binary matrix
sample2arg_data_Bmatrix <- as.matrix((sample2arg_data_matrix > 0) + 0)

# Create a list of sets
sets <- apply(sample2arg_data_Bmatrix, 2, function(col) which(col == 1))
names(sets) <- colnames(sample2arg_data_Bmatrix)
print(sets)
class(sets)

### Plot the Venn diagram
venn.plot <- venn.diagram(
  x = sets,
  category.names = c("ARGs", "VFs", "MGEs"),
  fill = c(card="#b02d02", vfdb="#c47a02", plasmidfinder="#5c0404"), 
  alpha = 0.5,
  height = 50, width = 50,
  filename = NULL,
  cat.cex = 1.5,  # Adjust this value to change the font size
  cex = 2       # Adjust this value to change the font size of the numbers
)
grid.newpage()
grid.draw(venn.plot)



# Clear environment
rm(list = ls())

# Clear packages
detach("package:datasets", unload = TRUE)

# Clear plots
dev.off()  # But only if there IS a plot

# Clear console
cat("\014")  # ctrl+L

# Clear mind :)
