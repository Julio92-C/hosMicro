# Load packages
library(VennDiagram)
library(dplyr)
library(readr)
library(tidyr)
library(tidyverse)
library(stringr)




# Load dataset
# abri_kraken2_merged <- read_csv("C:/Users/peros/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/abri_kraken2_merged.csv")
# abri_kraken2_merged <- read_csv("C:/Users/peros/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/abricate_kraken2_merged.csv")
abri_kraken2_merged <- read_csv("OneDrive - University of West London/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/Datasets/abricate_kraken2_merged.csv")

# Vector of samples to remove
samples_to_remove <- c("A60B", "A61B", "A62B", "A63B", "A25R")


# Filter the dataset by database
abri_kraken2_filtered <- abri_kraken2_merged %>%
  filter(grepl("card|vfdb|plasmidfinder", DATABASE)) %>%
  filter(!sample %in% samples_to_remove) %>%
  mutate(sample = ifelse(grepl("A37R", sample), "A37", sample)) %>%
  arrange(name)

# Pivot the Sample values into columns
df_wide <- abri_kraken2_filtered %>%
  arrange(sample) %>%
  pivot_wider(names_from = sample, values_from = c(taxid)) %>%
  group_by(name) %>%
  summarise(across(16:75, ~ paste(., collapse = ", "))) %>%
  mutate(across(2:61, ~ str_replace_all(., "(NA,|,NA|,NA,|NA|,| )", "")))

colnames(df_wide)

# Assign the values of col1 as row names
# Air dataset
air_dataset <- df_wide %>%
  unite("Air", 2:28, sep = ",") %>%
  select(name, Air) %>%
  mutate(across(2, ~ str_replace_all(., "(,| )", "")))


# Surface dataset
Surface_dataset <- df_wide %>%
  unite("Surface", 29:61, sep = ",") %>%
  select(name, Surface) %>%
  mutate(across(2, ~ str_replace_all(., "(,| )", "")))

 # Merge Air and Surface dataset
air_Sur_dataset <- merge(air_dataset, Surface_dataset, by="name", all = TRUE)

# Assign the values of col1 as row names
air_Sur_dataset <- air_Sur_dataset %>%
  remove_rownames %>% 
  column_to_rownames(var="name") %>%
  mutate(across(everything(), ~ if_else(. == "", NA_character_, .)))


# Convert from a dataframe to a matrix
sample2arg_data_matrix <- data.matrix(air_Sur_dataset, rownames.force = NA)

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

# ?venn
# ?venn.diagram

### Plot the Venn diagram
venn.plot <- venn.diagram(
  x = sets,
  category.names = c("Air", "Surface"),
  fill = 1:2, 
  alpha = 0.3,
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



