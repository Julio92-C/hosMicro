
# Load libraries
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)


# Import Abricate summary report
Summary <- read_csv("Hospital_microbiome/Datasets/Summary.csv")
View(Summary)


# Surface Kraken report
H1S_kraken2_db1_combined_reports <- read_delim("Hospital_microbiome/Datasets/H1S_kraken2_db1_combined_reports.txt", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE)

# Air Kraken Report
H1A_kraken2_db1_combined_reports <- read_delim("Hospital_microbiome/Datasets/H1A_kraken2_db1_combined_reports.txt", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE)

# Recentrifuge contaminats list
contaminats <- read_csv("Hospital_microbiome/Datasets/NX_AH4_contam.csv")


# Rename the column to "sample"
names(Summary)[names(Summary) == "#FILE"] <- "sample"


# Remove the ".fasta" extension
Summary$sample <- gsub(".fasta", "", Summary$sample)

# Step 1: Split the string into two parts at the underscore
Summary <- Summary %>%
  separate(SEQUENCE, into = c("sequence", "taxid"), sep = "_")

# Step 2: Extract the part after the '|' and create a new column called 'taxid'
Summary <- Summary %>%
  mutate(taxid = sub(".*\\|", "", taxid))

# Subset taxID and and name from kraken2 output
# taxid_name <- kraken2_db1_combined_reports[, c("taxid", "name")]
taxid_name_Sur <- H1S_kraken2_db1_combined_reports[, c("taxid", "name", "#perc", "tot_all")]
taxid_name_Air <- H1A_kraken2_db1_combined_reports[, c("taxid", "name", "#perc", "tot_all")]


# Remove contaminants from a subset of the air samples (NX_AH4_contamint.csv)
# Merge the taxid_name_Air and contaminats datasets
contaminats_taxid_name_Air <- merge(contaminats, taxid_name_Air, by="taxid", all.x = TRUE)

#  Pivoting samples and control into a single column
contaminants_pivoted <- contaminats_taxid_name_Air %>%
  pivot_longer(cols = A45:A57, names_to = "sample", values_to = "count") %>%
  mutate(NX_difference = NX - count) %>%
  filter(NX_difference <= 0)

# Merge Air and Sur tax id
taxid_name <- bind_rows(taxid_name_Air, taxid_name_Sur)


# Merge abricate and kraken2 output
abri_kraken2_merged <- merge(Summary, taxid_name, by ="taxid", all.x = TRUE)


# Vector of samples to remove
samples_to_remove <- c("A60B", "A61B", "A62B", "A63B", "A25R")

# Perform the anti-join to remove contaminants:
abri_kraken2_cleaned <- abri_kraken2_merged %>%
  mutate(taxid = as.double(taxid)) %>%
  anti_join(contaminants_pivoted, by = c("taxid", "sample")) %>%
  filter(!grepl("root|Homo sapiens|cellular organisms|Bacteria", name)) %>%
  filter(!sample %in% samples_to_remove) %>%
  mutate(sample = ifelse(grepl("A37R", sample), "A37", sample)) %>%
  filter(grepl("card|vfdb|plasmidfinder", DATABASE)) %>%
  rename(Percentage = "#perc") %>%
  rename(count = "tot_all") %>%
  arrange(name)


# Save abri_kraken2_merged dadaset as csv file
write.csv(abri_kraken2_cleaned, "Hospital_microbiome/Datasets/abri_kraken2_cleaned1.csv", row.names = FALSE)


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

