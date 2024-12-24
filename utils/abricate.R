
# Load libraries
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)


# Import dataset
Summary <- read_csv("C:/Users/peros/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/Summary.csv")
View(Summary)

# Tax ID
# kraken2_db1_combined_reports <- read_delim("C:/Users/peros/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/Kraken2/kraken2_db1_combined_reports.txt", 
#                                            delim = "\t", escape_double = FALSE, 
#                                            trim_ws = TRUE)
# colnames(kraken2_db1_combined_reports)

# Surface Kraken report
H1S_kraken2_db1_combined_reports <- read_delim("C:/Users/peros/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/H1S_kraken2_db1_combined_reports.txt", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE)

# Air Kraken Report
H1A_kraken2_db1_combined_reports <- read_delim("C:/Users/peros/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/H1A_kraken2_db1_combined_reports.txt", 
                                               delim = "\t", escape_double = FALSE, 
                                               trim_ws = TRUE)



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
taxid_name_Sur <- H1S_kraken2_db1_combined_reports[, c("taxid", "name")]
taxid_name_Air <- H1A_kraken2_db1_combined_reports[, c("taxid", "name")]

# Merge Air and Sur tax id
taxid_name <- merge(taxid_name_Air, taxid_name_Sur, by=c("taxid","name"))


# Merge abricate and kraken2 output
abri_kraken2_merged <- merge(Summary, taxid_name, by ="taxid") 

# Clean up the abri_kraken2_merged dataframe
abri_kraken2_merged <- abri_kraken2_merged %>%
  filter(!grepl("root|Homo sapiens|cellular organisms", name)) %>%
  filter(name != "Bacteria") 


# Save abri_kraken2_merged dadaset as csv file
df <- write.csv(abri_kraken2_merged, "C:/Users/peros/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/abricate_kraken2_merged.csv", row.names = FALSE)

# Filter by Species and counts greater than 3
df_species <- abri_kraken2_merged %>%
  filter(!grepl("root|Homo sapiens|cellular organisms", name)) %>%
  filter(name != "Bacteria") %>%
  #filter(DATABASE == "card") %>%
  #distinct(name, GENE) %>%
  group_by(name) %>%
  arrange(name) %>%
  filter(n() > 30)

df_species[631:645, 17] <- "Acinetobacter radioresistens"

# Species count
Species_count <- ggplot(df_species, aes(x = reorder(name, -table(name)[name]))) +
  geom_bar(fill = "Skyblue") +  # Change the color of the bars
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +  # Add count labels on top of the bars
  labs(title = "",
       x = "Species",
       y = "Count") +
  theme_classic() +  # Apply classic theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))   # Rotate x-axis labels

# Plot species count
Species_count



# Filter by Card DB and counts greater than 3
df_filtered <- abri_kraken2_merged %>%
  filter(!grepl("root|Homo sapiens|cellular organisms", name)) %>%
  filter(name != "Bacteria") %>%
  filter(DATABASE == "card") %>%
  #distinct(name, GENE) %>%
  group_by(GENE) %>%
  arrange(GENE) %>%
  filter(n() > 10)


df_filtered[1712:1747, 7] <- "CpxR"
df_filtered[136:160, 7] <- "Af_CAT"
 
# # Plot AMR card database
AMR_count <- ggplot(df_filtered, aes(x = reorder(GENE, -table(GENE)[GENE]))) +
  geom_bar(fill = "#b02d02") +  # Change the color of the bars
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +  # Add count labels on top of the bars
  labs(title = "",
       x = "AMR",
       y = "Count") +
  theme_classic() +  # Apply classic theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))   # Rotate x-axis labels

# Plot AMR count 
AMR_count

# Filter by vfdb and counts greater than 3
vfdb_filtered <- abri_kraken2_merged %>%
  filter(!grepl("root|Homo sapiens|cellular organisms", name)) %>%
  filter(name != "Bacteria") %>%
  filter(DATABASE == "vfdb") %>%
  #distinct(name, GENE) %>%
  group_by(GENE) %>%
  arrange(GENE) %>%
  filter(n() > 20)  

# # Plot VFs vfdb database
VFs_count <- ggplot(vfdb_filtered, aes(x = reorder(GENE, -table(GENE)[GENE]))) +
  geom_bar(fill = "#b02d02") +  # Change the color of the bars
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +  # Add count labels on top of the bars
  labs(title = "",
       x = "VFs",
       y = "Count") +
  theme_classic() +  # Apply classic theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))   # Rotate x-axis labels


# # Plot VFs
VFs_count

# Filter by plasmidfinder and counts greater than 3
plasmid_filtered <- abri_kraken2_merged %>%
  filter(!grepl("root|Homo sapiens|cellular organisms", name)) %>%
  filter(name != "Bacteria") %>%
  filter(DATABASE == "plasmidfinder") %>%
  #distinct(name, GENE) %>%
  group_by(GENE) %>%
  arrange(GENE) %>%
  filter(n() > 1)  

plasmid_filtered[1186:1560, 7] <- "rep14a_4"

# # Plot MGEs counts
MGEs_count <- ggplot(plasmid_filtered, aes(x = reorder(GENE, -table(GENE)[GENE]))) +
  geom_bar(fill = "#c47a02") +  # Change the color of the bars
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +  # Add count labels on top of the bars
  labs(title = "",
       x = "MGEs",
       y = "Count") +
  theme_classic() +  # Apply classic theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))   # Rotate x-axis labels

# Plot MGEs count
MGEs_count

# Arrange plots in a 2x2 grid
combined_plot <- ggarrange(AMR_count, VFs_count, Species_count, MGEs_count, 
                           ncol = 2, nrow = 2)
combined_plot


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

