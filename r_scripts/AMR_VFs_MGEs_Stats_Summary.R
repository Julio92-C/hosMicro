
# Load libraries
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(paletteer)

# Load dataset
abri_kraken2_merged <- read_csv("Hospital_microbiome/Datasets/abri_kraken2_merged1.csv")


# # Filter the most abundant species with counts greater than 10K
df_species <- abri_kraken2_filtered %>%
  group_by(name) %>%
  filter(!grepl("Terrabacteria group|Bacteroidota/Chlorobiota group|FCB group|Rhizobium/Agrobacterium group|Klebsiella/Raoultella group", name)) %>%
  filter(str_count(name, "\\S+") > 1) %>%
  mutate(
    name = str_replace(name, "^(\\w)\\w+\\s(\\w+)", "\\1. \\2"),
    name = if_else(
      str_length(name) > 30,
      str_replace(name, "^((?:\\S+\\s+){2}).*", "\\1"),
      name
    )
  ) %>%
  summarise(count = sum(count), .groups = "drop") %>%
  filter(count > 10000)


# # Species count
Species_count <- ggplot(df_species, aes(x = reorder(name, -count), y = count)) +
  geom_bar(stat = "identity", fill = "skyblue") +  # Use stat = "identity" for pre-calculated counts
  geom_text(aes(label = round(count/10000,1)), vjust = -1.5, size = 2.8) +  # Explicitly define y for labels
  labs(title = "",
       x = "Taxonomy Level",
       y = "log10(Taxa Count)") +
  scale_y_log10() +   #Change Y axis to logarithmic scale
  theme_classic() +  # Apply classic theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
        text = element_text(size = 14))  # Set text size

# Plot species count
Species_count

# Filter by Card DB and counts greater than 3
df_AMR <- abri_kraken2_filtered %>%
  filter(DATABASE == "card") %>%
  group_by(GENE) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  filter(count > 20) %>%
  arrange(desc(count)) # Arrange by count in descending order

df_AMR[14, 1] <- "CpxR"
df_AMR[17, 1] <- "Af_CAT"
# df_AMR[31, 1] <- "APH(2'')-Ia"
# df_AMR[33, 1] <- "mdfA"


# # Plot AMR card database
AMR_count <- ggplot(df_AMR, aes(x = reorder(GENE, -count), y = count)) +
  geom_bar(stat = "identity", fill = "#b02d02") +  # Change the color of the bars
  geom_text(aes(label = round(count/1,1)), vjust = -0.5, size = 2.8) +  # Add count labels on top of the bars
  labs(title = "",
       x = "Antibiotics Resistance Genes (ARGs)",
       y = "ARG Count") +
  theme_classic() +  # Apply classic theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  theme(text = element_text(size = 16))

# Plot AMR count 
AMR_count


# Define the function to extract the Virulence factor "Functions"
extract_productFunction <- function(string) {
  match <- regmatches(string, regexpr("\\- ([^\\(]+) \\(", string))
  result <- gsub("\\- | \\(", "", match)
  return(result)
}

# Apply the function to the dataframe
refactored_df <- abri_kraken2_merged %>%
  filter(DATABASE == "vfdb") %>%
  mutate(Functions = sapply(PRODUCT, extract_productFunction))

# Filter by vfdb and counts greater than 3
vfdb_filtered <- refactored_df %>%
  #filter(DATABASE == "vfdb") %>%
  group_by(GENE, Functions, type) %>%
  #filter(!grepl("Exotoxin|Immune modulation|Antimicrobial activity/Competitive advantage|Effector delivery system", Functions)) %>%
  filter(grepl("Effector delivery system|Adherence|Regulation|Nutritional/Metabolic factor|Invasion", Functions)) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  filter(count > 10) %>%
  arrange(desc(count)) # Arrange by count in descending order  


# # Plot VFs vfdb database
VFs_count <- ggplot(vfdb_filtered, aes(x = reorder(GENE, -count), y = count)) +
  geom_bar(stat = "identity", aes(fill = type)) +  # Change the color of the bars
  geom_text(aes(label = round(count/1,1)), vjust = -0.5, size = 2.8) +  # Add count labels on top of the bars
  scale_fill_manual(values = c(Air = "#828385", Surface = "#ed8cb6")) +
  labs(title = "",
       x = "Virulence Factors (VFs)",
       y = "VFs Count",
       fill = "Sample") +
  theme_classic() +  # Apply classic theme
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  theme(text = element_text(size = 14)) +
  facet_wrap(~ Functions, scales = "free_x", nrow = 3)


# # Plot VFs
VFs_count

# Filter by plasmidfinder and counts greater than 3
plasmid_filtered <- abri_kraken2_filtered %>%
  filter(DATABASE == "plasmidfinder") %>%
  group_by(GENE) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  filter(count > 40) %>%
  arrange(desc(count)) # Arrange by count in descending order  

# # Plot MGEs counts
MGEs_count <- ggplot(plasmid_filtered, aes(x = reorder(GENE, -count), y = count)) +
  geom_bar(stat = "identity", fill = "#5c0404") +  # Change the color of the bars
  geom_text(aes(label = round(count/1,1)), vjust = -0.5, size = 2.8) +  # Add count labels on top of the bars
  #scale_fill_manual(values = t_palette) +
  labs(title = "",
       x = "Treatment groups",
       y = "MGE Count") +
  #scale_y_log10() +  # Change Y axis to logarithmic scale
  theme_classic() +  # Apply classic theme
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  theme(text = element_text(size = 15))

# Plot MGEs count
MGEs_count

# Count occurrences
data_summary <- abri_kraken2_filtered %>%
  group_by(sample, DATABASE) %>%
  summarise(Count = n())

# Calculate percentages
data_summary <- data_summary %>%
  group_by(sample) %>%
  mutate(Percentage = Count / sum(Count) * 100)



#  stacked barplot where each bar represents a sample, 
#  and the segments within each bar represent the counts of different genes
ggplot(data_summary, aes(x = factor(sample), y = Percentage, fill = DATABASE)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#b02d02", "#5c0404", "#c47a02")) +
  labs(x = "Sample", y = "Percentage", fill = "DATABASE") +
  theme_classic() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Filter the most abundant species with counts greater than 15K
df_species <- abri_kraken2_merged %>%
  group_by(name) %>%
  filter(count > 15000) %>%
  group_by(sample) %>%
  mutate(Percentage = count / sum(count) * 100)

# Annotation col colors for taxa
# Description variable
taxa <- as.factor(df_species$name)

# Generate a color palette based on the number of levels in gene_family
t <- length(levels(taxa))
t_palette <- paletteer_d("ggsci::default_igv", n = t)
t_palette


#  stacked barplot where each bar represents a sample, 
#  and the segments within each bar represent the relative abundance of different taxa
ggplot(df_species, aes(x = factor(sample), y = Percentage, fill = name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = t_palette) +
  labs(x = "Air and Surface sample groups", y = "Percentage", fill = "name") +
  theme_classic() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 14)) +
  facet_wrap(~ type, scales = "free_x", nrow = 1)


# Genetic elements distribution 
Gene_count <- table(abri_kraken2_filtered$DATABASE)

# PLot the Gene count by database 
pie_chart <- pie(Gene_count, labels = paste0(names(Gene_count), " / ", Gene_count, " / ", round(100 * Gene_count/sum(Gene_count), 2), "%"),
                 main="Gene count by database",
                 col = c("#b02d02", "#5c0404", "#c47a02"),
                 cex = 0.8, # Adjust label size
                 radius = 1, # adjust the size
                 xlim = c(-1.5, 1.5) # Adjust x-axis limits to create more space for labels
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