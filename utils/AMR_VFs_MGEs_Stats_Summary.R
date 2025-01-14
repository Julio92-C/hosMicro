
# Load libraries
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(paletteer)

# Load dataset
# abri_kraken2_merged <- read_csv("C:/Users/peros/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/abricate_kraken2_merged.csv")
# abri_kraken2_merged <- read_csv("OneDrive - University of West London/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/Datasets/abricate_kraken2_merged.csv")

## Load clean dataset
abri_kraken2_filtered <- read_csv("OneDrive - University of West London/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/Datasets/abri_kraken2_cleaned.csv")


sample <- length(unique(abri_kraken2_filtered$sample))
taxa <- length(unique(abri_kraken2_filtered$name))
gene <- length(unique(abri_kraken2_filtered$GENE))
resistance <- length(unique(abri_kraken2_filtered$RESISTANCE))

# Write abri_kraken2 to a csv file
# write.csv(abri_kraken2_filtered, "OneDrive - University of West London/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/Datasets/abricate_kraken2_filtered.csv", row.names = FALSE)

# pseudomonas <- filter(abri_kraken2_filtered, abri_kraken2_filtered$name == "Pseudomonas")
# Stutzerimonas_stutzeri <- filter(abri_kraken2_filtered, abri_kraken2_filtered$name == "Stutzerimonas stutzeri")
# Enterococcus_faecium <- filter(abri_kraken2_filtered, abri_kraken2_filtered$name == "Enterococcus faecium")
# Acinetobacter_lwoffii <- filter(abri_kraken2_filtered, abri_kraken2_filtered$name == "Acinetobacter lwoffii")
# Enterobacterales <- filter(abri_kraken2_filtered, abri_kraken2_filtered$name == "Enterobacterales")
# Enterobacterale_Air <- Enterobacterales %>%
#   filter(grepl("A", sample)) 
# 
# Enterobacterale_Air_sample <- as.factor(Enterobacterale_Air$sample)
# 
# Bacilli <- filter(abri_kraken2_filtered, abri_kraken2_filtered$name == "Bacilli")
# Bacilli_Air <- Bacilli %>%
#   filter(grepl("A", sample)) 
# 
# Bacilli_Air_sample <- as.factor(Bacilli_Air$sample)
# 
# Rahnella <- filter(abri_kraken2_filtered, abri_kraken2_filtered$name == "Rahnella")
# Rahnella_Air <- Rahnella %>%
#   filter(grepl("A", sample)) 
# 
# Rahnella_Air_sample <- as.factor(Rahnella_Air$sample)
# 
# ARGs_type <- filter(abri_kraken2_filtered, abri_kraken2_filtered$DATABASE == "card")
# 
# ARGs_type_category <- as.factor(unique(ARGs_type$GENE))
# 
# VFs_type <- filter(abri_kraken2_filtered, abri_kraken2_filtered$DATABASE == "vfdb")
# 
# VFs_type_category <- as.factor(unique(VFs_type$GENE))
# 
# MGEs_type <- filter(abri_kraken2_filtered, abri_kraken2_filtered$DATABASE == "plasmidfinder")
# 
# MGEs_type_category <- as.factor(unique(MGEs_type$GENE))
# 
# # Filter the dataset by T6SS gene family
# T6SS_family <- abri_kraken2_filtered %>%
#   filter(grepl("tss", GENE))
# 
# T6SS_family_category <- as.factor(unique(T6SS_family$GENE))
#   
# 
# 
# # Unique species
# unique_species <- abri_kraken2_filtered %>%
#    distinct(name)

# Filter by Species and counts greater than 30
df_species <- abri_kraken2_filtered %>%
  group_by(name) %>%
  arrange(name) %>%
  filter(n() > 50)

df_species[631:645, 17] <- "Acinetobacter radioresistens"
df_species[7021:7052, 17] <- "Rahnella aquatilis"

# Species count
Species_count <- ggplot(df_species, aes(x = reorder(name, -table(name)[name]))) +
  geom_bar(fill = "skyblue") +  # Change the color of the bars
  #scale_fill_manual(values = t_palette) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +  # Add count labels on top of the bars
  labs(title = "",
       x = "Taxonomy Level",
       y = "Count") +
  theme_classic() +  # Apply classic theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  theme(text = element_text(size = 16))
  
# Plot species count
Species_count



# Filter by Card DB and counts greater than 3
df_AMR <- abri_kraken2_filtered %>%
  filter(DATABASE == "card") %>%
  group_by(GENE) %>%
  arrange(GENE) %>%
  filter(n() > 20)


df_AMR[2536:2607, 7] <- "CpxR"
df_AMR[273:322, 7] <- "Af_CAT"

# # Plot AMR card database
AMR_count <- ggplot(df_AMR, aes(x = reorder(GENE, -table(GENE)[GENE]))) +
  geom_bar(fill = "#b02d02") +  # Change the color of the bars
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +  # Add count labels on top of the bars
  labs(title = "",
       x = "ARGs",
       y = "Count") +
  theme_classic() +  # Apply classic theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  theme(text = element_text(size = 16))
  
# Plot AMR count 
AMR_count

# Filter by vfdb and counts greater than 3
vfdb_filtered <- abri_kraken2_filtered %>%
  filter(DATABASE == "vfdb") %>%
  group_by(GENE) %>%
  arrange(GENE) %>%
  filter(n() > 30)  

# # Plot VFs vfdb database
VFs_count <- ggplot(vfdb_filtered, aes(x = reorder(GENE, -table(GENE)[GENE]))) +
  geom_bar(fill = "#c47a02") +  # Change the color of the bars
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +  # Add count labels on top of the bars
  labs(title = "",
       x = "VFs",
       y = "Count") +
  theme_classic() +  # Apply classic theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  + # Rotate x-axis labels
  theme(text = element_text(size = 16))

# # Plot VFs
VFs_count

# Filter by plasmidfinder and counts greater than 3
plasmid_filtered <- abri_kraken2_filtered %>%
  filter(DATABASE == "plasmidfinder") %>%
  group_by(GENE) %>%
  arrange(GENE) %>%
  filter(n() > 9)  

#plasmid_filtered[1186:1560, 7] <- "rep14a_4"

# # Plot MGEs counts
MGEs_count <- ggplot(plasmid_filtered, aes(x = reorder(GENE, -table(GENE)[GENE]))) +
  geom_bar(fill = "#5c0404") +  # Change the color of the bars
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +  # Add count labels on top of the bars
  labs(title = "",
       x = "MGEs",
       y = "Count") +
  theme_classic() +  # Apply classic theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  theme(text = element_text(size = 14))

# Plot MGEs count
MGEs_count

# Arrange plots in a 2x2 grid
combined_plot <- ggarrange(AMR_count, VFs_count, MGEs_count, pie_chart,
                           ncol = 2, nrow = 2)
combined_plot


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
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 14))

# Filter by Species and counts greater than 50
df_relativeAbundance <- df_species %>%
  group_by(sample, name) %>%
  summarise(Count = n()) %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Annotation col colors for taxa
# Description variable
taxa <- as.factor(df_species$name)

# Generate a color palette based on the number of levels in gene_family
t <- length(levels(taxa))
t_palette <- paletteer_d("ggsci::default_igv", n = t)
t_palette



#  stacked barplot where each bar represents a sample, 
#  and the segments within each bar represent the relative abundance of different taxa
ggplot(df_relativeAbundance, aes(x = factor(sample), y = Percentage, fill = name)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = t_palette) +
  labs(x = "Air and Surface sample groups", y = "Percentage", fill = "name") +
  theme_classic() +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(text = element_text(size = 16))


# Determine gene count in filter dataset
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