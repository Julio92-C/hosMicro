# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyverse)
library(stringr)
library(ggsignif)  # For adding statistical significance



# Load abri_kraken2_merged dataset
# abri_kraken2_merged <- read_csv("C:/Users/peros/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/abricate_kraken2_merged.csv")
abri_kraken2_merged <- read_csv("OneDrive - University of West London/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/Datasets/abricate_kraken2_merged.csv")
MetadataLocations <- read_csv("OneDrive - University of West London/Desktop/PhD Proposal/Bioinformatics/Oscar_metagenomics/Datasets/MetadataLocations1.csv")



# Vector of samples to remove
samples_to_remove <- c("A60B", "A61B", "A62B", "A63B", "A25R")


# Filter the data set by AMR, VFs and MGEs
abri_kraken2_filtered <- abri_kraken2_merged %>%
  filter(grepl("card|vfdb|plasmidfinder", DATABASE)) %>%
  arrange(name) %>%
  filter(!sample %in% samples_to_remove) %>%
  mutate(sample = ifelse(grepl("A37R", sample), "A37", sample)) %>%
  group_by(sample, name) %>%
  summarise(Taxa_count = n(), .groups = 'drop') %>%
  ungroup() %>% 
  mutate(Taxa_log_count = log(Taxa_count + 1))  # Adding 1 to avoid log(0)


# Rename samples
abri_kraken2_filtered[1:536, 1] <- "Air"
abri_kraken2_filtered[537:1669, 1] <- "Surface"

# Calculate p-value for significance
p_value <- t.test(Taxa_log_count ~ sample, data = abri_kraken2_filtered)$p.value

# Generate the violin plot with statistical significance for taxa
ggplot(abri_kraken2_filtered, aes(x = sample, y = Taxa_log_count, fill = sample)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  labs(title = "",
       x = "Sample groups",
       y = "Number of Taxa Log") +
  theme_minimal() +
  theme(legend.position = "top") +
  theme(text = element_text(size = 16)) +
  geom_signif(comparisons = list(c("Air", "Surface")), 
              map_signif_level = TRUE, 
              textsize = 3.5) +
  annotate("text", x = 1.5, y = max(abri_kraken2_filtered$Taxa_log_count) + 1, 
           label = paste("p =", format(p_value, digits = 2)), 
           size = 5, color = "black")

# Filter the data set by AMR
abri_kraken2_filtered <- abri_kraken2_merged %>%
  filter(grepl("card", DATABASE)) %>%
  arrange(name) %>%
  group_by(sample, GENE) %>%
  summarise(Gene_count = n(), .groups = 'drop') %>%
  ungroup() %>% 
  mutate(Gene_log_count = log(Gene_count + 1))  # Adding 1 to avoid log(0)

# Rename samples
abri_kraken2_filtered[1:93, 1] <- "Air"
abri_kraken2_filtered[94:342, 1] <- "Surface"


# Calculate p-value for significance
p_value <- t.test(Gene_log_count ~ sample, data = abri_kraken2_filtered)$p.value

t_test_result <- t.test(Gene_log_count ~ sample, data = abri_kraken2_filtered)
print(t_test_result)


# Generate the violin plot with statistical significance for AMR
ggplot(abri_kraken2_filtered, aes(x = sample, y = Gene_log_count, fill = sample)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.6) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#04b5d9", "#aeb6b8")) +
  labs(title = "",
       x = "Sample groups",
       y = "Log(Number of ARGs)") +
  theme_classic() +
  theme(legend.position = "top") +
  theme(text = element_text(size = 16)) +
  geom_signif(comparisons = list(c("Air", "Surface")), 
              map_signif_level = TRUE, 
              textsize = 3.5) +
  annotate("text", x = 1.5, y = max(abri_kraken2_filtered$Gene_log_count) + 1.5, 
           label = paste("p =", format(p_value, digits = 2)), 
           size = 4, color = "black")



#### Filter the data set by VFs #####
abri_kraken2_filtered <- abri_kraken2_merged %>%
  filter(grepl("vfdb", DATABASE)) %>%
  arrange(name) %>%
  group_by(sample, GENE) %>%
  summarise(Gene_count = n(), .groups = 'drop') %>%
  ungroup() %>% 
  mutate(Gene_log_count = log(Gene_count + 1))  # Adding 1 to avoid log(0)

# Rename samples
abri_kraken2_filtered[1:81, 1] <- "Air"
abri_kraken2_filtered[82:742, 1] <- "Surface"


# Calculate p-value for significance
p_value <- t.test(Gene_log_count ~ sample, data = abri_kraken2_filtered)$p.value

t_test_result <- t.test(Gene_log_count ~ sample, data = abri_kraken2_filtered)
print(t_test_result)


# Generate the violin plot with statistical significance for AMR
ggplot(abri_kraken2_filtered, aes(x = sample, y = Gene_log_count, fill = sample)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.6) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#04b5d9", "#aeb6b8")) +
  labs(title = "",
       x = "Sample",
       y = "Log(Number of VFs)") +
  theme_classic() +
  theme(legend.position = "top") +
  geom_signif(comparisons = list(c("Air", "Surface")), 
              map_signif_level = TRUE, 
              textsize = 3.5) +
  annotate("text", x = 1.5, y = max(abri_kraken2_filtered$Gene_log_count) + 1.5, 
           label = paste("p =", format(p_value, digits = 2)), 
           size = 4, color = "black")

#### Filter the data set by MGEs #####
abri_kraken2_filtered <- abri_kraken2_merged %>%
  filter(grepl("plasmidfinder", DATABASE)) %>%
  arrange(name) %>%
  group_by(sample, GENE) %>%
  summarise(Gene_count = n(), .groups = 'drop') %>%
  ungroup() %>% 
  mutate(Gene_log_count = log(Gene_count + 1))  # Adding 1 to avoid log(0)

# Rename samples
abri_kraken2_filtered[1:66, 1] <- "Air"
abri_kraken2_filtered[67:107, 1] <- "Surface"


# Calculate p-value for significance
p_value <- t.test(Gene_log_count ~ sample, data = abri_kraken2_filtered)$p.value

t_test_result <- t.test(Gene_log_count ~ sample, data = abri_kraken2_filtered)
print(t_test_result)


# Generate the violin plot with statistical significance for AMR
ggplot(abri_kraken2_filtered, aes(x = sample, y = Gene_log_count, fill = sample)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.6) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#04b5d9", "#aeb6b8")) +
  labs(title = "",
       x = "Sample",
       y = "Log(Number of MGEs)") +
  theme_classic() +
  theme(legend.position = "top") +
  geom_signif(comparisons = list(c("Air", "Surface")), 
              map_signif_level = TRUE, 
              textsize = 3.5) +
  annotate("text", x = 1.5, y = max(abri_kraken2_filtered$Gene_log_count) + 1.5, 
           label = paste("p =", format(p_value, digits = 2)), 
           size = 4, color = "black")

# CLEAN UP #################################################

# Clear environment
rm(list = ls()) 

# Clear plots
if (length(dev.list()) > 0) {
  dev.off()  # Only if there is an active plot
}

# Clear console
cat("\014")  # ctrl+L





