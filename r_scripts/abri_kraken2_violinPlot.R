# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)
library(tidyverse)
library(stringr)
library(ggsignif)  # For adding statistical significance



# Load abri_kraken2_merged dataset
abri_kraken2_merged <- read_csv("Hospital_microbiome/Datasets/abri_kraken2_merged1.csv")
colnames(abri_kraken2_merged)

# Filter the data set by AMR, VFs and MGEs
abri_kraken2_filtered <- abri_kraken2_merged %>%
  filter(grepl("card|vfdb|plasmidfinder", DATABASE)) %>%
  arrange(name) %>%
  group_by(sample, type, name) %>%
  summarise(Taxa_count = n(), .groups = 'drop') %>%
  ungroup() %>% 
  mutate(Taxa_log_count = log(Taxa_count + 1))  # Adding 1 to avoid log(0)


# Calculate p-value for significance
p_value <- t.test(Taxa_log_count ~ type, data = abri_kraken2_filtered)$p.value

# Generate the violin plot with statistical significance for taxa
ggplot(abri_kraken2_filtered, aes(x = type, y = Taxa_log_count, fill = type)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  #geom_jitter(width = 0.1, size = 2, alpha = 0.8) + # Jittered points for individual samples
  scale_fill_manual(values = c(Air = "#828385", Surface = "#ed8cb6")) +
  labs(title = "",
       x = "Sample groups",
       y = "Log(Number of Taxa)",
       fill = "Sample") +
  theme_classic() +
  theme(legend.position = "top") +
  theme(text = element_text(size = 14)) +
  geom_signif(comparisons = list(c("Air", "Surface")),
              map_signif_level = TRUE,
              textsize = 3.5) +
  annotate("text", x = 1.5, y = max(abri_kraken2_filtered$Taxa_log_count) + 0.5,
           label = paste("p =", format(p_value, digits = 2)),
           size = 5, color = "black")

# Filter the data set by AMR
abri_kraken2_filtered <- abri_kraken2_merged %>%
  filter(grepl("card", DATABASE)) %>%
  arrange(name) %>%
  group_by(sample, GENE, type) %>%
  summarise(Gene_count = n(), .groups = 'drop') %>%
  ungroup() %>% 
  mutate(Gene_log_count = log(Gene_count + 1))  # Adding 1 to avoid log(0)

# Calculate p-value for significance
p_value <- t.test(Gene_log_count ~ type, data = abri_kraken2_filtered)$p.value

t_test_result <- t.test(Gene_log_count ~ type, data = abri_kraken2_filtered)
print(t_test_result)


# Generate the violin plot with statistical significance for AMR
ggplot(abri_kraken2_filtered, aes(x = type, y = Gene_log_count, fill = type)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.6) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  #geom_jitter(width = 0.1, size = 2, alpha = 0.8) + # Jittered points for individual samples
  scale_fill_manual(values = c(Air = "#828385", Surface = "#ed8cb6")) +
  labs(title = "",
       x = "Sample groups",
       y = "Log(Number of ARGs)",
       fill = "Sample") +
  theme_classic() +
  theme(legend.position = "top") +
  theme(text = element_text(size = 14)) +
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
  group_by(sample, GENE, type) %>%
  summarise(Gene_count = n(), .groups = 'drop') %>%
  ungroup() %>% 
  mutate(Gene_log_count = log(Gene_count + 1))  # Adding 1 to avoid log(0)

# Calculate p-value for significance
p_value <- t.test(Gene_log_count ~ type, data = abri_kraken2_filtered)$p.value

t_test_result <- t.test(Gene_log_count ~ type, data = abri_kraken2_filtered)
print(t_test_result)


# Generate the violin plot with statistical significance for VFs
ggplot(abri_kraken2_filtered, aes(x = type, y = Gene_log_count, fill = type)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.6) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  #geom_jitter(width = 0.1, size = 2, alpha = 0.8) + # Jittered points for individual samples
  scale_fill_manual(values = c(Air = "#828385", Surface = "#ed8cb6")) +
  labs(title = "",
       x = "Sample groups",
       y = "Log(Number of VFs)",
       fill = "Sample") +
  theme_classic() +
  theme(legend.position = "top") +
  theme(text = element_text(size = 14)) +
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
  group_by(sample, GENE, type) %>%
  summarise(Gene_count = n(), .groups = 'drop') %>%
  ungroup() %>% 
  mutate(Gene_log_count = log(Gene_count + 1))  # Adding 1 to avoid log(0)


# Calculate p-value for significance
p_value <- t.test(Gene_log_count ~ type, data = abri_kraken2_filtered)$p.value

t_test_result <- t.test(Gene_log_count ~ type, data = abri_kraken2_filtered)
print(t_test_result)


# Generate the violin plot with statistical significance for MGEs
ggplot(abri_kraken2_filtered, aes(x = type, y = Gene_log_count, fill = type)) +
  geom_violin(trim = FALSE, scale = "width", width = 0.6) +
  geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
  # geom_jitter(width = 0.1, size = 2, alpha = 0.8) + # Jittered points for individual samples
  scale_fill_manual(values = c(Air = "#828385", Surface = "#ed8cb6")) +
  labs(title = "",
       x = "Sample groups",
       y = "Log(Number of MGEs)",
       fill = "Sample") +
  theme_classic() +
  theme(text = element_text(size = 14)) +
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





