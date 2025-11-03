# Load libraries
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)
library(paletteer)
library(gtsummary)
library(gt)
library(glue)



# Load dataset
abri_kraken2_merged <- read_csv("hospital_Microbiome/Datasets/abri_kraken2_merged1.csv")


# Perform one-way ANOVA to assess statistical differences between treatment groups
#### Filter the data set by MGEs #####
abri_kraken2_filtered <- abri_kraken2_merged %>%
  #filter(grepl("plasmidfinder", DATABASE)) %>%
  #filter(grepl("card", DATABASE)) %>%
  filter(grepl("vfdb", DATABASE)) %>%
  arrange(name) %>% 
   group_by(type, GENE) %>%
  summarise(Gene_count = n(), .groups = 'drop') %>%
  ungroup() %>% 
  mutate(Gene_log_count = log(Gene_count + 1))  # Adding 1 to avoid log(0)

# one-way ANOVA analysis
anova_result <- aov(Gene_log_count ~ type, data = abri_kraken2_filtered)
summary(anova_result)


# Extract and format the p-value
anova_p <- format(summary(anova_result)[[1]][["Pr(>F)"]][1], digits = 2)

# Arrange the dataset before plotting the table
abri_kraken2_arranged <- abri_kraken2_merged %>%
  #filter(grepl("plasmidfinder", DATABASE)) %>%
  #filter(grepl("card", DATABASE)) %>%
  filter(grepl("vfdb", DATABASE)) %>%
  group_by(type, GENE, name) %>%
  summarise(Gene_count = n(), .groups = 'drop') %>%
  ungroup() %>% 
  #filter(Gene_count > 10)
  # filter(Gene_count > 25)
  filter(Gene_count > 50)

# Summarise gene counts per treatment
data_stats = (abri_kraken2_arranged %>% 
                tbl_summary(
                  by = type,
                  include = c("GENE"),
                  statistic = list(
                    all_continuous() ~ "{mean} ({sd})",
                    all_categorical() ~ "{n} / {N} ({p}%)"
                  ),
                  digits = all_continuous() ~ 2,
                  #label = c(GENE ~ "MGEs"),
                  label = c(GENE ~ "Virulence Factors (VFs)"),
                  sort = list(everything() ~ "frequency"),
                  missing_text = "Missing",
                  missing = "no"
                )) %>%
  #add_p() %>%
  add_overall() %>%
  # --- This is the key step to group rows visually ---
  # Use the 'Taxa' variable to create header bands/groups in the table output
  modify_spanning_header(c("stat_0", "stat_1", "stat_2") ~ "**Location**") %>% # Optional: adjust main header
  # This function turns a variable into row group headers
  # CONVERT TO A {gt} TABLE! VERY IMPORTANT STEP!
  as_gt() %>%
  tab_header(md(glue("**Table 1. Virulence Factors (VFs) statistics results by location (_ANOVA p = {anova_p}_)**")))


# Print the summary table
print(data_stats)

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
