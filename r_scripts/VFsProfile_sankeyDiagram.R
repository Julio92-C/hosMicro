
# Load required libraries
library(plotly)
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(networkD3)
library(htmlwidgets)



# Load clean dataset
abri_kraken2_merged <- read_csv("hospital_Microbiome/Datasets/abri_kraken2_merged1.csv")


# Filter dataset by virulence factors
refactored_df <- abri_kraken2_merged %>%
  filter(DATABASE == "vfdb") %>%
  filter(type == "Air") %>%
  #filter(type == "Surface") %>%
  filter(`%IDENTITY` > 85) %>%
  filter(!grepl("Enterobacteriaceae|Rhizobium/Agrobacterium group|Bacillota|Bacteria|Bacteroidales|Bacteroidota", name)) %>%
  filter(str_count(name, "\\S+") > 1) %>%
  mutate(
    name = str_replace(name, "^(\\w)\\w+\\s(\\w+)", "\\1. \\2"),
    name = if_else(
      str_length(name) > 30,
      str_replace(name, "^((?:\\S+\\s+){2}).*", "\\1"),
      name
    )
  ) %>%
  filter(count > 0) %>%
  distinct() 

# Unique taxa carrying VFs  
uniTaxa <- length(unique(refactored_df$name))

# Define the function to extract the words
extract_productFunction <- function(string) {
  match <- regmatches(string, regexpr("\\- ([^\\(]+) \\(", string))
  result <- gsub("\\- | \\(", "", match)
  return(result)
}

# Apply the function to the dataframe
refactored_df <- refactored_df %>%
  mutate(Functions = sapply(PRODUCT, extract_productFunction))

# Function to create nodes dataframe
add_nodes <- function(data, nodes, col_name, color, border) {
  new_nodes <- data %>%
    select(name = {{col_name}}) %>%
    distinct() %>%
    arrange(name) %>%
    mutate(
      id = seq(max(nodes$id) + 1, max(nodes$id) + n()),
      color = color,
      border = border
    )
  
  nodes <- bind_rows(nodes, new_nodes)
  return(nodes)
}


# Create the nodes dataframe
nodes <- data.frame(
  name = sort(unique(refactored_df$sample)),
  id = seq(0, length(unique(refactored_df$sample)) - 1)
) %>%
  mutate(
    color = "#175709", 
    border = "black"
  )

# Add genus nodes to the nodes dataframe
nodes <- add_nodes(refactored_df, nodes, name, "#825cdb", "black")

# Add gene nodes to the nodes dataframe
nodes <- add_nodes(refactored_df, nodes, GENE, "#fc3503", "white")

# Add class nodes to the nodes dataframe
nodes <- add_nodes(refactored_df, nodes, Functions, "#b5b5b5", "black")


# Function to create links dataframe
create_links <- function(data, nodes, source_col, target_col, color) {
  links <- data %>%
    arrange(sample) %>%
    left_join(nodes %>% select(!!source_col := name, source = id), by = c(source_col)) %>%
    left_join(nodes %>% select(!!target_col := name, target = id), by = c(target_col)) %>%
    mutate(color = color) %>%
    select(source, target, color) %>%
    distinct(source, target, .keep_all = TRUE)
  
  return(links)
}

# Create the links dataframes
links0 <- create_links(refactored_df, nodes, "sample", "name", "#98ed85")
links1 <- create_links(refactored_df, nodes, "name", "GENE", "#4fc1e3")
links2 <- create_links(refactored_df, nodes, "GENE", "Functions", "#d1d0c24D")

# Merge the edges and remove duplicates
links <- bind_rows(links0, links1, links2) %>%
  distinct(source, target, .keep_all = TRUE)

# Add a new column 'value' with a value of 1
links <- links %>%
  arrange(source) %>%
  mutate(value = 1) 
# slice (2:n())


#Create the Sankey diagram
sankey_plot <- sankeyNetwork(Links = links, Nodes = nodes, Source = "source",
                             Target = "target", Value = "value", NodeID = "name",
                             fontSize = 14, fontFamily = "arial", width = 1500, height= 750, sinksRight=TRUE,
                             margin = list(top = 10, right = 10, bottom = 10, left = 10),
                             # iterations = 0,
                             #nodePadding = 30
)

# Plot
sankey_plot

# Save html plot
#saveWidget(sankey_plot, "hospital_Microbiome/Figures/Surface_VFsProfile_sankeyNetwork.html")
saveWidget(sankey_plot, "hospital_Microbiome/Figures/Air_VFsProfile_sankeyNetwork.html")
?sankeyNetwork

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


