# Load library
library(tibble)
library(dplyr)
library(igraph)
library(ggraph)
library(readr)
library(stringr)
library(vegan)
library(tidygraph) # Assuming your graph 'g' is a tbl_graph object or similar
library(ggrepel)   # For geom_node_text(repel = TRUE)
library(ggforce)   # For geom_mark_hull (if you want to mark clusters)
library(funrar) # transform an abundance matrix into a relative abundance matrix


# Load dataset
microbiome_df <- read_csv("Hospital_microbiome/Datasets/abri_kraken2_merged1.csv")
# colnames(abri_kraken2_merged)

# # Filter by Species and counts greater than 30
df_species <- microbiome_df %>%
  group_by(name) %>%
  filter(!grepl("Terrabacteria group|Bacteroidota/Chlorobiota group|FCB group", name)) %>%
  filter(str_count(name, "\\S+") > 1) %>%
  mutate(
    name = str_replace(name, "^(\\w)\\w+\\s(\\w+)", "\\1. \\2"),
    name = if_else(
      str_length(name) > 30,
      str_replace(name, "^((?:\\S+\\s+){2}).*", "\\1"),
      name
    )
  ) %>%
  #summarise(count = sum(count), .groups = "drop") %>%
  filter(count > 100) %>%
  distinct() %>%
  rename(taxa = name)

# Unique taxa
uniqueTaxa <- length(unique(df_species$taxa))

# 1. Prepare Edges (including new links between Samples and Treatments)

# Weighted Sample → Taxa edges
sample_taxa_edges <- df_species %>%
  select(from = sample, to = taxa, weight = count)

# Unweighted Taxa → Gene edges (can be binary or frequency-based later)
taxa_gene_edges <- df_species %>%
  select(from = taxa, to = GENE) %>%
  mutate(weight = 1)

# Combine edges
edges <- bind_rows(sample_taxa_edges, taxa_gene_edges)

# 2. Prepare Node Metadata (including Treatment Groups as their own nodes)

# Sample metadata
sample_meta <- df_species %>%
  distinct(sample, type, location, ward)

# Gene metadata
gene_meta <- df_species %>%
  group_by(GENE) %>%
  select(GENE, DATABASE) %>%
  rename(gene = GENE) %>%
  rename(geneCategory = DATABASE) %>%
  distinct()  %>% 
  mutate(geneCategory = case_when(
    grepl("card", geneCategory) ~ "AMR",
    grepl("vfdb", geneCategory) ~ "VFs",
    grepl("plasmidfinder", geneCategory) ~ "MGEs",
    TRUE ~ geneCategory
  ))
  

# Build node list
nodes <- unique(c(edges$from, edges$to)) %>%
  tibble(name = .) %>%
  mutate(category = case_when(
    name %in% df_species$sample ~ "sample",
    name %in% df_species$taxa   ~ "taxa",
    name %in% df_species$GENE   ~ "gene"
  )) %>%
  left_join(sample_meta, by = c("name" = "sample")) %>%
  left_join(gene_meta,   by = c("name" = "gene")) %>%
  distinct(name, .keep_all = TRUE) # remove duplicates

# 3. (Optional) Your Bray-Curtis Clustering analysis
# This analysis is valid for grouping SAMPLES, but is separate from the network itself
# and just adds metadata to the node table.

# Create a sample × taxa matrix from your microbiome_df
abundance_matrix <- df_species %>%
  select(sample, taxa, count) %>%
  # mutate(
  #   count = (count / sum(count)) * 100 # Convert the count to relative abundance towards to normal distributions
  # ) %>%
  tidyr::pivot_wider(
    names_from = taxa,
    values_from = count,
    values_fill = 0,
    values_fn = sum # Use `sum` to aggregate duplicate values
  ) %>%
  # distinct(sample, .keep_all = TRUE) %>%
  column_to_rownames("sample")

# Convert the data frame to a matrix
abundance_matrix <- as.matrix(abundance_matrix)

# transform an abundance matrix into a relative abundance matrix
relative_abundance_matrix <- abundance_matrix %>%
  make_relative()

# Compute Bray-Curtis distance
dist_matrix <- vegdist(relative_abundance_matrix, method = "bray")

# Hierarchical clustering
hc <- hclust(dist_matrix, method = "average")

# Cut tree into clusters (e.g. 2 groups)
sample_clusters <- cutree(hc, k = 4)

# Add cluster info to node table
nodes$distance_cluster <- sample_clusters[match(nodes$name, names(sample_clusters))]


## 4. Prepare tables for Gephi Export

# Clean edge table for Gephi
edges_gephi <- edges %>%
  rename(Source = from, Target = to, Weight = weight)

summary(edges_gephi)

# Clean node table for Gephi
nodes_gephi <- nodes %>%
  rename(Id = name) %>%
  select(Id, category, type, location, ward, geneCategory, distance_cluster)


# Replace NA values with an empty string ""
nodes_gephi$type[is.na(nodes_gephi$type)] <- ""
nodes_gephi$location[is.na(nodes_gephi$location)] <- ""
nodes_gephi$ward[is.na(nodes_gephi$ward)] <- ""
nodes_gephi$geneCategory[is.na(nodes_gephi$geneCategory)] <- ""
nodes_gephi$distance_cluster[is.na(nodes_gephi$distance_cluster)] <- ""

# Class 
summary(nodes_gephi)

# Export node and edges to CSV files
write.csv(edges_gephi, "Hospital_microbiome/Datasets/gephi_edges2.csv", row.names = FALSE)
write.csv(nodes_gephi, "Hospital_microbiome/Datasets/gephi_nodes2.csv", row.names = FALSE)


# Create graph object
g <- graph_from_data_frame(edges, vertices = nodes, directed = FALSE)

# --- Start of your modified plotting code ---

# Plot with gene nodes colored by GeneCategory
ggraph(g, layout = "fr",
       niter = 1000, # Increase iterations for more stable layout, potentially better spread
       area = vcount(g)^2 * 1.5) + # Increase 'area' parameter for Fruchterman-Reingold
  # to potentially spread nodes further apart
  
  geom_edge_link(aes(width = weight), 
                 alpha = 0.2) + # Reduce edge opacity to thin out dense areas
  # (original was 0.4, try 0.2 or even lower like 0.1)
  
  geom_node_point(aes(
    shape = type,
    color = case_when(
      type == "sample" ~ as.factor(distance_cluster),
      type == "gene"   ~ geneCategory,
      TRUE             ~ type
    )
  ), 
  size = 5,
  alpha = 0.7) + # Add alpha to node points for transparency, especially for overlap
  # Adjust 'alpha' to your preference (e.g., 0.7-0.9)
  
  geom_node_text(aes(label = name), 
                 repel = TRUE, 
                 size = 3, 
                 max.overlaps = Inf,
                 segment.color = 'grey50', # Add a segment line for clarity when labels are repelled
                 min.segment.length = 0.5) + # Only draw segments if they are long enough
  
  scale_edge_width(range = c(0.5, 3)) +
  scale_color_manual(values = c(
    "1" = "#1f77b4", "2" = "#ff7f0e", "3" = "yellow", "4"="green", # Sample clusters
    AMR         = "#d62728",
    VFs         = "#9467bd",
    MGEs        = "#2ca02c",
    taxa        = "#8c564b"
  )) +
  scale_shape_manual(values = c(sample = 16, taxa = 17, gene = 15)) +
  
  # --- Optional: Add geom_mark_hull to highlight clusters if needed ---
  # Assuming 'distance_cluster' or 'geneCategory' defines your central clusters
  # You might want to filter this to specific clusters if the entire center is a mix
  # For example, to highlight 'gene' nodes by their 'geneCategory':
  # geom_mark_hull(aes(x = x, y = y, fill = geneCategory, filter = type == "gene"),
  #                alpha = 0.1, show.legend = FALSE) +
  
  theme_void() +
  theme(legend.position = "none")

# --- End of your modified plotting code ---


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
