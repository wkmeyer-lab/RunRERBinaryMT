####################################################################################################################
####    Initial script written by Alexander Seaver                                                              ####
####    Select sections have been revised utilizing Google Gemini 3.6, manually reviewed by Alexander Seaver    ####
####################################################################################################################



# -------------------------------------------------------------------------
# 1. LOAD REQUIRED PACKAGES
# -------------------------------------------------------------------------
library(tidyverse)
library(igraph)
library(ggraph)
library(msigdbr)

# -------------------------------------------------------------------------
# 2. CONFIGURATION & LOADING
# -------------------------------------------------------------------------
raw_data <- read_csv("Simplify2GeneSets.csv")

# >>> SET YOUR DYNAMIC COLUMNS HERE <<<
target_column <- "Fig_2"
ivstat_column <- "IV_stat"

# -------------------------------------------------------------------------
# 3. FILTER AND PROCESS DATA DYNAMICALLY
# -------------------------------------------------------------------------
node_data <- raw_data %>%
  rename(pathway = 1) %>% 
  
  # 1. Dynamically keep rows active in your specified target column
  filter(.data[[target_column]] == 1) %>%
  
  # 2. Convert your specified iv_stat column to a strict binary category
  mutate(
    # Force to numeric to prevent text-based evaluation bugs
    iv_stat_num = as.numeric(.data[[ivstat_column]]),
    
    diet_direction = case_when(
      iv_stat_num < 0  ~ "Invertivore", 
      iv_stat_num > 0  ~ "Vertivore",   
      iv_stat_num == 0 ~ "Neutral",     # Catch exact zeros
      TRUE             ~ "Unknown"      # Catch NAs or unparseable text
    ),
    
    # Keep it as a plain character string so igraph doesn't corrupt it
    diet_direction = as.character(diet_direction),
    
    # Extract prefix database category (KEGG, REACTOME, PID)
    category = str_extract(pathway, "^[^_]+")
  ) %>%
  select(pathway, diet_direction, category, everything())

# -------------------------------------------------------------------------
# 4. FETCH CANONICAL GENES FROM MSigDB
# -------------------------------------------------------------------------
human_msigdb <- msigdbr(species = "Homo sapiens")

filtered_ref <- human_msigdb %>%
  filter(gs_name %in% node_data$pathway) %>%
  select(pathway = gs_name, gene = gene_symbol)

node_sizes <- filtered_ref %>%
  group_by(pathway) %>%
  summarise(total_reference_genes = n(), .groups = "drop")

node_data <- node_data %>%
  left_join(node_sizes, by = "pathway") %>%
  drop_na(total_reference_genes)

# -------------------------------------------------------------------------
# 5. CALCULATE OVERLAP (EDGES) BETWEEN PATHWAYS
# -------------------------------------------------------------------------
pathway_gene_lists <- filtered_ref %>%
  group_by(pathway) %>%
  summarise(genes = list(gene), .groups = "drop")

num_pathways <- nrow(pathway_gene_lists)
pairs <- combn(num_pathways, 2)

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

edge_data <- tibble(
  from = pathway_gene_lists$pathway[pairs[1, ]],
  to = pathway_gene_lists$pathway[pairs[2, ]],
  similarity = map2_dbl(pathway_gene_lists$genes[pairs[1, ]], 
                        pathway_gene_lists$genes[pairs[2, ]], 
                        jaccard)
)

edge_data_filtered <- edge_data %>% 
  filter(similarity > 0.15)

# -------------------------------------------------------------------------
# 6. BUILD AND PLOT THE NETWORK GRAPH
# -------------------------------------------------------------------------
# 1. Build the initial network object
initial_net <- graph_from_data_frame(d = edge_data_filtered, 
                                     vertices = node_data, 
                                     directed = FALSE)

# 2. Identify distinct clusters (islands) in the graph
cluster_assignments <- components(initial_net)$membership

# 3. Figure out which clusters contain at least one node
valid_clusters <- node_data %>%
  mutate(cluster_id = cluster_assignments[pathway]) %>%
  group_by(cluster_id) %>%
  summarise(has_large_node = any(total_reference_genes > 0)) %>%
  filter(has_large_node == TRUE) %>%
  pull(cluster_id)

# 4. Filter the node dataset to keep only nodes belonging to valid clusters
filtered_nodes <- node_data %>%
  mutate(cluster_id = cluster_assignments[pathway]) %>%
  filter(cluster_id %in% valid_clusters) %>%
  # Dynamically strip out the prefix (e.g. "KEGG_", "REACTOME_", "PID_")
  mutate(clean_name = str_remove(pathway, "^[^_]+_")) %>%
  # Group by the cluster ID to evaluate names locally within each island
  group_by(cluster_id) %>%
  mutate(
    # row_number() == which.min(...) flags only the single shortest string in the group
    is_shortest = row_number() == which.min(nchar(clean_name)),
    cluster_label = if_else(is_shortest, clean_name, "")
  ) %>%
  ungroup()

# 5. Filter the edge dataset to match only remaining nodes
filtered_edges <- edge_data_filtered %>%
  filter(from %in% filtered_nodes$pathway & to %in% filtered_nodes$pathway)

# 6. Rebuild the final filtered network object
net <- graph_from_data_frame(d = filtered_edges,
                             vertices = filtered_nodes,
                             directed = FALSE)

# >>> CRITICAL FIX: PRE-CALCULATE THE LAYOUT <<<
# Set a seed so you get the exact same shape if you run the script again tomorrow
set.seed(42) 

# Calculate the 'fr' layout once and save the coordinates
fixed_layout <- create_layout(net, layout = "fr")

# --- FIRST PLOT (With labels) ---
# Notice we pass `fixed_layout` instead of `net`
p <- ggraph(fixed_layout) + 
  geom_edge_link(aes(alpha = 1, width = similarity), 
                 color = "grey40", show.legend = TRUE) + 
  scale_edge_width(range = c(0.5, 2.5)) +
  
  geom_node_point(aes(size = total_reference_genes, 
                      color = diet_direction)) + 
  
  geom_node_text(aes(label = clean_name), 
                 repel = TRUE, 
                 size = 2.5, 
                 fontface = "bold",
                 max.overlaps = 100) +
  
  scale_color_manual(
    values = c(
      "Invertivore" = "#7570B3", 
      "Vertivore"   = "#E72A8B",
      "Neutral"     = "#E6AB02",
      "Unknown"     = "#999999" 
    ),
    name = "Diet"
  ) +
  
  scale_size_continuous(range = c(5, 16)) + 
  guides(edge_alpha = "none") +
  theme_void() +
  labs(
    size = "Total Genes in Reference",
    edge_width = "Gene Overlap"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 10)
  ) +
  guides(
    edge_alpha = "none",                  
    color = guide_legend(order = 1),      
    size = guide_legend(order = 2),       
    edge_width = guide_legend(order = 3)  
  )

ggsave("Figure2AllLegends.jpg", plot = p, width = 20, height = 20, dpi = 600)

# --- SECOND PLOT (No labels) ---
# Pass the exact same `fixed_layout` here
p2 <- ggraph(fixed_layout) + 
  geom_edge_link(aes(alpha = 1, width = similarity), 
                 color = "grey40", show.legend = TRUE) + 
  scale_edge_width(range = c(0.5, 2.5)) +
  
  geom_node_point(aes(size = total_reference_genes, 
                      color = diet_direction)) + 
  
  # label="" keeps the text layer empty but preserves any spatial mapping if needed
  geom_node_text(aes(label = ""), 
                 repel = TRUE, 
                 size = 2.5, 
                 fontface = "bold",
                 max.overlaps = 100) +
  
  scale_color_manual(
    values = c(
      "Invertivore" = "#7570B3", 
      "Vertivore"   = "#E72A8B",
      "Neutral"     = "#E6AB02",
      "Unknown"     = "#999999"
    ),
    name = "Diet"
  ) +
  
  scale_size_continuous(range = c(5, 16)) + 
  guides(edge_alpha = "none") +
  theme_void() +
  labs(
    size = "Total Genes in Reference",
    edge_width = "Gene Overlap"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 10)
  ) +
  guides(
    edge_alpha = "none",                  
    color = guide_legend(order = 1),      
    size = guide_legend(order = 2),       
    edge_width = guide_legend(order = 3)  
  )

ggsave("Figure2NoLabs.jpg", plot = p2, width = 20, height = 20, dpi = 600)