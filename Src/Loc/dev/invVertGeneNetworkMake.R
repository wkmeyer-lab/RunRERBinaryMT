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
library(ggforce)

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
# 4. GMT FILE & REFERENCE GENE IMPORT
# -------------------------------------------------------------------------
gmt_file <- "KeggReactome.gmt"

# Helper function to read GMT files into a tidy dataframe
read_gmt_tidy <- function(file_path) {
  lines <- readLines(file_path)
  map_dfr(lines, function(line) {
    tokens <- unlist(strsplit(line, "\t"))
    pathway_name <- tokens[1]
    genes <- tokens[3:length(tokens)]
    genes <- genes[genes != "" & !is.na(genes)]
    tibble(pathway = pathway_name, gene = genes)
  })
}

# 1. Parse GMT file
ref_gmt <- read_gmt_tidy(gmt_file)

filtered_ref <- ref_gmt %>%
  filter(pathway %in% node_data$pathway) %>%
  select(pathway, gene)

# 2. Calculate gene set sizes
node_sizes <- filtered_ref %>%
  group_by(pathway) %>%
  summarise(total_reference_genes = n(), .groups = "drop")

# 3. Store target pathways before merging
expected_pathways <- node_data$pathway

node_data <- node_data %>%
  left_join(node_sizes, by = "pathway")

# 4. TryCatch validation block: Exit immediately if any target pathways are missing
tryCatch({
  missing_pathways <- setdiff(expected_pathways, node_data$pathway[!is.na(node_data$total_reference_genes)])
  
  if (length(missing_pathways) > 0) {
    stop(
      paste0(
        "Mismatch detected! The following ", length(missing_pathways),
        " pathway(s) marked for ", target_column, " were not found in the .gmt file or yielded 0 genes:\n",
        paste(" -", missing_pathways, collapse = "\n")
      )
    )
  }
}, error = function(e) {
  cat("\n=========================================================\n")
  cat("ERROR IN SECTION 4: PATHWAY VALIDATION FAILED\n")
  cat("=========================================================\n")
  cat(e$message, "\n\n")
  stop("Execution halted due to missing pathways in the GMT file.", call. = FALSE)
})

# Proceed with removing NAs after successful validation
node_data <- node_data %>% drop_na(total_reference_genes)

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
# 6. BUILD GRAPH & LAYOUT
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

# Calculate the 'fr' layout once and save the coordinates
set.seed(42) 
fixed_layout <- create_layout(net, layout = "fr")

# -------------------------------------------------------------------------
# 7. SIZING & OVERLAP RESOLVER
# -------------------------------------------------------------------------
min_point_size <- 20
max_point_size <- 40

x_span <- max(fixed_layout$x) - min(fixed_layout$x)
y_span <- max(fixed_layout$y) - min(fixed_layout$y)
coord_span <- max(x_span, y_span)

min_gene <- min(fixed_layout$total_reference_genes, na.rm = TRUE)
max_gene <- max(fixed_layout$total_reference_genes, na.rm = TRUE)

fixed_layout <- fixed_layout %>%
  mutate(
    gene_norm = if (max_gene == min_gene) 1 else (sqrt(total_reference_genes) - sqrt(min_gene)) / (sqrt(max_gene) - sqrt(min_gene)),
    size_mm = min_point_size + gene_norm * (max_point_size - min_point_size),
    radius = (size_mm / 2) * (coord_span / 600)
  )

# Nudge overlapping nodes apart while keeping layout topology intact
min_gap <- coord_span * 0.01  # Spacing buffer between node borders

for (iter in 1:100) {
  moved <- FALSE
  n_nodes <- nrow(fixed_layout)
  for (i in 1:(n_nodes - 1)) {
    for (j in (i + 1):n_nodes) {
      dx <- fixed_layout$x[j] - fixed_layout$x[i]
      dy <- fixed_layout$y[j] - fixed_layout$y[i]
      dist <- sqrt(dx^2 + dy^2)
      min_dist <- fixed_layout$radius[i] + fixed_layout$radius[j] + min_gap
      
      if (dist < min_dist) {
        if (dist == 0) {
          dx <- runif(1, -0.01, 0.01)
          dy <- runif(1, -0.01, 0.01)
          dist <- sqrt(dx^2 + dy^2)
        }
        overlap <- min_dist - dist
        shift_x <- (dx / dist) * (overlap / 2)
        shift_y <- (dy / dist) * (overlap / 2)
        
        fixed_layout$x[i] <- fixed_layout$x[i] - shift_x
        fixed_layout$y[i] <- fixed_layout$y[i] - shift_y
        fixed_layout$x[j] <- fixed_layout$x[j] + shift_x
        fixed_layout$y[j] <- fixed_layout$y[j] + shift_y
        moved <- TRUE
      }
    }
  }
  if (!moved) break
}

# Define palette & shared theme settings
diet_colors <- c(
  "Invertivore" = "#7570B3", 
  "Vertivore"   = "#E72A8B",
  "Neutral"     = "#E6AB02",
  "Unknown"     = "#999999"
)

base_theme <- theme_void() + 
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 10)
  )

# -------------------------------------------------------------------------
# 8. PLOTS
# -------------------------------------------------------------------------

# --- FIRST PLOT (With labels) ---
p <- ggraph(fixed_layout) + 
  coord_fixed() + 
  
  geom_edge_link(aes(width = similarity), color = "grey20", show.legend = TRUE) + 
  scale_edge_width(range = c(2.0, 4.0)) +
  
  # Invisible points layer to maintain size legend scale
  geom_node_point(aes(size = total_reference_genes), alpha = 0) +
  
  # 1. Single-segment solid circle filled by diet direction
  geom_circle(data = fixed_layout, 
              aes(x0 = x, y0 = y, r = radius, fill = diet_direction), 
              color = NA) +
  
  # 2. Outer single black stroke overlay (no fill)
  geom_circle(data = fixed_layout,
              aes(x0 = x, y0 = y, r = radius),
              color = "black", fill = NA, linewidth = 0.8) +
  
  geom_node_text(aes(label = clean_name), 
                 repel = TRUE, 
                 size = 2.5, 
                 fontface = "bold",
                 max.overlaps = 100) +
  
  scale_fill_manual(values = diet_colors, name = "Diet") +
  scale_size_continuous(range = c(5, 16), name = "Total Genes in Reference") +
  
  base_theme +
  labs(
    fill = "Diet",
    size = "Total Genes in Reference",
    edge_width = "Gene Overlap"
  ) +
  guides(
    fill = guide_legend(order = 1),      
    size = guide_legend(order = 2),       
    edge_width = guide_legend(order = 3)  
  )

ggsave("Figure2AllLegends.jpg", plot = p, width = 20, height = 20, dpi = 300)

# --- SECOND PLOT (No labels) ---
p2 <- ggraph(fixed_layout) + 
  coord_fixed() + 
  
  geom_edge_link(aes(width = similarity), color = "grey20", show.legend = TRUE) + 
  scale_edge_width(range = c(2.0, 4.0)) +
  
  geom_node_point(aes(size = total_reference_genes), alpha = 0) +
  
  # 1. Single-segment solid circle filled by diet direction
  geom_circle(data = fixed_layout, 
              aes(x0 = x, y0 = y, r = radius, fill = diet_direction), 
              color = NA) +
  
  # 2. Outer single black stroke overlay (no fill)
  geom_circle(data = fixed_layout,
              aes(x0 = x, y0 = y, r = radius),
              color = "black", fill = NA, linewidth = 0.8) +
  
  geom_node_text(aes(label = ""), 
                 repel = TRUE, 
                 size = 2.5, 
                 fontface = "bold",
                 max.overlaps = 100) +
  
  scale_fill_manual(values = diet_colors, name = "Diet") +
  scale_size_continuous(range = c(5, 16), name = "Total Genes in Reference") +
  
  base_theme +
  labs(
    fill = "Diet",
    size = "Total Genes in Reference",
    edge_width = "Gene Overlap"
  ) +
  guides(
    fill = guide_legend(order = 1),      
    size = guide_legend(order = 2),       
    edge_width = guide_legend(order = 3)  
  )

ggsave("Figure2NoLabs.jpg", plot = p2, width = 20, height = 20, dpi = 300)
