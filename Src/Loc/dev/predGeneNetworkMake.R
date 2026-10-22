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
library(ggforce)

# -------------------------------------------------------------------------
# 2. CONFIGURATION & LOADING
# -------------------------------------------------------------------------
raw_data <- read_csv("Simplify2GeneSets.csv")
target_column <- "Fig_4"

p_thresh <- 0.05         # Adjusted p-value threshold
num_sig_thresh <- 50     # Minimum number of significant alternates

# -------------------------------------------------------------------------
# 3. FILTER AND PROCESS DATA
# -------------------------------------------------------------------------
node_data <- raw_data %>%
  rename(pathway = 1) %>%
  filter(.data[[target_column]] == 1) %>%
  
  # Ensure all numerical statistical columns and count columns are parsed correctly
  mutate(across(
    any_of(c("CH_stat", "HI_stat", "HV_stat",
             "CH_padj", "HI_padj", "HV_padj",
             "CH_PadjNumSignificant", "HI_PadjNumSignificant", "HV_PadjNumSignificant")),
    as.numeric
  )) %>%
  mutate(
    # Dual significance criteria check per comparison (padj < 0.05 AND PadjNumSignificant >= 50)
    is_sig_CH = !is.na(CH_padj) & CH_padj < p_thresh &
      !is.na(CH_PadjNumSignificant) & CH_PadjNumSignificant >= num_sig_thresh &
      (is.na(CH_stat) | CH_stat < 0),

    is_sig_HI = !is.na(HI_padj) & HI_padj < p_thresh &
      !is.na(HI_PadjNumSignificant) & HI_PadjNumSignificant >= num_sig_thresh &
      (is.na(HI_stat) | HI_stat > 0),

    is_sig_HV = !is.na(HV_padj) & HV_padj < p_thresh &
      !is.na(HV_PadjNumSignificant) & HV_PadjNumSignificant >= num_sig_thresh &
      (is.na(HV_stat) | HV_stat > 0),

    # Calculate log significance for node alpha
    CH_sig = if_else(is_sig_CH, -log10(CH_padj + 1e-10), 0),
    HI_sig = if_else(is_sig_HI, -log10(HI_padj + 1e-10), 0),
    HV_sig = if_else(is_sig_HV, -log10(HV_padj + 1e-10), 0),
    dom_sig = pmax(-log10(CH_padj + 1e-10), -log10(HI_padj + 1e-10), -log10(HV_padj + 1e-10), na.rm = TRUE),

    # Combined overlap category for multi-group nodes
    overlap_category = case_when(
      is_sig_HI & is_sig_HV & is_sig_CH ~ "HI + HV + CH",
      is_sig_HI & is_sig_HV             ~ "HI + HV",
      is_sig_HI & is_sig_CH             ~ "HI + CH",
      is_sig_HV & is_sig_CH             ~ "HV + CH",
      is_sig_HI                         ~ "HI",
      is_sig_HV                         ~ "HV",
      is_sig_CH                         ~ "CH",
      TRUE                              ~ "None"
    ),

    # Single category shading based on the LOWEST p-value comparison among passing thresholds
    dominant_comparison = case_when(
      !is_sig_HI & !is_sig_HV & !is_sig_CH ~ "None",
      HI_padj == pmin(if_else(is_sig_HI, HI_padj, Inf),
                      if_else(is_sig_HV, HV_padj, Inf),
                      if_else(is_sig_CH, CH_padj, Inf), na.rm = TRUE) ~ "HI",
      HV_padj == pmin(if_else(is_sig_HI, HI_padj, Inf),
                      if_else(is_sig_HV, HV_padj, Inf),
                      if_else(is_sig_CH, CH_padj, Inf), na.rm = TRUE) ~ "HV",
      CH_padj == pmin(if_else(is_sig_HI, HI_padj, Inf),
                      if_else(is_sig_HV, HV_padj, Inf),
                      if_else(is_sig_CH, CH_padj, Inf), na.rm = TRUE) ~ "CH",
      TRUE ~ "None"
    )
  )



# -------------------------------------------------------------------------
# 4. GMT FILE & OVERLAP (EDGES)
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
        " pathway(s) marked for Fig_4 were not found in the .gmt file or yielded 0 genes:\n",
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

# 5. Extract gene lists and calculate pairwise Jaccard similarity (edges)
pathway_gene_lists <- filtered_ref %>%
  group_by(pathway) %>%
  summarise(genes = list(gene), .groups = "drop")

pairs <- combn(nrow(pathway_gene_lists), 2)
jaccard <- function(a, b) { length(intersect(a, b)) / length(union(a, b)) }

edge_data <- tibble(
  from = pathway_gene_lists$pathway[pairs[1, ]],
  to = pathway_gene_lists$pathway[pairs[2, ]],
  similarity = map2_dbl(pathway_gene_lists$genes[pairs[1, ]],
                        pathway_gene_lists$genes[pairs[2, ]], jaccard)
) %>% filter(similarity > 0.15)

# -------------------------------------------------------------------------
# 5. BUILD GRAPH & LAYOUT
# -------------------------------------------------------------------------
initial_net <- graph_from_data_frame(d = edge_data, vertices = node_data, directed = FALSE)
cluster_assignments <- components(initial_net)$membership

filtered_nodes <- node_data %>%
  mutate(cluster_id = cluster_assignments[pathway]) %>%
  group_by(cluster_id) %>%
  mutate(cluster_size = n()) %>%
  filter(cluster_size >= 1, any(total_reference_genes > 0)) %>% 
  mutate(clean_name = str_remove(pathway, "^[^_]+_")) %>%
  ungroup()

filtered_edges <- edge_data %>%
  filter(from %in% filtered_nodes$pathway & to %in% filtered_nodes$pathway)

net <- graph_from_data_frame(d = filtered_edges, vertices = filtered_nodes, directed = FALSE)

# Use 'stress' layout for improved, uniform node placement
# Revert to 'fr' layout (classic network shape)
fixed_layout <- create_layout(net, layout = "fr")

# -------------------------------------------------------------------------
# LOCK FIXED VIEWPORT LIMITS
# -------------------------------------------------------------------------
# Add 20% padding around the node coordinates so labels fit inside locked bounds
x_pad <- (max(fixed_layout$x) - min(fixed_layout$x)) * 0.20
y_pad <- (max(fixed_layout$y) - min(fixed_layout$y)) * 0.20

fixed_xlim <- c(min(fixed_layout$x) - x_pad, max(fixed_layout$x) + x_pad)
fixed_ylim <- c(min(fixed_layout$y) - y_pad, max(fixed_layout$y) + y_pad)

# -------------------------------------------------------------------------
# SIZING & OVERLAP RESOLVER
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

# Nudge overlapping nodes apart while keeping 'fr' topology intact
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

# -------------------------------------------------------------------------
# PIE CHART & SINGLE NODE DATA PREPARATION
# -------------------------------------------------------------------------
node_categories <- fixed_layout %>%
  select(name, x, y, radius, total_reference_genes, dom_sig, clean_name, overlap_category, dominant_comparison,
         is_sig_CH, is_sig_HI, is_sig_HV) %>%
  pivot_longer(
    cols = c(is_sig_CH, is_sig_HI, is_sig_HV),
    names_to = "comparison",
    values_to = "is_active"
  ) %>%
  mutate(
    group = case_when(
      comparison == "is_sig_CH" ~ "CH",
      comparison == "is_sig_HI" ~ "HI",
      comparison == "is_sig_HV" ~ "HV"
    )
  ) %>%
  filter(is_active == TRUE)

split_pie_data <- node_categories %>%
  group_by(name) %>%
  mutate(num_slices = n()) %>%
  filter(num_slices >= 2) %>%
  mutate(
    slice_idx = row_number(),
    start_angle = (slice_idx - 1) * (2 * pi / num_slices),
    end_angle = slice_idx * (2 * pi / num_slices)
  ) %>%
  ungroup()

multi_node_names <- unique(split_pie_data$name)

single_nodes <- fixed_layout %>%
  filter(!name %in% multi_node_names)

base_theme <- theme_void() + 
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16),
    legend.title = element_text(face = "bold", size = 10)
  )

# -------------------------------------------------------------------------
# 6. PLOTS
# -------------------------------------------------------------------------
plot_master_pie <- ggraph(fixed_layout) + 
  coord_fixed() + 
  geom_edge_link(aes(width = similarity), color = "grey20") + 
  scale_edge_width(range = c(2.0, 4.0)) +
  
  # 1. Single-category solid nodes (using coordinate-space geom_circle)
  geom_circle(data = single_nodes, 
              aes(x0 = x, y0 = y, r = radius, fill = overlap_category), 
              color = NA) +
  
  # 2. Multi-category pie slices
  geom_arc_bar(data = split_pie_data,
               aes(x0 = x, y0 = y, r0 = 0, r = radius, 
                   start = start_angle, end = end_angle, 
                   fill = group, alpha = 1),
               color = "white", linewidth = 0.4) +
  
  # 3. Outer dark halo underlay (ALL nodes)
  geom_circle(data = fixed_layout,
              aes(x0 = x, y0 = y, r = radius * 1.03),
              color = "black", fill = NA, linewidth = 4.5) +
  
  # 4. Outer dominant comparison stroke (ALL nodes)
  geom_circle(data = fixed_layout,
              aes(x0 = x, y0 = y, r = radius * 1.03, color = dominant_comparison),
              fill = NA, linewidth = 2.5) +
  
  geom_node_text(aes(label = if_else(total_reference_genes >= 1, clean_name, "")), 
                 repel = TRUE, size = 3, fontface = "bold") +
  
  scale_color_manual(values = category_colors, guide = "none") +
  scale_fill_manual(values = category_colors, name = "Comparison") +
  scale_alpha_continuous(range = 1) +
  scale_size_continuous(range = c(20, 40), name = "Category Size") +
  base_theme

# Colorblind-Friendly Palette
category_colors <- c(
  # Individual comparisons (specified)
  "HI"           = "#7570B3",
  "HV"           = "#E72A8B",
  "CH"           = "#E6AB02",
  
  # Overlaps
  "HI + HV"      = "#F4A582",
  "HI + CH"      = "#00A1F2",
  "HV + CH"      = "#D95F02",
  "HI + HV + CH" = "#E63946",
  
  # Non-significant
  "None"         = "grey80"
)

# -------------------------------------------------------------------------
# VERSION 1A: PIE CHART WITH LABELS
# -------------------------------------------------------------------------
plot_master_pie <- ggraph(fixed_layout) + 
  coord_fixed() + 
  geom_edge_link(aes(width = similarity), color = "grey20") + 
  scale_edge_width(range = c(2.0, 4.0)) +
  
  # 1. Single-category solid nodes (using coordinate-space geom_circle)
  geom_circle(data = single_nodes, 
              aes(x0 = x, y0 = y, r = radius, fill = overlap_category), 
              color = NA) +
  
  # 2. Multi-category pie slices
  geom_arc_bar(data = split_pie_data,
               aes(x0 = x, y0 = y, r0 = 0, r = radius, 
                   start = start_angle, end = end_angle, 
                   fill = group, alpha = 1),
               color = "white", linewidth = 0.4) +
  
  # 3. Outer dark halo underlay (ALL nodes)
  geom_circle(data = fixed_layout,
              aes(x0 = x, y0 = y, r = radius * 1.03),
              color = "black", fill = NA, linewidth = 4.5) +
  
  # 4. Outer dominant comparison stroke (ALL nodes)
  geom_circle(data = fixed_layout,
              aes(x0 = x, y0 = y, r = radius * 1.03, color = dominant_comparison),
              fill = NA, linewidth = 2.5) +
  
  geom_node_text(aes(label = if_else(total_reference_genes >= 1, clean_name, "")), 
                 repel = TRUE, size = 3, fontface = "bold") +
  
  scale_color_manual(values = category_colors, guide = "none") +
  scale_fill_manual(values = category_colors, name = "Comparison") +
  scale_alpha_continuous(range = 1) +
  scale_size_continuous(range = c(20, 40), name = "Category Size") +
  base_theme

# -------------------------------------------------------------------------
# VERSION 1B: PIE CHART WITHOUT LABELS
# -------------------------------------------------------------------------
plot_master_pie_nolab <- ggraph(fixed_layout) + 
  coord_fixed() + 
  geom_edge_link(aes(width = similarity), color = "grey20") + 
  scale_edge_width(range = c(2.0, 4.0)) +
  
  # 1. Single-category solid nodes (using coordinate-space geom_circle)
  geom_circle(data = single_nodes, 
              aes(x0 = x, y0 = y, r = radius, fill = overlap_category), 
              color = NA) +
  
  # 2. Multi-category pie slices
  geom_arc_bar(data = split_pie_data,
               aes(x0 = x, y0 = y, r0 = 0, r = radius, 
                   start = start_angle, end = end_angle, 
                   fill = group, alpha = 1),
               color = "white", linewidth = 0.4) +
  
  # 3. Outer dark halo underlay (ALL nodes)
  geom_circle(data = fixed_layout,
              aes(x0 = x, y0 = y, r = radius * 1.03),
              color = "black", fill = NA, linewidth = 4.5) +
  
  # 4. Outer dominant comparison stroke (ALL nodes)
  geom_circle(data = fixed_layout,
              aes(x0 = x, y0 = y, r = radius * 1.03, color = dominant_comparison),
              fill = NA, linewidth = 2.5) +
  
  
  scale_color_manual(values = category_colors, guide = "none") +
  scale_fill_manual(values = category_colors, name = "Comparison") +
  scale_alpha_continuous(range = 1) +
  scale_size_continuous(range = c(20, 40), name = "Category Size") +
  base_theme

# -------------------------------------------------------------------------
# VERSION 2A: OVERLAP CATEGORY COMBINATIONS WITH LABELS
# -------------------------------------------------------------------------
plot_master_overlap <- ggraph(fixed_layout) +
  coord_fixed() +
  geom_edge_link(aes(width = similarity), color = "grey20") +
  scale_edge_width(range = c(2.0, 4.0)) +

  geom_node_point(aes(size = total_reference_genes, color = overlap_category, alpha = 1)) +
  geom_node_text(aes(label = if_else(total_reference_genes >= 1, clean_name, "")),
                 repel = TRUE, size = 3, fontface = "bold") +

  scale_color_manual(values = category_colors, name = "Overlap Category") +
  scale_alpha_continuous(range = 1) +
  scale_size_continuous(range = c(20, 40), name = "Category Size") +
  base_theme

# -------------------------------------------------------------------------
# VERSION 2B: OVERLAP CATEGORY COMBINATIONS WITHOUT LABELS
# -------------------------------------------------------------------------
plot_master_overlap_nolab <- ggraph(fixed_layout) +
  coord_fixed() +
  geom_edge_link(aes(width = similarity), color = "grey20") +
  scale_edge_width(range = c(2.0, 4.0)) +

  geom_node_point(aes(size = total_reference_genes, color = overlap_category, alpha = 1)) +
  geom_node_text(aes(label = "")) +

  scale_color_manual(values = category_colors, name = "Overlap Category") +
  scale_alpha_continuous(range = 1) +
  scale_size_continuous(range = c(20, 40), name = "Category Size") +
  base_theme

# -------------------------------------------------------------------------
# VERSION 3A: DOMINANT COMPARISON (SINGLE SHADING BY LOWEST P-VALUE) WITH LABELS
# -------------------------------------------------------------------------
plot_master_dominant <- ggraph(fixed_layout) +
  coord_fixed() +
  geom_edge_link(aes(width = similarity), color = "grey20") +
  scale_edge_width(range = c(2.0, 4.0)) +

  geom_node_point(aes(size = total_reference_genes, color = dominant_comparison, alpha = 1)) +
  geom_node_text(aes(label = if_else(total_reference_genes >= 1, clean_name, "")),
                 repel = TRUE, size = 3, fontface = "bold") +

  scale_color_manual(values = category_colors, name = "Dominant Comparison\n(Lowest p-value)") +
  scale_alpha_continuous(range = 1) +
  scale_size_continuous(range = c(20, 40), name = "Category Size") +
  base_theme

# -------------------------------------------------------------------------
# VERSION 3B: DOMINANT COMPARISON (SINGLE SHADING BY LOWEST P-VALUE) WITHOUT LABELS
# -------------------------------------------------------------------------
plot_master_dominant_nolab <- ggraph(fixed_layout) +
  coord_fixed() +
  geom_edge_link(aes(width = similarity), color = "grey20") +
  scale_edge_width(range = c(2.0, 4.0)) +

  geom_node_point(aes(size = total_reference_genes, color = dominant_comparison, alpha = 1)) +
  geom_node_text(aes(label = "")) +

  scale_color_manual(values = category_colors, name = "Dominant Comparison\n(Lowest p-value)") +
  scale_alpha_continuous(range = 1) +
  scale_size_continuous(range = c(20, 40), name = "Category Size") +
  base_theme

# Save all 6 plot variations
ggsave("Fig4_PieChart.jpg", plot_master_pie, height = 20, width = 20, dpi = 300)
ggsave("Fig4_PieChart_nolab.jpg", plot_master_pie_nolab, height = 20, width = 20, dpi = 300)
# ggsave("Fig4_OverlapCategory.jpg", plot_master_overlap, height = 20, width = 20, dpi = 300)
# ggsave("Fig4_OverlapCategory_nolab.jpg", plot_master_overlap_nolab, height = 20, width = 20, dpi = 300)
# ggsave("Fig4_Dominant.jpg", plot_master_dominant, height = 20, width = 20, dpi = 300)
# ggsave("Fig4_Dominant_nolab.jpg", plot_master_dominant_nolab, height = 20, width = 20, dpi = 300)