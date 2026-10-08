# 1. SETUP PATHS
input_dir <- "/share/ceph/wym219group/shared/projects/seaverProjects/RunRERBinaryMT/Output/CategoricalInsVertivoreTree/Hyphy/"
file_list <- list.files(path = input_dir, pattern = "\\.json$", full.names = TRUE)

if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite", repos="https://cloud.r-project.org/")
if (!requireNamespace("pheatmap", quietly = TRUE)) install.packages("pheatmap", repos="https://cloud.r-project.org/")
library(jsonlite)
library(pheatmap)

# 2. EXTRACTION FUNCTION
process_relax_json <- function(file_path) {
  data <- tryCatch(fromJSON(file_path), error = function(e) return(NULL))
  if (is.null(data)) return(NULL)
  
  stats <- data[["test results"]]
  k_val <- stats[["relaxation or intensification parameter"]]
  if (is.null(k_val)) return(NULL)
  
  fname <- basename(file_path)
  clean_name <- gsub("CategoricalInsVertivoreTree-Hyphy-relax-", "", fname)
  clean_name <- gsub("\\.json$", "", clean_name)
  parts <- strsplit(clean_name, "-Foreground_")[[1]]
  
  fg_id <- if(length(parts) > 1) parts[2] else "Unknown"
  
  # Mapping with explicit numeric ordering
  fg_name <- switch(fg_id,
                    "1" = "Herbivory",
                    "2" = "Invertivory",
                    "3" = "Omnivory",
                    "4" = "Vertivory",
                    paste0("Group_", fg_id))
  
  return(data.frame(
    Gene = parts[1],
    Foreground = fg_name,
    K = as.numeric(k_val),
    p = as.numeric(stats[["p-value"]]),
    stringsAsFactors = FALSE
  ))
}

# 3. PROCESS DATA
cat("Processing", length(file_list), "files...\n")
results_list <- lapply(file_list, process_relax_json)
df <- do.call(rbind, results_list)

# 4. CALCULATE SCORES & ORDERING
df$p_adj <- p.adjust(df$p, method = "fdr")
df$score <- sign(log10(df$K)) * -log10(df$p_adj + 1e-10)

# DEFINE EXPLICIT ORDER (Updated: Herb -> Omni -> Invert -> Vert)
diet_order <- c("Herbivory", "Omnivory", "Invertivory", "Vertivory")

# 5. PIVOT TO MATRIX
genes <- unique(df$Gene)
mat <- matrix(0, nrow = length(genes), ncol = length(diet_order), 
              dimnames = list(genes, diet_order))

for(i in 1:nrow(df)) {
  # Only populate if the category is in our ordered list
  if (df$Foreground[i] %in% diet_order) {
    mat[df$Gene[i], df$Foreground[i]] <- df$score[i]
  }
}

# 6. GENERATE PLOTS DIRECTLY TO INPUT DIRECTORY
output_all <- paste0(input_dir, "RELAX_Heatmap_All_Genes.pdf")
output_top <- paste0(input_dir, "RELAX_Heatmap_Top_Hits.pdf")
output_csv <- paste0(input_dir, "RELAX_Final_Summary_Table.csv")

my_colors <- colorRampPalette(c("dodgerblue4", "white", "firebrick3"))(100)
# Create a descriptive title to explain the unit values on the key
legend_desc <- "Score = sign(log10 K) * -log10(p_adj) | Red > 0: Intensified, Blue < 0: Relaxed"

# FIGURE A: All Genes
pdf(output_all, width = 8, height = 11)
pheatmap(mat, 
         cluster_cols = FALSE, # Keep Herbivory -> Vertivory order
         cluster_rows = TRUE,  # Cluster genes by similar patterns
         show_rownames = FALSE, 
         color = my_colors, 
         main = paste("Global Evolutionary Landscape\n", legend_desc), 
         angle_col = 45)
dev.off()

# FIGURE B: Top Hits
sig_genes <- unique(df$Gene[df$p_adj < 0.05])
if(length(sig_genes) > 1) {
  mat_top <- mat[sig_genes, , drop = FALSE]
  pdf(output_top, width = 10, height = max(6, length(sig_genes) * 0.18))
  pheatmap(mat_top, 
           cluster_cols = FALSE, # Keep Herbivory -> Vertivory order
           cluster_rows = TRUE, 
           show_rownames = TRUE, 
           fontsize_row = 7,
           color = my_colors, 
           main = paste("Significant Evolutionary Shifts\n", legend_desc), 
           angle_col = 45)
  dev.off()
}

write.csv(df, output_csv, row.names = FALSE)
cat("Success! Files saved to:", input_dir, "\n")
