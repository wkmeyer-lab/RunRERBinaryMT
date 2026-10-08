library(jsonlite)

# 1. Setup paths
input_dir <- "/share/ceph/wym219group/shared/projects/seaverProjects/RunRERBinaryMT/Output/CategoricalInsVertivoreTree/Hyphy/"
file_list <- list.files(path = input_dir, 
                        pattern = "CategoricalInsVertivoreTree-Hyphy-relax-.*\\.json$", 
                        full.names = TRUE)

if (length(file_list) == 0) stop("No JSON files found in the directory.")

# 2. Extraction Function
process_relax_json <- function(file_path) {
  # Load JSON safely
  data <- tryCatch(fromJSON(file_path), error = function(e) return(NULL))
  if (is.null(data)) return(NULL)
"/share/ceph/wym219group/shared/projects/seaverProjects/RunRERBinaryMT/Output/CategoricalInsVertivoreTree/Hyphy/"  
  # Navigate to "test results"
  stats <- data[["test results"]]
  
  # Check if the specific "relaxation..." key exists
  k_val <- stats[["relaxation or intensification parameter"]]
  if (is.null(k_val)) return(NULL)
  
  # Extract metadata from filename
  fname <- basename(file_path)
  clean_name <- gsub("CategoricalInsVertivoreTree-Hyphy-relax-", "", fname)
  clean_name <- gsub("\\.json$", "", clean_name)
  
  # Split Gene and Foreground
  parts <- strsplit(clean_name, "-Foreground_")[[1]]
  gene <- parts[1]
  fg_val <- if(length(parts) > 1) parts[2] else "Unknown"
  
  return(data.frame(
    Gene = gene,
    Foreground = fg_val,
    K_parameter = as.numeric(k_val),
    LRT = as.numeric(stats[["LRT"]]),
    p_value = as.numeric(stats[["p-value"]]),
    stringsAsFactors = FALSE
  ))
}

# 3. Execution
cat("Processing", length(file_list), "files...\n")
results_list <- lapply(file_list, process_relax_json)

# Filter NULLs and combine
results_list <- results_list[!vapply(results_list, is.null, logical(1))]

if (length(results_list) > 0) {
  final_table <- do.call(rbind, results_list)
  
  # Multiple testing correction
  final_table$p_adj <- p.adjust(final_table$p_value, method = "fdr")
  
  # Interpretation
  final_table$Result <- "NS"
  final_table$Result[final_table$p_adj < 0.05 & final_table$K_parameter < 1] <- "Relaxed"
  final_table$Result[final_table$p_adj < 0.05 & final_table$K_parameter > 1] <- "Intensified"
  
  # Save
  write.csv(final_table, "/share/ceph/wym219group/shared/projects/seaverProjects/RunRERBinaryMT/Output/CategoricalInsVertivoreTree/Hyphy/RELAX_Summary_Final.csv", row.names = FALSE)
  cat("Success! Generated table with", nrow(final_table), "rows.\n")
} else {
  cat("Error: Still no data extracted. Double-check the key names in one JSON file.\n")
}
