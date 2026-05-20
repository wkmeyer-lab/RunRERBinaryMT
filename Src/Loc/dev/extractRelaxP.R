#load jsonlite
if (!require("jsonlite", quietly = TRUE)) {
  local_lib <- Sys.getenv("R_LIBS_USER")
  if (!dir.exists(local_lib)) dir.create(local_lib, recursive = TRUE)
  install.packages("jsonlite", lib = local_lib, repos = "http://cran.us.r-project.org")
  library("jsonlite", lib.loc = local_lib)
}

#set target directory
target_dir <- "/share/ceph/wym219group/shared/projects/seaverProjects/RunRERBinaryMT/Output/CategoricalInsVertivoreTree/Hyphy/"

#set pattern to match JSON files
file_pattern <- "-relax-.*\\.json$"

#set output base path
output_file_base <- "/share/ceph/wym219group/shared/projects/seaverProjects/RunRERBinaryMT/Output/CategoricalInsVertivoreTree/Hyphy/hyphy_relax_p"

#get list of files
all_files <- list.files(path = target_dir, pattern = file_pattern, full.names = TRUE)

if (length(all_files) == 0) {
  stop(paste("No JSON files found in", target_dir))
}

message(paste("Found", length(all_files), "JSON files. Starting extraction..."))

results_list <- list()
skipped_count <- 0

#loop through files
for (file_path in all_files) {
  
  filename <- basename(file_path)

  #check file size
  info <- file.info(file_path)
  if (is.na(info$size) || info$size < 10) {
      skipped_count <- skipped_count + 1
      next
  }

  tryCatch({
    
    #regex genename
    gene_name <- gsub(".*-relax-(.*?)-Foreground_.*", "\\1", filename)
    
    #regex foreground value
    fg_value <- gsub(".*-Foreground_(\\d+)\\.json$", "\\1", filename)
    if (fg_value == filename) fg_value <- "Unknown"

    #read json
    json_data <- fromJSON(file_path)
    
    #extract values from "test results"
    if ("test results" %in% names(json_data)) {
        tr <- json_data[["test results"]]
        
        lrt_val <- tr[["LRT"]]
        raw_p   <- tr[["p-value"]]
        raw_k   <- tr[["relaxation or intensification parameter"]]
        
        #store data with consistent naming
        results_list[[file_path]] <- data.frame(
          Gene = gene_name,
          LRT = lrt_val,
          p_val = raw_p,
          k_parameter = raw_k,
          Foreground = fg_value,
          stringsAsFactors = FALSE
        )
        
    } else {
        skipped_count <- skipped_count + 1
    }

  }, error = function(e) {
    skipped_count <<- skipped_count + 1
  })
}

#combine Results
if (length(results_list) == 0) {
  stop("All files were empty or invalid. No data extracted.")
}

final_results <- do.call(rbind, results_list)
rownames(final_results) <- NULL

#apply holm-bonferroni - creates p_val_adj
final_results$p_val_adj <- p.adjust(final_results$p_val, method = "holm")

#--- Extraction Summary ---#
message("--- Extraction Complete. Processing per Foreground... ---")
message(paste("Total Files Scanned:", length(all_files)))
message(paste("Successfully Extracted:", nrow(final_results)))
message(paste("Skipped/Corrupt Files:", skipped_count))

#write CSVs
unique_fgs <- unique(final_results$Foreground)

for (fg in unique_fgs) {
  #subset the data for the specific foreground
  subset_data <- final_results[final_results$Foreground == fg, ]
  
  #run Holm-Bonferroni correction ONLY on this subset
  subset_data$p_val_adj <- p.adjust(subset_data$p_val, method = "holm")
  
  #select and clean columns
  subset_data_clean <- subset_data[, c("Gene", "LRT", "p_val", "p_val_adj", "k_parameter")]
  
  #sort by the adjusted p-value
  subset_data_clean <- subset_data_clean[order(subset_data_clean$p_val_adj), ]
  
  #save the file
  this_output_file <- paste0(output_file_base, "_FG_", fg, ".csv")
  write.csv(subset_data_clean, file = this_output_file, row.names = FALSE)
  
  message(paste("Saved:", this_output_file, "(", nrow(subset_data_clean), "genes corrected independently)"))
}

message("--- Done ---")
