#Set the directory where your .Rds files are located
target_dir <- "/share/ceph/wym219group/shared/projects/seaverData/RunRERBinaryMT/Output/CategoricalInsVertivoreTree/Hyphy/"
1
#Set the pattern to match your files
file_pattern <- "-relax-.*\\.Rds$"

#Set a base name for your final output CSV file
output_file_base <- "/share/ceph/wym219group/shared/projects/seaverData/RunRERBinaryMT/Output/CategoricalInsVertivoreTree/Hyphy/hyphy_relax_p"


#Get a list of all files in the directory that match the pattern
all_files <- list.files(path = target_dir,
                        pattern = file_pattern,
                        full.names = TRUE)

#Check if any files are found
if (length(all_files) == 0) {
  stop(paste("No files found in", target_dir, "matching the pattern", file_pattern))
}

message(paste("Found", length(all_files), "files to process..."))

#Determine the output filename
first_filename <- basename(all_files[1])
fg_value <- gsub(".*-Foreground_(\\d+)\\.Rds$", "\\1", first_filename)

output_file <- ""
if (fg_value != first_filename) {
  output_file <- paste0(output_file_base, "_FG_", fg_value, ".csv")
  message(paste("Detected Foreground value", fg_value, "from first file."))
} else {
  output_file <- paste0(output_file_base, ".csv")
  message("Could not detect foreground value, using default output name.")
}
message(paste("Output will be saved to:", output_file))


#Create an empty list to store the results
results_list <- list()

#Loop through each file
for (file_path in all_files) {
  
  tryCatch({
    
    #Get just the filename (e.g., "Categorical...-GENENAME-Foreground_1.Rds")
    filename <- basename(file_path)
    
    #Regex to extract the gene name
    gene_name <- gsub(".*-relax-(.*?)-Foreground_.*", "\\1", filename)
    
    #Read Rds file and Extract p-value
    data <- readRDS(file_path)
    p_val <- data$p_value
    
    #Handle cases where the p-value wasn't found (is NULL)
    if (is.null(p_val)) {
      p_val <- NA # Set to NA (Not Available)
      message(paste("Warning: p.value not found in", filename, "(gene:", gene_name, ")"))
    }

    #Store as a data frame
    results_list[[file_path]] <- data.frame(
      Gene = gene_name,
      PValue = p_val
    )
    
  }, error = function(e) {
    #If anything fails (file read error, regex error, etc.),
    #print a warning and skip the file.
    message(paste("Error processing file:", file_path))
    message(paste("  Error was:", e$message))
  })
}

#Combine all the small data frames into one big data frame
final_results <- do.call(rbind, results_list)

# Clean up row names
rownames(final_results) <- NULL

# 6. Write the final data frame to a CSV file
write.csv(final_results, file = output_file, row.names = FALSE)

message("---")
message(paste("Successfully processed", nrow(final_results), "files."))
message(paste("Results saved to:", output_file))
