#Read your .rds file
combinedGenes <- readRDS("/share/ceph/wym219group/shared/projects/seaverData/RunRERBinaryMT/Output/CategoricalInsVertivoreTree/CategoricalInsvertivoreTreecombinedGeneResults.rds")
#Define all the columns you want to check
significant_cols <- c(
  "HI-significant", "CH-significant", "HO-significant",
  "HV-significant", "IO-significant", "IV-significant", "OV-significant"
)

#Find all row indices where *any* of those columns are TRUE
list_of_rows <- lapply(significant_cols, function(col_name) {
  
  #Make sure the column actually exists in the data
  if (col_name %in% names(combinedGenes)) {
    
    # Get row numbers where the column is TRUE
    which(combinedGenes[[col_name]] == TRUE)
    
  } else {
    # If the column doesn't exist, return an empty vector
    integer(0)
  }
})

#Combine all unique indices into a single vector
all_matching_rows <- Reduce(union, list_of_rows)

#Get the row names corresponding to those unique rows
unique_significant_rownames <- rownames(combinedGenes)[all_matching_rows]

#Write just this list of row names to a text file
writeLines(
  unique_significant_rownames,
  con = "/share/ceph/wym219group/shared/projects/seaverData/RunRERBinaryMT/Output/CategoricalInsVertivoreTree/combinedSignificantGenes.txt"
)

print(paste("Found", length(unique_significant_rownames), "unique row names."))
print("Written to /share/ceph/wym219group/shared/projects/seaverData/RunRERBinaryMT/Output/CategoricalInsVertivoreTree/combinedSignificantGenes.txt")
