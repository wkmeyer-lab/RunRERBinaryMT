#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required to process phylogenetic trees. Install it with install.packages('ape').")
  }
  library(ape)
})

# ==============================================================================
# Helper Function: Clean Diet Strings (Locale-Safe)
# ==============================================================================
clean_diet_string <- function(x) {
  if (is.null(x)) return(x)
  x <- enc2utf8(as.character(x))
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  x <- gsub("^[_[:space:]]+", "", x)
  x <- gsub("[_[:space:]]+$", "", x)
  return(x)
}

# ==============================================================================
# Helper Function: Robust Tree Loader (.rds, .nwk, .phy)
# ==============================================================================
load_tree_file <- function(file_path) {
  ext <- tolower(tools::file_ext(file_path))
  
  if (ext == "rds") {
    obj <- readRDS(file_path)
    # Handle direct phylo object, or list containing a tree element
    if (inherits(obj, "phylo")) {
      return(obj)
    } else if (is.list(obj) && "tree" %in% names(obj) && inherits(obj$tree, "phylo")) {
      return(obj$tree)
    } else {
      stop("RDS object is not a valid 'phylo' tree.")
    }
  } else {
    return(ape::read.tree(file_path))
  }
}

# ==============================================================================
# 1. Parse CLI Arguments
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)

if ("--help" %in% args || "-h" %in% args) {
  cat("\nUsage: Rscript categoryCompare.R [path/to/tree1.rds path/to/tree2.nwk ...]\n")
  cat("Example: Rscript categoryCompare.R Trees/phenotypeTree.rds\n\n")
  quit(status = 0)
}

# ==============================================================================
# 2. Load Data and Clean Columns
# ==============================================================================
csv_path <- "../../../Data/mergedData.csv"

# Fallback to local Data directory if run from repository root
if (!file.exists(csv_path) && file.exists("Data/mergedData.csv")) {
  csv_path <- "Data/mergedData.csv"
}

if (!file.exists(csv_path)) {
  stop(sprintf("Error: Cannot find mergedData.csv at '%s'. Check your working directory.", csv_path))
}

mergedData <- read.csv(csv_path, stringsAsFactors = FALSE)

for (col in colnames(mergedData)) {
  if (is.character(mergedData[[col]]) || is.factor(mergedData[[col]])) {
    mergedData[[col]] <- clean_diet_string(mergedData[[col]])
  }
}

# Determine diet column for classification breakdown
diet_col <- if ("insVertivoreDiet" %in% colnames(mergedData)) {
  "insVertivoreDiet"
} else if ("SimplifiedDietConvertEqualDrop" %in% colnames(mergedData)) {
  "SimplifiedDietConvertEqualDrop"
} else {
  "Meyer.Lab.Classification"
}

# Determine tip identifier column
tip_id_col <- if ("ZoonomiaTip" %in% colnames(mergedData)) "ZoonomiaTip" else "Scientific_Binomial"

# Valid dataset species tips
data_species <- mergedData[[tip_id_col]]
data_species <- data_species[!is.na(data_species) & data_species != "" & data_species != "vs_NA"]

cat("==============================================================================\n")
cat(sprintf(" Loaded: %s\n", csv_path))
cat(sprintf(" Data Summary: %d valid species tips found using column '%s'\n", length(data_species), tip_id_col))
cat("==============================================================================\n")

# ==============================================================================
# 3. Compare with Phenotype Trees Passed via CLI
# ==============================================================================
if (length(args) == 0) {
  cat("\n[Note] No tree files provided via command line arguments.")
  cat("\nTo compare against tree files, pass file paths: Rscript categoryCompare.R <tree1.rds> <tree2.nwk>\n\n")
} else {
  cat(sprintf("\nComparing dataset against %d tree file(s)...\n\n", length(args)))
  
  for (tree_path in args) {
    cat("------------------------------------------------------------------------------\n")
    cat(sprintf("Tree File: %s\n", tree_path))
    cat("------------------------------------------------------------------------------\n")
    
    if (!file.exists(tree_path)) {
      cat(sprintf("  [ERROR] File not found: %s\n\n", tree_path))
      next
    }
    
    tree <- tryCatch(
      load_tree_file(tree_path),
      error = function(e) {
        cat(sprintf("  [ERROR] Failed to load tree: %s\n", e$message))
        return(NULL)
      }
    )
    
    if (is.null(tree)) next
    
    tree_tips <- tree$tip.label
    
    matching_tips <- intersect(data_species, tree_tips)
    missing_in_data <- setdiff(tree_tips, data_species)
    missing_in_tree <- setdiff(data_species, tree_tips)
    
    coverage_pct <- (length(matching_tips) / length(tree_tips)) * 100
    
    cat(sprintf("  Total Tree Tips:            %d\n", length(tree_tips)))
    cat(sprintf("  Matching Species in Data:   %d (%.1f%% coverage)\n", length(matching_tips), coverage_pct))
    cat(sprintf("  Tree Tips Missing in Data:  %d\n", length(missing_in_data)))
    cat(sprintf("  Data Tips Missing in Tree:  %d\n", length(missing_in_tree)))
    
    matched_rows <- mergedData[[tip_id_col]] %in% matching_tips
    tree_diet_counts <- table(mergedData[[diet_col]][matched_rows])
    
    cat("\n  Diet Phenotype Distribution on Matched Tree Tips:\n")
    print(tree_diet_counts)
    
    if (length(missing_in_data) > 0) {
      cat("\n  First 5 Tree Tips Missing from Dataset:\n   ")
      cat(head(missing_in_data, 5), sep = ", ")
      cat("\n")
    }
    cat("\n")
  }
}
