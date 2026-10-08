#!/usr/bin/env Rscript

# Generated Using Gemini Pro 3.6

# Set graphics renderer for headless HPC Linux nodes (prevents Fontconfig errors)
options(bitmapType = "cairo")

suppressPackageStartupMessages({
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required. Install with install.packages('ape').")
  }
  library(ape)
})

# ==============================================================================
# Helper Functions
# ==============================================================================
clean_diet_string <- function(x) {
  if (is.null(x)) return(x)
  x <- enc2utf8(as.character(x))
  x <- gsub("\u00A0", " ", x, fixed = TRUE)
  x <- gsub("^[_[:space:]]+", "", x)
  x <- gsub("[_[:space:]]+$", "", x)
  return(x)
}

load_tree_obj <- function(path) {
  obj <- readRDS(path)
  if (inherits(obj, "phylo")) return(obj)
  if (is.list(obj) && "tree" %in% names(obj) && inherits(obj$tree, "phylo")) return(obj$tree)
  stop("Invalid tree object in file: ", path)
}

# Extract tip phenotype states from tree edges
get_tip_states <- function(tree) {
  tip_indices <- 1:length(tree$tip.label)
  edge_indices <- match(tip_indices, tree$edge[, 2])
  states <- tree$edge.length[edge_indices]
  names(states) <- tree$tip.label
  return(states)
}

# Decode numeric tip state integers (1..N) to trait strings using sorted dataset factor levels
decode_tree_states <- function(states, col_name, df) {
  if (is.numeric(states)) {
    ref_vals <- sort(unique(na.omit(df[[col_name]])))
    mapped_chars <- ref_vals[states]
  } else {
    mapped_chars <- as.character(states)
  }
  names(mapped_chars) <- names(states)
  return(mapped_chars)
}

# Helper functional mappings
map_to_pantheria <- function(v) {
  ifelse(v %in% c("Invertivore", "Insectivore", "Vertivore", "zMixedPredator", "MixedPredator"), "Carnivore", v)
}

map_to_walker <- function(v) {
  ifelse(v %in% c("Insectivore", "Planktivore"), "Invertivore",
  ifelse(v %in% c("Carnivore", "Piscivore", "Hematophagy"), "Vertivore", v))
}

# ==============================================================================
# 1. Parse CLI Arguments & Load Data
# ==============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("\nUsage: Rscript threeWayConcordance.R tree1.rds tree2.rds tree3.rds\n\n")
  quit(status = 1)
}

csv_path <- "../../../Data/mergedData.csv"
if (!file.exists(csv_path) && file.exists("Data/mergedData.csv")) {
  csv_path <- "Data/mergedData.csv"
}

if (!file.exists(csv_path)) {
  stop(sprintf("Cannot find mergedData.csv at '%s'. Check working directory.", csv_path))
}

mergedData <- read.csv(csv_path, stringsAsFactors = FALSE)
for (col in colnames(mergedData)) {
  if (is.character(mergedData[[col]]) || is.factor(mergedData[[col]])) {
    mergedData[[col]] <- clean_diet_string(mergedData[[col]])
  }
}

tip_col <- if ("ZoonomiaTip" %in% colnames(mergedData)) "ZoonomiaTip" else "Scientific_Binomial"
pantheria_col <- "panTheriaTrophicLevelCharacter"
walker_col <- if ("Meyer.Lab.Classification.Compressed" %in% colnames(mergedData)) {
  "Meyer.Lab.Classification.Compressed"
} else {
  "Meyer.Lab.Classification"
}

# Load input trees and simplify display names
tree_names <- gsub("CategoricalTree\\.rds$", "", basename(args))
tree_names <- gsub("^ComplexDietCentralAnalysis", "", tree_names)
trees <- lapply(args, load_tree_obj)
names(trees) <- tree_names

col_mapping <- list(
  SimplifyStrictPred  = "SimplifiedDietConvertStrictPred",
  SimplifyEqualDrop   = "SimplifiedDietConvertEqualDrop",
  Simplify2           = "insVertivoreDiet"
)

# ==============================================================================
# 2. Pairwise Inter-Tree Concordance
# ==============================================================================
cat("==============================================================================\n")
cat(" 1. THREE-WAY PAIRWISE INTER-TREE CONCORDANCE\n")
cat("==============================================================================\n\n")

pair_names <- c()
pair_matches <- c()
pair_mismatches <- c()
pair_pcts <- c()

for (i in 1:(length(trees) - 1)) {
  for (j in (i + 1):length(trees)) {
    t1 <- get_tip_states(trees[[i]])
    t2 <- get_tip_states(trees[[j]])
    
    valid <- !is.na(t1) & !is.na(t2) & t1 != "" & t2 != ""
    match_cnt <- sum(t1[valid] == t2[valid])
    total_cnt <- sum(valid)
    mismatch_cnt <- total_cnt - match_cnt
    pct <- if (total_cnt > 0) (match_cnt / total_cnt) * 100 else 0
    
    p_label <- sprintf("%s\nvs %s", names(trees)[i], names(trees)[j])
    pair_names <- c(pair_names, p_label)
    pair_matches <- c(pair_matches, match_cnt)
    pair_mismatches <- c(pair_mismatches, mismatch_cnt)
    pair_pcts <- c(pair_pcts, pct)
    
    cat(sprintf("  %s  vs  %s: %d / %d (%.2f%% agreement)\n", 
                names(trees)[i], names(trees)[j], match_cnt, total_cnt, pct))
  }
}

# ==============================================================================
# 3. External Database Concordance (PanTHERIA & Walker's Mammals)
# ==============================================================================
cat("\n==============================================================================\n")
cat(" 2. EXTERNAL DATABASE CONCORDANCE\n")
cat("==============================================================================\n\n")

first_tree_tips <- names(get_tip_states(trees[[1]]))
df_idx <- match(first_tree_tips, mergedData[[tip_col]])

pantheria_raw <- mergedData[[pantheria_col]][df_idx]
walker_raw    <- mergedData[[walker_col]][df_idx]
walker_mapped <- map_to_walker(walker_raw)

model_names <- names(trees)
pan_matches <- c(); pan_mismatches <- c(); pan_pcts <- c()
walk_matches <- c(); walk_mismatches <- c(); walk_pcts <- c()

for (name in model_names) {
  raw_states <- get_tip_states(trees[[name]])
  col_name <- if (name %in% names(col_mapping)) col_mapping[[name]] else "insVertivoreDiet"
  
  decoded <- decode_tree_states(raw_states, col_name, mergedData)
  
  # --- A. PanTHERIA Concordance ---
  decoded_pan <- map_to_pantheria(decoded)
  pan_valid <- !is.na(decoded_pan) & !is.na(pantheria_raw) & decoded_pan != "" & pantheria_raw != ""
  p_match <- sum(decoded_pan[pan_valid] == pantheria_raw[pan_valid])
  p_total <- sum(pan_valid)
  p_mismatch <- p_total - p_match
  p_pct <- if (p_total > 0) (p_match / p_total) * 100 else 0
  
  pan_matches <- c(pan_matches, p_match)
  pan_mismatches <- c(pan_mismatches, p_mismatch)
  pan_pcts <- c(pan_pcts, p_pct)
  
  # --- B. Walker's Mammals Concordance ---
  decoded_walk <- map_to_walker(decoded)
  w_valid <- !is.na(decoded_walk) & !is.na(walker_mapped) & decoded_walk != "" & walker_mapped != "" & decoded_walk != "zMixedPredator"
  w_match <- sum(decoded_walk[w_valid] == walker_mapped[w_valid])
  w_total <- sum(w_valid)
  w_mismatch <- w_total - w_match
  w_pct <- if (w_total > 0) (w_match / w_total) * 100 else 0
  
  walk_matches <- c(walk_matches, w_match)
  walk_mismatches <- c(walk_mismatches, w_mismatch)
  walk_pcts <- c(walk_pcts, w_pct)
  
  cat(sprintf("=== Tree Model: %s ===\n", name))
  cat(sprintf("  vs PanTHERIA (%s):         %d / %d (%.2f%% concordance)\n", pantheria_col, p_match, p_total, p_pct))
  cat(sprintf("  vs Walker's Mammals (%s): %d / %d (%.2f%% concordance)\n\n", walker_col, w_match, w_total, w_pct))
}

# ==============================================================================
# 4. Plot & Save 3-Panel Stacked Bar Chart
# ==============================================================================
out_dir <- "../../../Output"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}
png_path <- file.path(out_dir, "concordance.png")

png(png_path, width = 3200, height = 1100, res = 200, type = "cairo-png")
par(mfrow = c(1, 3), mar = c(6, 5, 4, 2))

# Panel A: Pairwise Inter-Tree Stacked Bar
pairwise_mat <- rbind(Match = pair_matches, Mismatch = pair_mismatches)
bp1 <- barplot(
  pairwise_mat,
  names.arg = pair_names,
  col = c("#2B5C8F", "#D95F02"),
  main = "A. Pairwise Inter-Tree Concordance",
  ylab = "Number of Species / Edges",
  ylim = c(0, max(pair_matches + pair_mismatches) * 1.15),
  las = 1,
  cex.names = 0.85
)
legend("bottomright", legend = c("Concordant (Match)", "Discordant (Mismatch)"),
       fill = c("#2B5C8F", "#D95F02"), bty = "n", cex = 0.9)

for (k in seq_along(pair_pcts)) {
  text(bp1[k], pair_matches[k] / 2, sprintf("%.1f%%", pair_pcts[k]), col = "white", font = 2, cex = 1.1)
}

# Panel B: PanTHERIA Concordance Stacked Bar
pan_mat <- rbind(Match = pan_matches, Mismatch = pan_mismatches)
bp2 <- barplot(
  pan_mat,
  names.arg = model_names,
  col = c("#2B5C8F", "#D95F02"),
  main = "B. PanTHERIA Benchmark",
  ylab = "Number of Evaluated Species",
  ylim = c(0, max(pan_matches + pan_mismatches) * 1.15),
  las = 1,
  cex.names = 0.95
)
legend("bottomright", legend = c("Concordant (Match)", "Discordant (Mismatch)"),
       fill = c("#2B5C8F", "#D95F02"), bty = "n", cex = 0.9)

for (k in seq_along(pan_pcts)) {
  text(bp2[k], pan_matches[k] / 2, sprintf("%.1f%%", pan_pcts[k]), col = "white", font = 2, cex = 1.1)
}

# Panel C: Walker's Mammals Concordance Stacked Bar
walk_mat <- rbind(Match = walk_matches, Mismatch = walk_mismatches)
bp3 <- barplot(
  walk_mat,
  names.arg = model_names,
  col = c("#2B5C8F", "#D95F02"),
  main = "C. Walker's Mammals Benchmark",
  ylab = "Number of Evaluated Species",
  ylim = c(0, max(walk_matches + walk_mismatches) * 1.15),
  las = 1,
  cex.names = 0.95
)
legend("bottomright", legend = c("Concordant (Match)", "Discordant (Mismatch)"),
       fill = c("#2B5C8F", "#D95F02"), bty = "n", cex = 0.9)

for (k in seq_along(walk_pcts)) {
  text(bp3[k], walk_matches[k] / 2, sprintf("%.1f%%", walk_pcts[k]), col = "white", font = 2, cex = 1.1)
}

dev.off()
cat(sprintf("[SUCCESS] 3-panel concordance figure saved to: %s\n\n", png_path))
