#!/usr/bin/env Rscript

# Usage: 
# Rscript violinPlotMake.R <rer_matrix.rds> <trees_object.rds> <phenotype_tree.rds_or_txt> <genes_list.txt>

# 1. Load Required Libraries
suppressPackageStartupMessages({
  library(phangorn)
  library(phytools)
  library(ggplot2)
  library(RERconverge)
})

# 2. Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Missing arguments.\nUsage: Rscript violinPlotMake.R <rer_matrix.rds> <trees_object.rds> <phenotype_tree.rds_or_txt> <genes_list.txt>")
}

rer_file      <- args[1]
trees_file    <- args[2]
phentree_file <- args[3]
genes_file    <- args[4]

# 3. Load Data Objects
message("Loading data objects...")
mamrer   <- readRDS(rer_file)
mamtrees <- readRDS(trees_file)

# Support loading phenotype tree from either RDS or raw text/Newick format
if (grepl("\\.rds$", phentree_file, ignore.case = TRUE)) {
  phentree <- readRDS(phentree_file)
} else {
  phentree <- read.tree(phentree_file)
}

# Read gene list and drop empty lines
genes_to_plot <- readLines(genes_file)
genes_to_plot <- genes_to_plot[trimws(genes_to_plot) != ""]

message(sprintf("Found %d genes to plot.", length(genes_to_plot)))

# 4. Calculate phenotype paths
message("Mapping phenotype paths...")
phenvdiet <- tree2Paths(phentree, mamtrees)

# Define explicit categorical order
diet_levels <- c("Herbivore", "Omnivore", "Vertivore", "Invertivore")

# Define custom color palette
diet_colors <- c(
  "Herbivore"   = "#1B9E77",
  "Omnivore"    = "black",
  "Vertivore"   = "#E72A8B",
  "Invertivore" = "#7570B3"
)

# 5. Loop over each gene to extract and plot
for (gene in genes_to_plot) {
  message(sprintf("Processing %s...", gene))
  
  # Validate gene exists in both objects
  if (!(gene %in% names(mamtrees$trees))) {
    warning(sprintf("Gene %s not found in trees object. Skipping.", gene))
    next
  }
  if (!(gene %in% rownames(mamrer))) {
    warning(sprintf("Gene %s not found in RER matrix. Skipping.", gene))
    next
  }
  
  gene_tree <- mamtrees$trees[[gene]]
  
  # Prune phenotype tree to match gene tree tips
  pphentree <- drop.tip(phentree, setdiff(phentree$tip.label, gene_tree$tip.label))
  
  # Get RERs mapped to tree branches
  rertree <- returnRersAsTree(mamtrees, mamrer, index = gene, phenv = phenvdiet)
  
  # Set up dummy edge lengths to find foreground/background indices
  gene_tree_temp <- gene_tree
  gene_tree_temp$edge.length <- rep(2, nrow(gene_tree_temp$edge))
  
  # Call internal RERconverge helper via ::: namespace operator
  ee <- RERconverge:::edgeIndexRelativeMaster(gene_tree_temp, mamtrees$masterTree)
  ii <- mamtrees$matIndex[ee[, c(2,1)]]
  
  RelativeRate <- rertree$edge.length
  Phenotype <- rep(NA, length(RelativeRate))
  
  # Map multi-class phenotype to labels
  Phenotype[phenvdiet[ii] == 1] <- "Herbivore"
  Phenotype[phenvdiet[ii] == 2] <- "Invertivore"
  Phenotype[phenvdiet[ii] == 3] <- "Omnivore"
  Phenotype[phenvdiet[ii] == 4] <- "Vertivore"
  
  rtoplot <- data.frame(RelativeRate, Phenotype)
  rtoplot <- rtoplot[!is.na(rtoplot$Phenotype), ]
  
  # Set explicit factor levels to enforce plotting order: Herbivore -> Omnivore -> Vertivore -> Invertivore
  rtoplot$Phenotype <- factor(rtoplot$Phenotype, levels = diet_levels)
  
  # Plotting
  p <- ggplot(rtoplot, aes(x = Phenotype, y = RelativeRate, col = Phenotype)) +
    geom_violin(adjust = 1/3) +
    geom_jitter(position = position_jitter(0.2)) +
    theme_classic() +
    scale_color_manual(values = diet_colors, limits = diet_levels) +
    theme(text = element_text(size = 20)) +
    labs(title = paste("RER Distribution:", gene),
         y = "Relative Rate",
         x = "Diet Phenotype")
  
  output_pdf <- paste0(gene, "_RelativeRatesByPhenotype.pdf")
  ggsave(output_pdf, plot = p, width = 7, height = 6)
  message(sprintf(" -> Saved %s", output_pdf))
}

message("All plots generated successfully.")
