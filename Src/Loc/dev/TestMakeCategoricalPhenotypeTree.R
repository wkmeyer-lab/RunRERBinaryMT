# ================================================================================
# INITIALIZATION & LIBRARIES
# ================================================================================
clusterRun = T
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} 
library(RERconverge)
library(tools)

setwd("/share/ceph/wym219group/shared/projects/seaverProjects/RunRERBinaryMT/")
source("Src/Reu/ZoonomTreeNameToCommon.R")

# ================================================================================
# 1. HARDCODED ARGUMENTS
# ================================================================================
filePrefix <- "CategoricalInsVertivoreTreeLiamInference"
forceUpdate <- TRUE
mainTreesLocation <- "Data/zoonomiaAllMammalsTrees.rds"
spreadSheetLocation <- "Data/mergedData.csv"
annotColumn <- "DerekDietClassification90InsVertivoreSorting"
nameColumn <- "ZoonomiaTip"
modelType <- "ER"
ancestralTrait <- NULL
useScreen <- FALSE
minimumBranchLength <- 0.01

categoryList <- c(
  "C-Invertebrate-eater", "C-Endotherm-Carnivore", "C-Herpetivore", "C-Piscivore", "C-Nonspecific-Vertebrate-eater", "C-Scavenger", 
  "O-For Examination", "O-Scavenger", 
  "H-Frugivore", "H-Nectarivore", "H-Granivore", "H-Nonspecific-Herbivore", 
  "C-Terrestrial-vertebrates-eater", "C-All-vertebrate-eater", "C-All-Animals-Eater", 
  "H-High-sugar-plants-Eater", "H-Low-sugar-plants-Eater", "H-All-plants-Eater", 
  "O-Generalist", 
  "C-InsVertivore-Mixed", "C-InsVertivore-Piscivore", "C-InsVertivore-Insectivore","C-InsVertivore-Carnivore",
  "Insectivore", "Herpetivore", "Piscivore", "Vertivore", "InsVertivore", "Omnivore", "Frugivore", "Nectarivore", "Glucivore", "Herbivore", "Generalist"
)

substitutions <- list(
  c("C-Invertebrate-eater", "Insectivore"), c("C-InsVertivore-Insectivore", "Insectivore"),
  c("C-Herpetivore", "Vertivore"),
  c("C-Piscivore", "Vertivore"), c("C-InsVertivore-Piscivore", "Vertivore"),
  c("C-Endotherm-Carnivore", "Vertivore"), c("C-Scavenger", "Vertivore"), c("C-Nonspecific-Vertebrate-eater", "Vertivore"),
  c("C-Terrestrial-vertebrates-eater", "Vertivore"), c("C-All-vertebrate-eater", "Vertivore"), c("C-InsVertivore-Carnivore", "Vertivore"),
  c("C-InsVertivore-Mixed", "Omnivore"), 
  c("O-For Examination", "Omnivore"), c("O-Scavenger", "Omnivore"),
  c("H-Frugivore", "Herbivore"), 
  c("H-Nectarivore", "Herbivore"), 
  c("H-High-sugar-plants-Eater", "Herbivore"),
  c("H-Granivore", "Herbivore"), c("H-Nonspecific-Herbivore", "Herbivore"), 
  c("H-Low-sugar-plants-Eater", "Herbivore"), c("H-All-plants-Eater", "Herbivore"),
  c("O-Generalist", "Omnivore")
)

manualUnprunedTips <- c(
  "vs_HLornAna3", "vs_HLtacAcu1",  "PreserveMonotremeBranches",
  "vs_HLdidVir1", "vs_HLgymLea1", "vs_HLpseCup1", "MarsupialTransitionPreservation",
  "vs_HLmyrTri1", "vs_HLchoDid1", "vs_HLchoHof3", "vs_HLproCap3", "AfrotheriaPreserveTransitions",
  "vs_HLpanLeo1", "LionClade", "vs_HLpanOnc1", "vs_HLaciJub2", "CheetahClade", "PreserveBigCats", 
  "vs_HLursThi1", "vs_ursMar1", "vs_HLursArc1", "vs_HLailMel2", "UrsaPreserveTransition",
  "vs_lepWed1", "vs_HLmirAng2", "vs_HLphoVit1", "vs_HLeriBar1", "SealPreserveTransitions", 
  "vs_HLodoRos1", "vs_HLcalUrs1", "vs_HLzalCal1", "SealSeaLionPreserveTransitions",
  "vs_HLmelCap1", "vs_HLgulGul1", "vs_HLneoVis1", "MustelidPreserveTransitions",
  "vs_HLpteBra1", "vs_HLlutLut1", "vs_enhLutKen1", "OtterPreserveTransitions",
  "vs_HLlycPic2", "CanidPreserveTransision",
  "vs_HLgloMel1", "vs_HLpepEle1", "vs_HLturAdu1", "DolphinClade", "vs_orcOrc1", "vs_HLescRob1", "vs_HLlniGeo1", "amazonRiverDolphinFromYeast","vs_HLbalEde1", "vs_HLmegNov1", "vs_HLcynGun1", "CetaceaPreserveTransitions",
  "vs_HLmerUng1", "BankVoleTransition",
  "vs_HLeulMon1", "vs_HLeulFul1", "LemurTransition",
  "vs_eulMac1", "vs_HLeulFla1", "vs_ponAbe3", "PrimateTransitions",
  "vs_panTro6", "vs_HLrhiRox2", "LangurClade", "vs_HLallNig1", "PrimateTransitionsContinued",
  "vs_HLeryPat1", "vs_chlSab2", "geunonClade", "vs_HLtheGel1", "PrimateTransitionsContinuedAgain",
  "vs_HLpapAnu5", "vs_HLmanSph1", "DrillMandrillClade", "vs_cerAty1", "Drilltransitions",
  "vs_HLtraJav1", "MouseDeerUsedInOtherAnalysis",
  "vs_mm10", "humans",
  "vs_HLmarFla1", "marmotClade", "DoormouseTransition"
)

manualPrunedTips <- c(
  "vs_HLellTal1", "vs_HLellLut1", "vs_HLarvAmp1", "voleClade",
  "vs_HLmusSpi1", "vs_HLmusCar1", "vs_HLmasCou1", "vs_HLmusPah1", "mouseClade",
  "vs_HLhysCri1", "vs_HLthrSwi1", "vs_HLpetTyp1", "vs_hetGla2", "vs_chiLan1", "vs_HLdinBra1", "vs_HLcteSoc1", "vs_octDeg1", "vs_HLcoePre1", "vs_HLdasPun1", "vs_HLdolPat1", "gundiGuineaPigClade",
  "vs_HLoryGaz1", "vs_HLbeaHun1", "vs_HLkobLecLec1", "vs_HLkobLecLec1", "vs_HLmadKir1", "vs_HLneoPyg1", "vs_HLphiMax1", "vs_HLoreOre1", "vs_HLneoMos1", "vs_HLaepMel1", "vs_HLtraImb1", "Bovidae",
  "vs_HLhydIne1", "vs_HLmunMun1", "Cervidae",
  "vs_HLtraKan1", "mouseDeerOtherIsKept",
  "vs_HLmurAurFea1", "outerVespert",
  "vs_HLmyoLuc1", "Nearctic",
  "vs_myoDav1", "Myotis",
  "vs_HLpipPip1", "vs_HLlasBor1", "vs_HLnycHum2", "Vespertilioninae",
  "vs_HLmacSob1", "FoxLongTounge",
  "vs_HLeidHel2", "outerPeropodidae",
  "vs_HLeonSpe1", "Roussetinae"
)

# ================================================================================
# 2. OUTPUT DIRECTORY SETUP
# ================================================================================
if(!dir.exists("Output")){ dir.create("Output", recursive = TRUE) }
outputFolderNameNoSlash = paste("Output/randomTrees/",filePrefix, sep = "")
if(!dir.exists(outputFolderNameNoSlash)){ dir.create(outputFolderNameNoSlash, recursive = TRUE) }

outputFolderName = paste("Output/randomTrees/",filePrefix,"/", sep = "")

# ================================================================================
# 3. DATA LOADING
# ================================================================================
if(file_ext(mainTreesLocation) == "rds"){
  if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
}else{
  if(!exists("mainTrees")){mainTrees = readTrees(mainTreesLocation)}
}

manualAnnots = read.csv(spreadSheetLocation)

# ================================================================================
# 4. PHENOTYPE MAPPING AND CLEANING 
# ================================================================================
relevantSpecies = manualAnnots[manualAnnots[[annotColumn]] %in% categoryList,]
relevantSpecies = relevantSpecies[!relevantSpecies[[nameColumn]] %in% "", ]

speciesNames = relevantSpecies[[nameColumn]] 
speciesCategories = relevantSpecies[[annotColumn]] 

phenotypeVector = speciesCategories 
names(phenotypeVector) = speciesNames 

# Run substitution loop
if(!is.null(substitutions)){
  for( i in 1:length(substitutions)){
    substitutePhenotypes = substitutions[[i]]
    message(paste("replacing", substitutePhenotypes[1], "with", substitutePhenotypes[2]))
    
    # Save and re-apply phenotype names
    saved_names <- names(phenotypeVector)
    phenotypeVector = gsub(substitutePhenotypes[1], substitutePhenotypes[2], phenotypeVector)
    names(phenotypeVector) <- saved_names
  }
}

# Remove species with no phenotype inference
invalid_or_na <- is.na(phenotypeVector) | (phenotypeVector == "") | (phenotypeVector == "NA")
if (any(invalid_or_na)) {
  cat("Found", sum(invalid_or_na), "unmapped or empty species phenotypes. Scrubbing them out...\n")
  phenotypeVector <- phenotypeVector[!invalid_or_na]
}

# ================================================================================
# 5. TREE PROCESSING AND PRUNING
# ================================================================================
workingTree = mainTrees$masterTree

# Drop species not in our clean phenotype vector
workingTree <- drop.tip(workingTree, setdiff(workingTree$tip.label, names(phenotypeVector)))

# Drop short branches
shortBranches = which(workingTree$edge.length < minimumBranchLength)
if(length(shortBranches) > 0){
  shortBranchTips = workingTree$edge[,2][shortBranches]
  shortBranchTipNames = workingTree$tip.label[shortBranchTips]
  shortBranchTipNames = shortBranchTipNames[!shortBranchTipNames %in% manualUnprunedTips]
  workingTree = drop.tip(workingTree, intersect(workingTree$tip.label, shortBranchTipNames))
}

# Finalize the alignment between vector and tree
phenotypeVector = phenotypeVector[names(phenotypeVector) %in% workingTree$tip.label]
speciesFilter = workingTree$tip.label

phenotypeTree <- workingTree

# ================================================================================
# 6. MULTINOMIAL PRUNING WITH VALIDATION (With Auto-Relaxation)
# ================================================================================
total_tips_to_keep <- 100
min_length_required <- 0.001 
max_retries <- 1000
valid_pruned_tree <- FALSE
retries <- 0

# Pull random seed from system time
set.seed(as.numeric(Sys.time()))

categories <- unique(phenotypeVector)
category_probs <- rep(1/length(categories), length(categories))

cat("Starting multinomial sampling for pruning...\n")
while(!valid_pruned_tree && retries < max_retries) {
  retries <- retries + 1
  
  # Parameter relaxation
  if (retries == 150 && min_length_required > 0.0001) {
    cat("\n!-- Struggle detected. Relaxing min_length_required to 0.0001 --!\n")
    min_length_required <- 0.0001
  }
  if (retries == 300) {
    cat("\n!-- High friction. Removing branch length constraint entirely to ensure completion --!\n")
    min_length_required <- 0
  }
  
  if(retries %% 50 == 0){
    cat("Evaluating permutation retry target:", retries, "...\n")
  }
  
  multinomial_counts <- rmultinom(1, size = total_tips_to_keep, prob = category_probs)
  row.names(multinomial_counts) <- categories
  
  species_to_keep <- unlist(lapply(categories, function(cat) {
    species_in_cat <- names(phenotypeVector[phenotypeVector == cat])
    num_to_sample <- min(multinomial_counts[cat, 1], length(species_in_cat))
    if(length(species_in_cat) > 1) {
      sample(species_in_cat, num_to_sample)
    } else {
      species_in_cat
    }
  }))
  
  # Validate
  sampled_categories <- unique(phenotypeVector[species_to_keep])
  categories_preserved <- all(categories %in% sampled_categories)
  
  if (!categories_preserved) next
  
  # Measure branch lengths
  pruned_tree <- drop.tip(phenotypeTree, setdiff(phenotypeTree$tip.label, species_to_keep))
  min_branch <- min(pruned_tree$edge.length)
  length_satisfied <- min_branch >= min_length_required
  
  # Final assessment
  if (length_satisfied) {
    prunedPaths <- tree2Paths(pruned_tree, mainTrees, binarize = T, speciesFilter)
    valid_pruned_tree <- TRUE
    message(paste("Successfully generated a valid pruned tree on attempt", retries))
  }
}

if (!valid_pruned_tree) {
  stop("Pruning limit reached: Failed to find a valid species subset.")
}

phenotypeTree <- pruned_tree
phenotypePaths <- prunedPaths

# ================================================================================
# 7. SAVE OUTPUTS
# ================================================================================
phenotypeTreeImageFilename = paste(outputFolderName, filePrefix, "phenotypeTree.pdf", sep="")
pdf(phenotypeTreeImageFilename, height = length(phenotypeVector)/7) 
plotTree(phenotypeTree)
dev.off() 

phenotypeTreeFilename = paste(outputFolderName, filePrefix, "phenotypeTree.rds", sep="")
saveRDS(phenotypeTree, phenotypeTreeFilename)

phenotypePathsFilename = paste(outputFolderName, filePrefix, "PathsFile.rds", sep= "")
saveRDS(phenotypePaths, phenotypePathsFilename)

phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") 
saveRDS(phenotypeVector, file = phenotypeVectorFilename)

message("Categorical Tree processing and sampling complete!")
