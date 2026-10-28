####################################################################################################################
####    Initial script written by Michael Tene and Alexander Seaver                                             ####
####    Select sections have been revised utilizing Google Gemini 3.6, manually reviewed by Alexander Seaver    ####
####################################################################################################################

# ==============================================================================
# DIAGNOSTIC / TROUBLESHOOTING VERSION OF hyphyFileMake.R
# ==============================================================================

cat("\n======================================================================\n")
cat(" [DEBUG 0] INITIALIZATION & ENVIRONMENT CHECK\n")
cat("======================================================================\n")

clusterRun = F
clusterRun = T

# MUST assign args so cmdArgImport reads commandArgs instead of the R base args() function
if(clusterRun){
  args = commandArgs(trailingOnly = TRUE)
  libPath = "/share/ceph/wym219group/shared/libraries/R4"
  if(dir.exists(libPath)){
    .libPaths(libPath)
  }
}

cat("[DEBUG 0] Working Directory:", getwd(), "\n")
cat("[DEBUG 0] R Version:", R.version.string, "\n")

# Load external libraries with failure reporting
load_pkg = function(pkg) {
  cat("[DEBUG 0] Loading package:", pkg, "... ")
  suppressPackageStartupMessages({
    loaded = require(pkg, character.only = TRUE)
  })
  if(loaded) {
    cat("OK (v", as.character(packageVersion(pkg)), ")\n", sep="")
  } else {
    cat("FAILED!\n")
    stop(paste("Failed to load required package:", pkg))
  }
}

load_pkg("seqinr")
load_pkg("RERconverge")

# Source required dependency files with existence verification
source_file = function(filePath) {
  cat("[DEBUG 0] Sourcing script:", filePath, "... ")
  if(file.exists(filePath)) {
    source(filePath)
    cat("OK\n")
  } else {
    cat("FILE NOT FOUND!\n")
    stop(paste("Required source file missing:", filePath))
  }
}

source_file("Src/Reu/cmdArgImport.R")
source_file("Src/Reu/paths2Tree.R")
source_file("Src/Reu/customSeqinrFunctions.R")
source_file("Src/Reu/makePhenMasterTree.R")

cat("\n======================================================================\n")
cat(" [DEBUG 1] COMMAND LINE ARGUMENTS ANALYSIS\n")
cat("======================================================================\n")

rawArgs = commandArgs(trailingOnly = TRUE)
cat("[DEBUG 1] Raw commandArgs(trailingOnly = TRUE):\n")
print(rawArgs)

if(length(rawArgs) == 0 && clusterRun) {
  cat("[WARN 1] No commandArgs received despite clusterRun = TRUE.\n")
}

# Helper to safely test and report imported command args
import_and_report = function(argChar, argDescription) {
  val = cmdArgImport(argChar)
  cat(sprintf("[DEBUG 1] Arg [-%s] (%s): %s\n", argChar, argDescription, ifelse(is.na(val), "NA (Not provided)", paste0("'", val, "'"))))
  return(val)
}

filePrefixVal            = import_and_report('r', "filePrefix")
forceUpdateVal          = import_and_report('v', "forceUpdate")
mainTreesVal            = import_and_report('m', "mainTreesLocation")
geneNameVal             = import_and_report('g', "geneName")
phenotypeTreeLocVal     = import_and_report('p', "phenotypeTreeLocation")
fastaLocVal             = import_and_report('a', "fastaLocation")
foregroundCategoryVal   = import_and_report('f', "foregroundCategory")
binarizeTreeVal         = import_and_report('b', "binarizeTree")

# Validate core required parameters
if(is.na(filePrefixVal)){
  stop("[FATAL 1] Required argument -r (filePrefix) is missing or NA.")
}
if(is.na(geneNameVal)){
  stop("[FATAL 1] Required argument -g (geneName) is missing or NA.")
}
if(is.na(fastaLocVal)){
  stop("[FATAL 1] Required argument -a (fastaLocation) is missing or NA.")
}

filePrefix = filePrefixVal
geneName = geneNameVal
fastaLocation = fastaLocVal

# Set output directories
outputFolderNameNoSlash = paste0("Output/", filePrefix)
outputFolderName = paste0(outputFolderNameNoSlash, "/")
outputHyphyFolderName = paste0(outputFolderName, "Hyphy/")

cat("[DEBUG 1] Output directory target:", outputFolderName, "\n")
cat("[DEBUG 1] HyPhy output directory target:", outputHyphyFolderName, "\n")

if(!dir.exists("Output")) dir.create("Output", recursive = TRUE)
if(!dir.exists(outputFolderNameNoSlash)) dir.create(outputFolderNameNoSlash, recursive = TRUE)
if(!dir.exists(outputHyphyFolderName)) dir.create(outputHyphyFolderName, recursive = TRUE)

# Parse optional arguments
if(!is.na(mainTreesVal)) {
  mainTreesLocation = mainTreesVal
} else {
  if(!clusterRun) {
    mainTreesLocation = '../RunRERBinaryMT/Data/zoonomiaAllMammalsTrees.rds'
  } else {
    mainTreesLocation = 'Data/zoonomiaAllMammalsTrees.rds'
  }
}
cat("[DEBUG 1] Resolved mainTreesLocation:", mainTreesLocation, "\n")
cat("[DEBUG 1] Checking if mainTreesLocation file exists:", file.exists(mainTreesLocation), "\n")

useManualTree = FALSE
phenotypeTreeLocation = NULL
if(!is.na(phenotypeTreeLocVal)) {
  phenotypeTreeLocation = phenotypeTreeLocVal
  useManualTree = TRUE
  cat("[DEBUG 1] Manual phenotype tree specified. Checking existence:", file.exists(phenotypeTreeLocation), "\n")
}

fgCat = ifelse(is.na(foregroundCategoryVal), "dshakgldskagkjshadgkhalgh", foregroundCategoryVal)
binTree = ifelse(is.na(binarizeTreeVal), FALSE, as.logical(binarizeTreeVal))

cat(sprintf("[DEBUG 1] foregroundCategory: '%s' | binarizeTree: %s\n", fgCat, as.character(binTree)))
foregroundCategoryValue = fgCat

cat("\n======================================================================\n")
cat(" [DEBUG 2] EXECUTING makePhenMasterTree()\n")
cat("======================================================================\n")

phenMasterTree = tryCatch({
  makePhenMasterTree(
    geneName = geneName,
    filePrefix = filePrefix,
    manualPhenotypeTreeLocation = phenotypeTreeLocation,
    foregroundCategory = fgCat,
    binarizeTree = binTree
  )
}, error = function(e) {
  cat("[FATAL 2] makePhenMasterTree() execution crashed!\n")
  cat("[FATAL 2] Error message:", conditionMessage(e), "\n")
  stop(e)
})

cat("[DEBUG 2] phenMasterTree object structure check:\n")
cat("[DEBUG 2]  - Class:", class(phenMasterTree), "\n")
cat("[DEBUG 2]  - Tip Count:", length(phenMasterTree$tip.label), "\n")
cat("[DEBUG 2]  - Node Count:", phenMasterTree$Nnode, "\n")

if(exists("originalTipValues", envir = .GlobalEnv)) {
  cat("[DEBUG 2] 'originalTipValues' found in GlobalEnv. Count:", length(originalTipValues), "\n")
  cat("[DEBUG 2] Sample originalTipValues (first 5):", paste(head(originalTipValues, 5), collapse=", "), "\n")
} else {
  cat("[WARN 2] 'originalTipValues' not found in GlobalEnv. Assigning from phenMasterTree$tip.label.\n")
  originalTipValues = phenMasterTree$tip.label
}

cat("\n======================================================================\n")
cat(" [DEBUG 3] READING & INSPECTING FASTA FILE\n")
cat("======================================================================\n")

cat("[DEBUG 3] FASTA Path:", fastaLocation, "\n")

if(!file.exists(fastaLocation)) {
  stop(paste("[FATAL 3] Target FASTA file does NOT exist at path:", fastaLocation))
}

fastaFileInfo = file.info(fastaLocation)
cat("[DEBUG 3] FASTA File Size (bytes):", fastaFileInfo$size, "\n")
if(fastaFileInfo$size == 0) {
  stop("[FATAL 3] Target FASTA file is 0 bytes (empty)!")
}

fasta = tryCatch({
  read.fasta(fastaLocation)
}, error = function(e) {
  cat("[FATAL 3] seqinr::read.fasta() failed to parse file!\n")
  cat("[FATAL 3] Error message:", conditionMessage(e), "\n")
  stop(e)
})

cat("[DEBUG 3] Raw FASTA sequence count loaded:", length(fasta), "\n")

lowQualityAlignment = grep("!", fasta)
cat("[DEBUG 3] Low quality alignment entries found (!):", length(lowQualityAlignment), "\n")
if(length(lowQualityAlignment) > 0) {
  fasta = fasta[-lowQualityAlignment]
  cat("[DEBUG 3] Sequence count after filtering '!':", length(fasta), "\n")
}

fastaTipHeaders = names(fasta)
cat("[DEBUG 3] Sample raw FASTA header names (first 5):\n")
print(head(fastaTipHeaders, 5))

# Sanitize headers
fastaTipHeaders = sub("\\t.*", "", fastaTipHeaders)
fastaTipHeaders[fastaTipHeaders == "REFERENCE"] = "vs_hg38"
names(fasta) = fastaTipHeaders

cat("\n======================================================================\n")
cat(" [DEBUG 4] FASTA VS TREE TIP OVERLAP ANALYSIS\n")
cat("======================================================================\n")

matchingSpecies = intersect(fastaTipHeaders, originalTipValues)
cat("[DEBUG 4] Number of matching species/tips between FASTA and Tree:", length(matchingSpecies), "\n")

if(length(matchingSpecies) == 0) {
  cat("[FATAL 4] ZERO overlap between FASTA headers and Tree tips!\n")
  cat("[DEBUG 4] Sample FASTA headers:", paste(head(fastaTipHeaders, 5), collapse=" | "), "\n")
  cat("[DEBUG 4] Sample Tree tips:", paste(head(originalTipValues, 5), collapse=" | "), "\n")
  stop("FATAL ERROR: Mismatch between FASTA headers and Tree tip names.")
}

noDataTips = phenMasterTree$tip.label[!originalTipValues %in% fastaTipHeaders]
cat("[DEBUG 4] Tree tips missing from FASTA (noDataTips count):", length(noDataTips), "\n")

if(length(noDataTips) > 0 && length(noDataTips) < length(phenMasterTree$tip.label)) {
  cat("[DEBUG 4] Dropping missing tips from tree...\n")
  phenMasterTree = drop.tip(phenMasterTree, noDataTips)
  cat("[DEBUG 4] Remaining tree tips count:", length(phenMasterTree$tip.label), "\n")
} else if(length(noDataTips) == length(phenMasterTree$tip.label)) {
  stop("[FATAL 4] Attempting to drop ALL tips from phenMasterTree!")
}

fastaToDrop = which(!fastaTipHeaders %in% originalTipValues)
cat("[DEBUG 4] FASTA sequences missing from Tree (fastaToDrop count):", length(fastaToDrop), "\n")

# GUARD: Check length > 0 to prevent fasta[-integer(0)] from wiping the list
if(length(fastaToDrop) > 0 && length(fastaToDrop) < length(fasta)) {
  cat("[DEBUG 4] Dropping non-tree sequences from FASTA...\n")
  fasta = fasta[-fastaToDrop]
  fastaTipHeaders = fastaTipHeaders[-fastaToDrop]
  names(fasta) = fastaTipHeaders
  cat("[DEBUG 4] Remaining FASTA sequence count:", length(fasta), "\n")
} else if(length(fastaToDrop) == length(fasta)) {
  stop("[FATAL 4] Attempting to drop ALL sequences from FASTA dataset!")
}

cat("\n======================================================================\n")
cat(" [DEBUG 5] FASTA & TREE SERIALIZATION\n")
cat("======================================================================\n")

if(!dir.exists("Results")) dir.create("Results", recursive = TRUE)

fastaLinesVar = NULL
cat("[DEBUG 5] Calling write.fastaToVar()...\n")
tryCatch({
  write.fastaToVar(fasta, names = names(fasta), file.out = "fastaLinesVar")
}, error = function(e) {
  cat("[WARN 5] write.fastaToVar execution error:", conditionMessage(e), "\n")
})

# Fallback generator if write.fastaToVar failed or did not set GlobalEnv variable
if(!exists("fastaLinesVar") || is.null(fastaLinesVar) || length(fastaLinesVar) == 0) {
  cat("[WARN 5] fastaLinesVar missing/null after write.fastaToVar. Generating manual FASTA strings...\n")
  fastaLinesVar = c()
  for(i in seq_along(fasta)) {
    fastaLinesVar = c(fastaLinesVar, paste0(">", names(fasta)[i]), paste(fasta[[i]], collapse = ""))
  }
}

cat("[DEBUG 5] 'fastaLinesVar' line count:", length(fastaLinesVar), "\n")

treeOut = write.tree(phenMasterTree)
cat("[DEBUG 5] 'treeOut' output string length:", nchar(treeOut), "\n")
cat("[DEBUG 5] Sample treeOut:", substr(treeOut, 1, 100), "...\n")

combinedContent = c(fastaLinesVar, treeOut)
cat("[DEBUG 5] 'combinedContent' total lines/elements:", length(combinedContent), "\n")

if(length(combinedContent) == 0 || all(nchar(combinedContent) == 0)) {
  stop("[FATAL 5] Final combinedContent is empty! Cannot write output file.")
}

cat("\n======================================================================\n")
cat(" [DEBUG 6] WRITING OUTPUT FILES\n")
cat("======================================================================\n")

fastaOutputFilename = paste0(outputHyphyFolderName, filePrefix, geneName, "_foreground_", foregroundCategoryValue, "_HyphyInputFile.fna")
tempOutputFilename = "Results/TempHyphyInputFile.fna"

cat("[DEBUG 6] Target output path 1:", fastaOutputFilename, "\n")
cat("[DEBUG 6] Target output path 2:", tempOutputFilename, "\n")

writeLines(combinedContent, fastaOutputFilename)
writeLines(combinedContent, tempOutputFilename)

# Verification checks
outInfo1 = file.info(fastaOutputFilename)
outInfo2 = file.info(tempOutputFilename)

cat("[DEBUG 6] Output file 1 size (bytes):", outInfo1$size, "\n")
cat("[DEBUG 6] Output file 2 size (bytes):", outInfo2$size, "\n")

if(is.na(outInfo1$size) || outInfo1$size == 0) {
  stop("[FATAL 6] Primary output file was created, but is 0 bytes!")
}

cat("\n======================================================================\n")
cat(" [SUCCESS] hyphyFileMake.R completed successfully for gene:", geneName, "\n")
cat("======================================================================\n")
