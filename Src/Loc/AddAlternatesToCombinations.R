# -- Usage:
# This script is used to Add alternate information to the combined results from the main analysis. 
# It uses a ratio (perentage of total altenrates) in which a gene or category must be singificant to be considered robust to alternate species set selection. 

# -- Libraries 
source("Src/Reu/loadLibraries.R")
loadLibraries() # Note that this library load includes cluster run detection. 

source("Src/Reu/trimGmtFilenameToSetname.R")

#-------------------------------
# --- Command arguments list ---
# - Required - 
# r = filePrefix                              REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# g = c("gene", "Sets")                       This is the gmtset(as shown in the filename), or the gmt file path of enrichment runs to be combined. 


# - Required, Defaults likely acceptable -
# f = forceUpdate <T or F>                    This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# o = alternateRatio <decimal>                This sets the percentage of alternates which must be significant for a result to be considered robust to alternate species selection. 
# w = writeData  <T or F>                     This sets if the data should be saved, potentially overwriting old combinations

# - Optional - 
# h = <T>                                     This ARG causes the script to return its ARG list, instead of running the script. 
# c = significance cutoff                     This sets the cutoff for inclusion in the significant columns
# - Overrides - 
# --- End arguments list ---
#-------------------------------
# -- Section to include manually-specified args --
if(!exists("args", envir = environment(), inherits = FALSE)){
  #This is a space where you can manually write in arguments for the script. 
  #The above if statement ensures these lines will not override any existing args if the script is run whole-sale. 
  
}
#-------------------------------

# --- Standard start-up code; identical in all scripts---
source("Src/Reu/StandardizedStartupCode.R") #This runs the standardized code block at this location, ensuring identity across scripts. 
#This code reads the "report args" help ARG, and reutrns the args if it is active. 
#It then reads the filePrefix and sets the output directory to that file prefix, and reads the forceUpdate arg. 



# --- Argument Imports ---
{  # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
  CurrentScript = "AddAlternatesToCombinations.R"
  #Defaults
  geneSets = NULL
  usingGo = F
  usingGene = T
  validGeneArguments = c("NULL", "gene", "Gene")
  alternateCutoffRatio = 0.5 
  writeData = T
  significanceCutoff = 0.05
  
  #Load from args 
  source("Src/Reu/ArgumentImport.R")
  
  geneSets = trimGmtFilenameToSetname(geneSets)
}




#----- 
#Main code

#Determine the number of alternates required to exceed the required ratio 
PairwiseCorrelationsFilenames = list.files(paste0(outputFolderName, "Alternates/"), pattern = "PairwiseCorrelationFile.rds")
PairwiseCorrelationsFilenames = PairwiseCorrelationsFilenames[-grep("Merged", PairwiseCorrelationsFilenames)]
totalAlternates = length(PairwiseCorrelationsFilenames)
alternateCutoff = totalAlternates * alternateCutoffRatio




for(j in 1:length(geneSets)){
  geneSet = geneSets[j]
  if(geneSet %in% validGeneArguments){
    isGene = T
  }else{
    isGene = F
  }

  
  if(isGene){
    #Read in the files
    combinedDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResults.rds")
    combinedData = readRDS(combinedDataFilename)
    
    alternatesFilename = paste0(outputFolderName, filePrefix, "AlternatesCombinedCorrelations.rds")
    alternatesData = readRDS(alternatesFilename)
    
  }else{
    combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResults-", geneSet, ".rds")
    combinedData = readRDS(combinedGODataFilename)
    
    alternatesFilename = paste0(outputFolderName, filePrefix, "AlternatesCombinedEnrichments-" , geneSet, ".rds")
    alternatesData = readRDS(alternatesFilename)
    
    
  }
  
  
  #Reformat the alternates   
  alternateSummaryColumns = map(alternatesData, ~ .x[, 4:ncol(.x), drop = FALSE])
  
  alternatesSingleDF <- imap(alternateSummaryColumns, function(df, nm) {
    parts <- strsplit(nm, "-")[[1]]
    prefix <- paste0(substr(parts[1], 1, 1),
                     substr(parts[2], 1, 1))
    
    df %>% rename_with(~ paste0(prefix, "-", .x))
  }) %>%
    bind_cols()
  
  
  if(!isGene){
    rownames(alternatesSingleDF) = names(alternatesData[[1]][[1]])
    
    #ensure that the GO are in the same order
    #Ensure that the GO categories are in the same order 
    longestGOLength = max(length(rownames(combinedData)), length(rownames(alternatesSingleDF)))
    if(length(rownames(combinedData)) == longestGOLength){
      newGoOrder = rownames(combinedData)
    }else{
      newGoOrder = rownames(alternatesSingleDF)
    }
    
    
    if(!all(rownames(combinedData) == newGoOrder)){
      message("Main data order not match longest GO order, reordering main data")
      matchedOrder = match(newGoOrder, rownames(combinedData))
      combinedData = combinedData[matchedOrder, ]
      rownames(combinedData) = newGoOrder
    }
    
    if(!all(rownames(alternatesSingleDF) == newGoOrder)){
      message("Alternate data order not match longest GO order, reordering Alternate data")
      matchedOrder = match(newGoOrder, rownames(alternatesSingleDF))
      alternatesSingleDF = alternatesSingleDF[matchedOrder, ]
      rownames(alternatesSingleDF) = newGoOrder
    }
    
  }
  
  #Add the alternate data to the main dataframe 
  combinedData = cbind(combinedData, alternatesSingleDF)
  
  
  
  
  {  #Reorder the columns
    
    cols <- colnames(combinedData)
    
    # get all prefixes (e.g., CH, CO, HI, etc.)
    prefixes <- unique(sub("-.*", "", cols))
    
    # function to order columns for one prefix
    
    if(isGene){
      order_prefix <- function(p) {
        c(
          paste0(p, "-significant"),
          paste0(p, "-Rho"),
          paste0(p, "-RhoMedian"),
          paste0(p, "-RhoMean"),
          paste0(p, "-RhoSD"),
          paste0(p, "-RhoCI95"),
          paste0(p, "-RhoIQR"),
          paste0(p, "-RhoMaxDiff"),
          paste0(p, "-P"),
          paste0(p, "-PNumSignificant"),
          paste0(p, "-PMedian"),
          paste0(p, "-PMean"),
          paste0(p, "-PSD"),
          paste0(p, "-PCI95"),
          paste0(p, "-PIQR"),
          paste0(p, "-PMaxDiff"),
          paste0(p, "-p.adj"),
          paste0(p, "-PadjNumSignificant"),
          paste0(p, "-PadjMedian"),
          paste0(p, "-PadjMean"),
          paste0(p, "-PadjSD"),
          paste0(p, "-PadjCI95"),
          paste0(p, "-PadjIQR"),
          paste0(p, "-PadjMaxDiff")
        )
      }
    }else{
      order_prefix <- function(p) {
        c(
          paste0(p, "-significant"),
          paste0(p, "-num.genes"),
          paste0(p, "-gene.vals"),
          paste0(p, "-stat"),
          paste0(p, "-statMedian"),
          paste0(p, "-statMean"),
          paste0(p, "-statSD"),
          paste0(p, "-statCI95"),
          paste0(p, "-statIQR"),
          paste0(p, "-statMaxDiff"),
          paste0(p, "-pval"),
          paste0(p, "-PNumSignificant"),
          paste0(p, "-PMedian"),
          paste0(p, "-PMean"),
          paste0(p, "-PSD"),
          paste0(p, "-PCI95"),
          paste0(p, "-PIQR"),
          paste0(p, "-PMaxDiff"),
          paste0(p, "-p.adj"),
          paste0(p, "-PadjNumSignificant"),
          paste0(p, "-PadjMedian"),
          paste0(p, "-PadjMean"),
          paste0(p, "-PadjSD"),
          paste0(p, "-PadjCI95"),
          paste0(p, "-PadjIQR"),
          paste0(p, "-PadjMaxDiff")
        )
      }
    }
    
    
    # build desired order for all prefixes
    ordered_main <- unlist(lapply(prefixes, order_prefix))
    
    # keep only columns that actually exist (important!)
    ordered_main <- ordered_main[ordered_main %in% cols]
    
    # everything else (overlaps, deltas, etc.)
    remaining <- setdiff(cols, ordered_main)
    
    # final order
    #combinedData <- combinedData[, c(ordered_main, remaining)]
    combinedData <- combinedData[, c(ordered_main)]
  }
  
  
  # --- Add column which includes the main analysis as an alternate
  
  prefixes <- unique(sub("-.*", "", names(combinedData)))
  
  for (prefix in prefixes) {
    
    p_col <- paste0(prefix, "-P")
    padj_col <- paste0(prefix, "-p.adj")
    
    pnum_col <- paste0(prefix, "-PNumSignificant")
    padjnum_col <- paste0(prefix, "-PadjNumSignificant")
    
    # skip if columns don't exist
    if (!all(c(p_col, padj_col, pnum_col, padjnum_col) %in% names(combinedData))) next
    
    # new column names
    pnum_new <- paste0(prefix, "-PNumSignificantInclusive")
    padjnum_new <- paste0(prefix, "-PadjNumSignificantInclusive")
    
    # create new values
    combinedData[[pnum_new]] <- combinedData[[pnum_col]] + 
      ifelse(combinedData[[p_col]] < 0.05, 1, 0)
    
    combinedData[[padjnum_new]] <- combinedData[[padjnum_col]] + 
      ifelse(combinedData[[padj_col]] < 0.05, 1, 0)
    
    # move columns to correct position (after originals)
    p_index <- match(pnum_col, names(combinedData))
    padj_index <- match(padjnum_col, names(combinedData))
    
    # reorder for PNum
    combinedData <- combinedData[, append(
      names(combinedData)[-which(names(combinedData) == pnum_new)],
      pnum_new,
      after = p_index
    )]
    
    # reorder for PadjNum
    combinedData <- combinedData[, append(
      names(combinedData)[-which(names(combinedData) == padjnum_new)],
      padjnum_new,
      after = padj_index
    )]
  }
  
  
  
  for (prefix in prefixes) {
    
    sigCol <- paste0(prefix, "-significant")
    
    padjnum_col <- paste0(prefix, "-PadjNumSignificant")
    
    # skip if columns don't exist
    if (!all(c(sigCol,  padjnum_col) %in% names(combinedData))) next
    
    # new column names
    newCol <- paste0(prefix, "-significantRobust")
    
    
    # create new values
    combinedData[[newCol]] = combinedData[[sigCol]] & combinedData[[padjnum_col]] > alternateCutoff
    
    
    # move columns to correct position (after originals)
    sigCol_index <- match(sigCol, names(combinedData))
    
    # reorder for PNum
    combinedData <- combinedData[, append(
      names(combinedData)[-which(names(combinedData) == newCol)],
      newCol,
      after = sigCol_index
    )]
    
  }
  
  
  
  if(writeData){
    
    if(isGene){
      combinedDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResultsWithAlternates")
      write.csv(combinedData, paste0(combinedDataFilename, ".csv"))
      saveRDS(combinedData, paste0(combinedDataFilename, ".rds"))    
    }else{
      combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResultsWithAlternates-", geneSet)
      write.csv(combinedData, paste0(combinedGODataFilename, ".csv"))
      saveRDS(combinedData, paste0(combinedGODataFilename, ".rds"))
    }
    
  }
  
  
  
}









