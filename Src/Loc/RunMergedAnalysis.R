# -- Usage:
# This script allows for adding a new-category in a categorical analysis which is made up of the combination of several of the main analysis' categories. 
# For example, combining "Vertivore" and "Invertivore" into "Predator", and then producing "Vertivore-Herbivore", "Invertivore-Herbivore", and "Predator-Herbivore" comparisons.
# This is done by directly combining any branches within the merged categories to the new category, and then running a second RERConverge analysis on the new merged-phenotype tree. 
# This script then adds the new pairwise comparisons produced by the merged phenotype as additional pairwise comparisons in the main analysis' output folder. 

# -- Libraries 
source("Src/Reu/loadLibraries.R")
loadLibraries() # Note that this library load includes cluster run detection. 
#-------------------------------
# --- Command arguments list ---
# - Required - 
# r = filePrefix                              REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# m = mainTreeFilename.txt or .rds            REQUIRED: This sets the location of the maintrees file
# e = c("phenotypes", "to", "merge)           This sets which phenotype should be combined into a single phenotype in a merged analysis. 
# t = c("nameOfMergedPhenotype")              This sets the name of the now merged phenotype

# - Required, Defaults likely acceptable -
# f = forceUpdate <T or F>                    This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# - Optional - 
# h = <T>                                     This ARG causes the script to return its ARG list, instead of running the script. 
# l = <min.sp value>                          This sets the min.sp value to be used in the correlation. 
# --- End arguments list ---
#-------------------------------
# -- Section to include manually-specified args --
if(!exists("args", envir = environment(), inherits = FALSE)){
  #This is a space where you can manually write in arguments for the script. 
  #The above if statement ensures these lines will not override any existing args if the script is run whole-sale. 
  
}
#-------------------------------

if(!exists("NoArgUpdate")){ #This line allows other scripts to run this code without regenerating args, if that script creates a "NoArgUpdate" variable.
  # --- Standard start-up code; identical in all scripts---
  source("Src/Reu/StandardizedStartupCode.R") #This runs the standardized code block at this location, ensuring identity across scripts. 
  #This code reads the "report args" help ARG, and reutrns the args if it is active. 
  #It then reads the filePrefix and sets the output directory to that file prefix, and reads the forceUpdate arg. 
  
  # ---------------------------
  # ----- SCRIPT MAIN CODE ----
  # ---------------------------
  
  # --- Argument Imports ---
  {   # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
    #Defaults
    mainTreesLocation = "Results/mainTreesFromAlignments.rds"
    minSpValue = 10
    toMerge = NULL
    mergeValue = NULL
    runMerge= F
    
    #Load from args 
    source("Src/Reu/ArgumentImport.R")
  }
}

#-------
# Make the merged phenotype
#-------
if(runMerge){ # This is a check that allows the script to report when it is doing nothing because no merge was supplied.
              # This is done so that it can be safely included in all categorical phenotypes, even when no merge is used. 
  
  
  #Handle the phenotyeVector
  mainPhenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
  
  mainPhenotypeVector = readRDS(mainPhenotypeVectorFilename)
  mainPhenotypeLabels = names(table(mainPhenotypeVector))
  
  replaceIndexes = which(mainPhenotypeLabels %in% toMerge)
  
  
  mergedPhenotypeLabels = mainPhenotypeLabels[-replaceIndexes]
  mergedPhenotypeLabels = append(mergedPhenotypeLabels, mergeValue)
  mergedPhenotypeLabels = mergedPhenotypeLabels[order(mergedPhenotypeLabels)]
  
  mergedPhenotypeVector = mainPhenotypeVector
  for(i in 1:length(toMerge)){
    mergedPhenotypeVector = gsub(toMerge[i], mergeValue, mergedPhenotypeVector)
  }
  
  mergedPhenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalMergedPhenotypeVector.rds",sep="") #make a filename based on the prefix
  saveRDS(mergedPhenotypeVector, mergedPhenotypeVectorFilename)
  
  #Handle the phenotype tree
  mainCategoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
  mainPhenotypeTree = readRDS(mainCategoricalTreeFilename)
  
  bufferValue = 10
  mergedTree = mainPhenotypeTree
  for(i in 1:length(mainPhenotypeLabels)){
    if(mainPhenotypeLabels[i] %in% toMerge){
      mergedPosition = which(mergedPhenotypeLabels %in% mergeValue)
      mergedTree$edge.length[mergedTree$edge.length == i] = mergedPosition+bufferValue
    }else{
      currentLabel = mainPhenotypeLabels[i]
      mergedPosition = which(mergedPhenotypeLabels == currentLabel)
      mergedTree$edge.length[mergedTree$edge.length == i] = mergedPosition+bufferValue
    }
  }
  mergedTree$edge.length = mergedTree$edge.length-bufferValue
  
  mergedCategoricalTreeFilename = paste(outputFolderName, filePrefix, "MergedCategoricalTree.rds", sep="") #make a filename based on the prefix
  saveRDS(mergedTree, mergedCategoricalTreeFilename)
  
  #-------
  # Run RER and correlation for it 
  #-------
  phenotypeTreeLocation = mergedCategoricalTreeFilename
  speciesFilterFileName = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #Make the name of the location a pre-made filter would have to test for it
  speciesFilter = readRDS(speciesFilterFileName)                       #if so, use it 
  phenotypeStyle = "Categorical"
  correlationFileNameValue = paste(outputFolderName, filePrefix, "MergedCorrelationFile", sep= "") #Make a correlation filename based on the prefix  
  phenotypeVectorFilenameValue = mergedPhenotypeVectorFilename 
  pairwiseCorrelationFileNameValue = paste(outputFolderName, filePrefix, "MergedPairwiseCorrelationFile", sep= "") #make a name for the pairwise comparisons based on prefix
  combinedCategoricalCorrelationFilenameValue = paste(outputFolderName, filePrefix, "MergedCombinedCategoricalCorrelationFile", sep= "") # make this file for later functions that want it in combo
  correlationsOverallFilenameValue = paste(outputSubdirectory, filePrefix, "MergedOverallCorrelationFile.rds", sep= "")
  
  #Run RER and Correlation script while block arg import from primary args. 
  alreadyBlockingArgImports = exists("NoArgUpdate")
  if(!alreadyBlockingArgImports){NoArgUpdate = "Exists"}
  
  source("Src/Loc/RunRERAndCorrelation.R")
  
  if(!alreadyBlockingArgImports){rm(NoArgUpdate)}
  
  #-------
  # Combine the existing Correlation file with the merged one 
  #-------
  
  #Figure out which folders are new 
  folderNames = basename(list.dirs(outputFolderName, full.names = TRUE, recursive = FALSE))  
  pairwiseFolders = folderNames[grep("-", folderNames)]
  mergeUniquePairwiseNames = pairwiseFolders[grep(mergeValue, pairwiseFolders)]
  
  #Add them to the combined correlation file 
  mainPairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "PairwiseCorrelationFile", sep= "") #make a name for the pairwise comparisons based on prefix
  mainPairwiseCategorical = readRDS(paste0(mainPairwiseCorrelationFileName, ".rds"))
  mergedPairwiseCategorical = readRDS(pairwiseCorrelationFileNameValue)
  
  names(mainPairwiseCategorical) = gsub(" ", "", names(mainPairwiseCategorical))
  names(mergedPairwiseCategorical) = gsub(" ", "", names(mergedPairwiseCategorical))
  
  for(i in 1:length(mergeUniquePairwiseNames)){
    currentName = mergeUniquePairwiseNames[i]
    
    if(currentName %in% names(mainPairwiseCategorical)){
      currentIndex = which(names(mainPairwiseCategorical) %in% currentName)
      mainPairwiseCategorical[currentIndex] = mergedPairwiseCategorical[which(names(mergedPairwiseCategorical) %in% currentName)]
    }else{
      mainPairwiseCategorical = append(mainPairwiseCategorical, mergedPairwiseCategorical[which(names(mergedPairwiseCategorical) %in% currentName)])
    }
    mergedPairwiseCategorical
  }
  
  write.csv(mainPairwiseCategorical, file= paste(mainPairwiseCorrelationFileName, ".csv", sep=""), row.names = T, quote = F) #save the correlations as a csv
  saveRDS(mainPairwiseCategorical, paste(mainPairwiseCorrelationFileName, ".rds", sep="")) #and as an rds 
}