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
# v = forceUpdate <T or F>                    This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
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


# --- Standard start-up code; identical in all scripts---
source("Src/Reu/StandardizedStartupCode.R") #This runs the standardized code block at this location, ensuring identity across scripts. 
#This code reads the "report args" help ARG, and reutrns the args if it is active. 
#It then reads the filePrefix and sets the output directory to that file prefix, and reads the forceUpdate arg. 

# ---------------------------
# ----- SCRIPT MAIN CODE ----
# ---------------------------

# --- Argument Imports ---
{   # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
  # Defaults
  mainTreesLocation = "data/zoonomiaAllMammalsTrees.rds"  #this is the standard location on Sol 
  minSpValue = 10
  toMerge = NULL
  mergeValue = NULL
  runMerge= F
  
  #MainTrees Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using default")
  }
  #min.sp Value
  if(!is.na(cmdArgImport('l'))){
    minSpValue = cmdArgImport('l')
  }else{
    message("min.sp not specified, using 10")
  }
  
  #toMerge
  if(!all(is.na(cmdArgImport('e')))){
    toMerge = cmdArgImport('e')
    runMerge = T
  }else{
    message("No phenotypes set to merge. Script will effectively do nothing; if you are not running a merged phenotype, this is fine.")
    runMerge = F
  }
  
  #mergeValue
  if(!is.na(cmdArgImport('t'))){
    mergeValue = cmdArgImport('t')
  }else{
    if(!runMerge){
      message("No phenotypes set to merge. Script will effectively do nothing; if you are not running a merged phenotype, this is fine.")  
    }else{
      stop("ERROR: Attempting to run merge, but the merged phenotype is unnamed. Exiting.")
    }
  }
}


#-------
# Make the merged phenotype
#-------
if(runMerge){ # This is a check that allows the script to reportwhen it is doing nothing because no merge was supplied.
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
  
  
  
  # -- Read in the trees --
  #MainTrees
  if(file_ext(mainTreesLocation) == "rds"){
    if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
  }else{
    if(!exists("mainTrees")){mainTrees = readTrees(mainTreesLocation)} 
  }
  #Phenotype tree
  phenotypeTree = mergedTree
  
  
  #Species filter
  speciesFilterFileName = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #Make the name of the location a pre-made filter would have to test for it
  speciesFilter = readRDS(speciesFilterFileName)                       #if so, use it 
  
  
  # --- RERs ---
  
  RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
  
  if(!file.exists(paste(RERFileName)) | forceUpdate){                             #if it does not exist, or update is forced 
    gc()
    RERObject = getAllResiduals(mainTrees, useSpecies = speciesFilter, plot = F)  #Calculate the RERs
    saveRDS(RERObject, file = RERFileName)                                        #Save them
  }else{                                                                          #Otherwise
    RERObject = readRDS(RERFileName)                                              #Use the existing ones
  }
  
  # --- PATHS ---
  
  pathsObject = tree2Paths(phenotypeTree, mainTrees, useSpecies = speciesFilter, categorical = TRUE) #do not binarize; the categorical data is already contained in the phenotype tree.
  
  
  
  # --- CORRELATION ---
  correlationFileName = paste(outputFolderName, filePrefix, "MergedCorrelationFile", sep= "") #Make a correlation filename based on the prefix
  
  {                                     #if categorical
    categoricalCorrelation = correlateWithCategoricalPhenotype(RERObject, pathsObject, min.sp = minSpValue, min.pos = 2) #Calculate with categorical, min 2 species per category 
    overalCategorical = categoricalCorrelation[[1]]                               #select the results relating to overall difference between all categories
    correlation = overalCategorical                                               # and classify it as the main correlation file
    
    #process the pairwise outputs
    pairwiseCategorical = categoricalCorrelation[[2]]                             #select the group of pairwise comparisons
    
    
    phenotypeVector = mergedPhenotypeVector                            #load in the phenotype vector 
    categories = map_to_state_space(phenotypeVector)                              #and use it to connect branch lengths to phenotype name
    categoryNames = categories$name2index                                         #store the length-phenotype connection
    
    pairwiseTableNames = names(pairwiseCategorical)                               #Prepare to repalce the number-number titles with phenotype-phenotype titles
    for(i in 1:length(categoryNames)){                                            #for each phenotype
      pairwiseTableNames= gsub(i, names(categoryNames)[i], pairwiseTableNames)                        #replace the number with the phenotype name  
    }
    names(pairwiseCategorical) = pairwiseTableNames                               #update the dataframe titles
    
    pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "MergedPairwiseCorrelationFile", sep= "") #make a name for the pairwise comparisons based on prefix
    write.csv(pairwiseCategorical, file= paste(pairwiseCorrelationFileName, ".csv", sep=""), row.names = T, quote = F) #save the correlations as a csv
    saveRDS(pairwiseCategorical, paste(pairwiseCorrelationFileName, ".rds", sep="")) #and as an rds 
    
    combinedCategoricalCorrelationFilename = pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "MergedCombinedCategoricalCorrelationFile", sep= "") # make this file for later functions that want it in combo
    saveRDS(categoricalCorrelation, paste(combinedCategoricalCorrelationFilename, ".rds", sep="")) #and as an rds 
    
    #save the outputs to subdirectories 
    outputSubdirectoryNoslash = paste(outputFolderName, "Overall", sep = "")
    if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
      dir.create(outputSubdirectoryNoslash)
    }
    outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
    
    correlationsOverallFilename = paste(outputSubdirectory, filePrefix, "MergedOverallCorrelationFile.rds", sep= "")
    saveRDS(categoricalCorrelation[[1]], correlationsOverallFilename)
    
    mergeUniquePairwiseNames = character()
    
    for(i in 1:length(pairwiseTableNames)){
      pairwiseTableNames= gsub(" ", "", pairwiseTableNames)
      
      outputSubdirectoryNoslash = paste(outputFolderName, pairwiseTableNames[i], sep = "")
      if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
          dir.create(outputSubdirectoryNoslash)
      }
      if(length(grep(mergeValue, outputSubdirectoryNoslash)) != 0 ){ #check that this is a merge-unique directory
        outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
        
        correlationsPairFilename = paste(outputSubdirectory, filePrefix, pairwiseTableNames[i], "CorrelationFile",".rds", sep= "")
        saveRDS(categoricalCorrelation[[2]][[i]], correlationsPairFilename)
        mergeUniquePairwiseNames = append(mergeUniquePairwiseNames, pairwiseTableNames[i])
      }
    }
  
  }
  
  write.csv(correlation, file= paste(correlationFileName, ".csv", sep=""), row.names = T, quote = F) #Save correlations as csv
  saveRDS(correlation, paste(correlationFileName, ".rds", sep=""))                          #and as an rds 
  
  
  #-------
  # Combine the existing Correlation file with the merged one 
  #-------
  
  mainPairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "PairwiseCorrelationFile", sep= "") #make a name for the pairwise comparisons based on prefix
  mainPairwiseCategorical = readRDS(paste0(mainPairwiseCorrelationFileName, ".rds"))
  
  names(mainPairwiseCategorical) = gsub(" ", "", names(mainPairwiseCategorical))
  names(pairwiseCategorical) = gsub(" ", "", names(pairwiseCategorical))
  
  for(i in 1:length(mergeUniquePairwiseNames)){
    currentName = mergeUniquePairwiseNames[i]
    
    if(currentName %in% names(mainPairwiseCategorical)){
      currentIndex = which(names(mainPairwiseCategorical) %in% currentName)
      mainPairwiseCategorical[currentIndex] = pairwiseCategorical[which(names(pairwiseCategorical) %in% currentName)]
    }else{
      mainPairwiseCategorical = append(mainPairwiseCategorical, pairwiseCategorical[which(names(pairwiseCategorical) %in% currentName)])
    }
    pairwiseCategorical
  }
  
  
  
  
  
  write.csv(mainPairwiseCategorical, file= paste(mainPairwiseCorrelationFileName, ".csv", sep=""), row.names = T, quote = F) #save the correlations as a csv
  saveRDS(mainPairwiseCategorical, paste(mainPairwiseCorrelationFileName, ".rds", sep="")) #and as an rds 
}