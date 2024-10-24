# -- Libraries 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge)
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/combineCategoricalPermulationsIntermediates.R")
library(data.table)

# -- Usage:
# Script is used to combine permulations into a single file. 

# -- Command arguments list
# r = filePrefix                      This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                        This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# i = runInstanceValue                This is used to generate unique filenames for each instance of the script. Used in parrallelization. 
# n = <number Of Permulations>        This is the number of permulation files the script will try to combine
# s = <Start number>                  This is the permulation number to start at. Used for parrallelization. 
# t = <PermulationPrefix>             This is used to add a prefix (eg. FastPruned) to permulations if present in the filename.
# c = <T or F>                        This is used to set if the script is being run to combine previous combinations. Called "metacombination". Used for parrallelization. 
# p = <T or F>                        This value determines whether or not to calculate the pValues.
# o = <T or F>                        This value is used to run only the p value calculation, and not the combination. 

#----------------
args = (c('r=CategoricalDiet3Phen', "n=2", "t=Dev"))
# --- Standard start-up code ---
args = commandArgs(trailingOnly = TRUE)
{  # Bracket used for collapsing purposes
  #File Prefix
  if(!is.na(cmdArgImport('r'))){                                                #This cmdArgImport script is a way to import arguments from the command line. 
    filePrefix = cmdArgImport('r')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
  }
  
  #  Output Directory 
  if(!dir.exists("Output")){                                      #Make output directory if it does not exist
    dir.create("Output")
  }
  outputFolderNameNoSlash = paste("Output/",filePrefix, sep = "") #Set the prefix sub directory
  if(!dir.exists(outputFolderNameNoSlash)){                       #create that directory if it does not exist
    dir.create(outputFolderNameNoSlash)
  }
  outputFolderName = paste("Output/",filePrefix,"/", sep = "")
  
  #  Force update argument
  forceUpdate = FALSE
  if(!is.na(cmdArgImport('v'))){                                 #Import if update being forced with argument 
    forceUpdate = cmdArgImport('v')
    forceUpdate = as.logical(forceUpdate)
  }else{
    paste("Force update not specified, not forcing update")
  }
}

# --- Argument Imports ---
# Defaults
runInstanceValue = NULL
permulationNumberValue = 100
startValue = 1
permulationPrefix = NULL
metacombineValue = FALSE
calulateValue = TRUE
onlyCalulateValue = FALSE

{ # Bracket used for collapsing purposes
  
  #Run Instance Number
  if(!is.na(cmdArgImport('i'))){
    runInstanceValue = cmdArgImport('i')
  }else{
    message("This script does not have a run instance value")
  }
  
  #Number of permulations to combine
  if(!is.na(cmdArgImport('n'))){
    permulationNumberValue = cmdArgImport('n')
    permulationNumberValue = as.numeric(permulationNumberValue)
  }else{
    message("Number of permulations not specified, using 100")
  }
  
  #Start Number
  if(!is.na(cmdArgImport('s'))){
    startValue = cmdArgImport('s')
    startValue = as.numeric(startValue)
  }else{
    message("Start value not specified, using 1")
  }
  
  #Permulation Prefix
  if(!is.na(cmdArgImport('t'))){
    permulationPrefix = cmdArgImport('t')
  }else{
    message("No filename prefix used, using 'PermulationsData'.")
  }
  
  #Metacombination 
  if(!is.na(cmdArgImport('c'))){
    metacombineValue = cmdArgImport('c')
    metacombineValue = as.logical(metacombineValue)
    if(is.na(metacombineValue)){
      metacombineValue = FALSE
      message("Metacombination value not interpretable as logical. Did you remember to capitalize? Using FALSE.")
    }
  }else{
    message("Metacombination value not specified, using FALSE. If you aren't parrallelizing, don't worry about this.")
  }
  
  #Calculate P Value
  if(!is.na(cmdArgImport('p'))){
    calulateValue = cmdArgImport('p')
    calulateValue = as.logical(calulateValue)
    if(is.na(calulateValue)){
      calulateValue = FALSE
      message("p-value calculation value not interpretable as logical. Did you remember to capitalize? Using FALSE.")
    }
  }else{
    message("p-value calculation value not specified, using TRUE. If you aren't parrallelizing, don't worry about this.")
  }
  
  #Only Calculate P Value
  if(!is.na(cmdArgImport('o'))){
    onlyCalulateValue = cmdArgImport('o')
    onlyCalulateValue = as.logical(calulateValue)
    if(is.na(calulateValue)){
      onlyCalulateValue = FALSE
      message("Only p-value calculation value not interpretable as logical. Did you remember to capitalize? Using FALSE.")
    }
  }else{
    message("Only p-value calculation value not specified, using FALSE. If you aren't parrallelizing, don't worry about this.")
  }
  
}



#                   ------- Code Body -------- 

# -- determine the filename -- 
if(metacombineValue == F){
  basePermulationsFilename = paste(outputFolderName, filePrefix, permulationPrefix, "CategoricalPermulationsIntermediates", sep= "")
}else{
  basePermulationsFilename = paste(outputFolderName, filePrefix, "Combined", permulationPrefix, "CategoricalPermulationsIntermediates",  sep="")
}

if(!onlyCalulateValue){
  # -- Do initial combination (before loop) --
  combinationSectionStart = Sys.time()
  # - First permulation file
  firstPermLoadStart = Sys.time()
  firstPermulationsFilename = paste(basePermulationsFilename, startValue, ".rds", sep="")
  firstPermulationsData = readRDS(firstPermulationsFilename)
  
  firstPermLoadEnd = Sys.time()
  firstPermLoadTime = firstPermLoadEnd - firstPermLoadStart
  message("First permulation load time: ", firstPermLoadTime, attr(firstPermLoadTime, "units"))
  
  
  # - Second permulation file - 
  secondPermLoadStart = Sys.time()
  secondPermulationsFilename = paste(basePermulationsFilename, (startValue+1), ".rds", sep="")
  secondPermulationsData = readRDS(secondPermulationsFilename)
  
  secondPermLoadEnd = Sys.time()
  secondPermLoadTime = secondPermLoadEnd - secondPermLoadStart
  message("Second permulation load time: ", secondPermLoadTime, attr(secondPermLoadTime, "units"))
  
  # - Initial combination (makes the CombinedPermulationsData object) - 
  
  firstCombinationStart = Sys.time()
  combinedPermulationsData = combineCategoricalPermulationIntermediates(firstPermulationsData, secondPermulationsData)
  firstCombinationEnd = Sys.time()
  firstCombinationTime = firstCombinationEnd - firstCombinationStart
  message("Initial permulation combination time: ", firstCombinationTime, attr(firstCombinationTime, "units"))
  
  # -- Do all subsequent combinations (loop) --
  
  if((startValue+2) < (startValue+permulationNumberValue-1)){                     #Sanity check that there are additional combinations to loop through
    for(i in (startValue+2):(startValue+permulationNumberValue-1)){
      tryCatch({
      message(i)
      iteratingPermulationsFilename = paste(basePermulationsFilename, i, ".rds", sep="")
      
      if(file.exists(iteratingPermulationsFilename)){                             #This check allows it to function even if the start+numberOfPermulations is larger than the actual number of permulation files, or a file is missing.
        iteratingPermulationStart = Sys.time()
        
        iteratingPermulationsData = readRDS(iteratingPermulationsFilename)
        iteratingPermulationLoadEnd= Sys.time()
        iteratingPermulationLoadTime = iteratingPermulationLoadEnd - iteratingPermulationStart
        message("Iterating permulation load time: ", iteratingPermulationLoadTime, attr(iteratingPermulationLoadTime, "units"))
        
        combinedPermulationsData = combineCategoricalPermulationIntermediates(combinedPermulationsData, iteratingPermulationsData)
        iteratingPermulationCombineEnd = Sys.time()
        iteratingPermulationCombineTime = iteratingPermulationCombineEnd - iteratingPermulationLoadEnd
        message("Iterating permulation combination time: ", iteratingPermulationCombineTime, attr(iteratingPermulationCombineTime, "units"))
        
        rm(iteratingPermulationsData)
        iteratingPermulationRemoveEnd = Sys.time()
        iteratingPermulationRemoveTime = iteratingPermulationRemoveEnd - iteratingPermulationCombineEnd
        message("Iterating permulation removal time: ", iteratingPermulationRemoveTime, attr(iteratingPermulationRemoveTime, "units"))
        
        message("Added file ", i, " to combination.")
      }else{
        message("Permulation file number ", i, " does not exist. Combining other files.")
      }
      }, error = function(i){
        message("This error occurred on this file:") 
        message(i)
        message("Skipping.")
      })
    }
  }
  combinationSectionEnd = Sys.time()
  totalCombinationTime = combinationSectionEnd- combinationSectionStart
  message(" Total Permulation Combination time: ", totalCombinationTime, attr(totalCombinationTime, "units"))
  
  # -- Save combined permulations file -- 
  if(metacombineValue == F){
    combinedDataFileName = paste(outputFolderName, filePrefix, "Combined", permulationPrefix,"PermulationsIntermediates", runInstanceValue, ".rds", sep="")
  }else{
    combinedDataFileName = paste(outputFolderName, filePrefix, "MetaCombined", permulationPrefix, "PermulationsIntermediates", runInstanceValue, ".rds", sep="")
  }
  fileSavingStart = Sys.time()
  saveRDS(combinedPermulationsData, file = combinedDataFileName)
  fileSavingEnd = Sys.time()
  fileSavingTime = fileSavingEnd - fileSavingStart
  message("Time to save combine permulations: ", fileSavingTime, attr(fileSavingTime, "units"))
}

# -- Calculate p values -- 

if(calulateValue){
  source("Src/Reu/CategoricalPermulationsParallelFunctions.R")
  
  #Correlations
  correlationFileName = paste(outputFolderName, filePrefix, "CombinedCategoricalCorrelationFile.rds", sep= "") #Make a correlation filename based on the prefix
  correlationsObject = readRDS(correlationFileName) 
  
  if(onlyCalulateValue){
    if(metacombineValue == F){
      combinedDataFileName = paste(outputFolderName, filePrefix, "Combined", permulationPrefix,"PermulationsIntermediates", runInstanceValue, ".rds", sep="")
    }else{
      combinedDataFileName = paste(outputFolderName, filePrefix, "MetaCombined", permulationPrefix, "PermulationsIntermediates", runInstanceValue, ".rds", sep="")
    }
    combinedPermulationsData = readRDS(combinedDataFileName)
  }
  permulationsPValues = CategoricalCalculatePermulationPValues(correlationsObject, combinedPermulationsData)
  permulationsPValuesOutput = permulationsPValues$res
  
  #Give pairwise outputs descriptive names
  phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #select the phenotype vector based on prefix
  phenotypeVector = readRDS(phenotypeVectorFilename)                            #load in the phenotype vector 
  categories = map_to_state_space(phenotypeVector)                              #and use it to connect branch lengths to phenotype name
  categoryNames = categories$name2index                                         #store the length-phenotype connection
  
  pairwiseTableNames = names(permulationsPValuesOutput[[2]])                               #Prepare to replace the number-number titles with phenotype-phenotype titles
  for(i in 1:length(categoryNames)){                                            #for each phenotype
    pairwiseTableNames= gsub(i, names(categoryNames)[i], pairwiseTableNames)    #replace the number with the phenotype name  
  }
  names(permulationsPValuesOutput[[2]]) = pairwiseTableNames                               #update the dataframe titles
  
  permulationsPValuesFilename = paste(outputFolderName, filePrefix, "PermulationsPValueCorrelations.rds", sep= "")
  saveRDS(permulationsPValuesOutput, permulationsPValuesFilename)
  
  outputSubdirectoryNoslash = paste(outputFolderName, "Overall", sep = "")
  if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
    dir.create(outputSubdirectoryNoslash)
  }
  outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
  
  permulationsPValuesOverallFilename = paste(outputSubdirectory, filePrefix, "OverallPermulationsCorrelationFile.rds", sep= "")
  saveRDS(permulationsPValuesOutput[[1]], permulationsPValuesOverallFilename)
  
  for(i in 1:length(pairwiseTableNames)){
    pairwiseTableNames= gsub(" ", "", pairwiseTableNames)
    
    outputSubdirectoryNoslash = paste(outputFolderName, pairwiseTableNames[i], sep = "")
    if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
      dir.create(outputSubdirectoryNoslash)
    }
    outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
    
    permulationsPValuesPairFilename = paste(outputSubdirectory, filePrefix, pairwiseTableNames[i], "PermulationsCorrelationFile",".rds", sep= "")
    saveRDS(permulationsPValuesOutput[[2]][[i]], permulationsPValuesPairFilename)
  }
}






