# -- Libraries 
source("Src/Reu/loadLibraries.R")
loadLibraries() # Note that this library load includes cluster run detection. 

source("Src/Reu/CategoricalPermulationsParallelFunctions.R")


# -- Usage:
# Script is used to combine permulations into a single file. 

# -- Command arguments list


#-------------------------------
# --- Command arguments list ---
# - Required - 
# r = filePrefix                              REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 

# - Required, Defaults likely acceptable -
# f = forceUpdate <T or F>                    This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# j = runInstanceValue                        This is used to generate unique filenames for each instance of the script. Used in parrallelization. Set as j to avoid overlap with alternates.
# nc = <number Of Permulations Files>         This is the number of permulation files the script will try to combine
# st = <Start number>                         This is the permulation number to start at. Used for parrallelization. 

# - Optional - 
# h = <T>                                     This ARG causes the script to return its ARG list, instead of running the script. 
# b = Metacombination <T or F>                This is used to set if the script is being run to combine previous combinations. Called "metacombination". Used for parrallelization. 
# d = CalculatePValue <T or F>                This value determines whether or not to calculate the pValues.
# v = onlyCalulateValue <T or F>              This value is used to run only the p value calculation, and not the combination. 
# - Overrides - 
# y = "PermulationPrefix"                     This is used to add a prefix (eg. FastPruned) to permulations if present in the filename.

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
{ # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
  runInstanceValue = NULL
  permulationNumberValue = NULL
  startValue = 1
  metacombineValue = FALSE
  calulateValue = TRUE
  onlyCalulateValue = FALSE
  permulationPrefix = NULL
  
  #Load from args 
  source("Src/Reu/ArgumentImport.R")
}



#                   ------- Code Body -------- 

# -- determine the filename -- 
if(!onlyCalulateValue){
  if(metacombineValue == F){
    basePermulationsFilename = paste(outputFolderName, filePrefix, permulationPrefix, "CategoricalPermulationsIntermediates", sep= "")
  
    correlationFileName = paste(outputFolderName, filePrefix, "CombinedCategoricalCorrelationFile.rds", sep= "") #Make a correlation filename based on the prefix
    correlationsObject = readRDS(correlationFileName) 
    
    
  
      # -- Do initial combination (before loop) --
      combinationSectionStart = Sys.time()
      # - First permulation file
      firstPermLoadStart = Sys.time()
      firstPermulationsFilename = paste(basePermulationsFilename, startValue, ".rds", sep="")
      firstPermulationsData = readRDS(firstPermulationsFilename)
      
      firstPermLoadEnd = Sys.time()
      firstPermLoadTime = firstPermLoadEnd - firstPermLoadStart
      message("First permulation load time: ", firstPermLoadTime, attr(firstPermLoadTime, "units"))
      
      CombinedIntermediates = CategoricalCollectIntermediateResults(correlationsObject, firstPermulationsData, initial = T)
      
      rm(firstPermulationsData)
      # -- Do all subsequent combinations (loop) --
      
      iteratingPermulationsData = NULL #make a dummy variable to be removed on first instance of loop 
      
      if((startValue) < (startValue+permulationNumberValue-1)){                     #Sanity check that there are additional combinations to loop through
        for(i in (startValue+1):(startValue+permulationNumberValue-1)){
          rm(iteratingPermulationsData)
          gc()
          tryCatch({
            message(i)
            iteratingPermulationsFilename = paste(basePermulationsFilename, i, ".rds", sep="")
            
            if(file.exists(iteratingPermulationsFilename)){                             #This check allows it to function even if the start+numberOfPermulations is larger than the actual number of permulation files, or a file is missing.
              iteratingPermulationStart = Sys.time()
              
              iteratingPermulationsData = readRDS(iteratingPermulationsFilename)
              iteratingPermulationLoadEnd= Sys.time()
              iteratingPermulationLoadTime = iteratingPermulationLoadEnd - iteratingPermulationStart
              message("Iterating permulation load end time:", Sys.time())
              message("Iterating permulation load time: ", iteratingPermulationLoadTime, attr(iteratingPermulationLoadTime, "units"))
              
              #Moving this function into the main code body to save on memory use 
              #combinedPermulationsData = combineCategoricalPermulationIntermediates(combinedPermulationsData, iteratingPermulationsData)
              
              CombinedIntermediates = CategoricalCollectIntermediateResults(CombinedIntermediates, iteratingPermulationsData)
              
              
              
              iteratingPermulationCombineEnd = Sys.time()
              iteratingPermulationCombineTime = iteratingPermulationCombineEnd - iteratingPermulationLoadEnd
              message("Iterating permulation combination end time:", Sys.time())
              message("Iterating permulation combination time: ", iteratingPermulationCombineTime, attr(iteratingPermulationCombineTime, "units"))
              
              rm(iteratingPermulationsData)
              iteratingPermulationRemoveEnd = Sys.time()
              iteratingPermulationRemoveTime = iteratingPermulationRemoveEnd - iteratingPermulationCombineEnd
              message("Iterating permulation removal end time:", Sys.time())
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
        combinedDataFileName = paste(outputFolderName, filePrefix, "Collected", permulationPrefix,"PermulationsIntermediatesExtemity", runInstanceValue, ".rds", sep="")
      }
      fileSavingStart = Sys.time()
      saveRDS(CombinedIntermediates, file = combinedDataFileName)
      fileSavingEnd = Sys.time()
      fileSavingTime = fileSavingEnd - fileSavingStart
      message("Time to save combine permulations: ", fileSavingTime, attr(fileSavingTime, "units"))
      
  }else{
    combinationSectionStart = Sys.time()
    basePermulationsFilename = paste(outputFolderName, filePrefix, "Collected", permulationPrefix, "PermulationsIntermediatesExtemity",  sep="")
    
    # - First permulation file
    firstPermLoadStart = Sys.time()
    firstPermulationsFilename = paste(basePermulationsFilename, startValue, ".rds", sep="")
    firstPermulationsData = readRDS(firstPermulationsFilename)
    
    firstPermLoadEnd = Sys.time()
    firstPermLoadTime = firstPermLoadEnd - firstPermLoadStart
    message("First permulation load time: ", firstPermLoadTime, attr(firstPermLoadTime, "units"))
    
    CombinedIntermediates = firstPermulationsData
    
    rm(firstPermulationsData)
    # -- Do all subsequent combinations (loop) --
    
    iteratingPermulationsData = NULL #make a dummy variable to be removed on first instance of loop 
    
    if((startValue) < (startValue+permulationNumberValue-1)){                     #Sanity check that there are additional combinations to loop through
      for(i in (startValue+1):(startValue+permulationNumberValue-1)){
        rm(iteratingPermulationsData)
        gc()
        tryCatch({
          message(i)
          iteratingPermulationsFilename = paste(basePermulationsFilename, i, ".rds", sep="")
          
          if(file.exists(iteratingPermulationsFilename)){                             #This check allows it to function even if the start+numberOfPermulations is larger than the actual number of permulation files, or a file is missing.
            iteratingPermulationStart = Sys.time()
            
            iteratingPermulationsData = readRDS(iteratingPermulationsFilename)
            iteratingPermulationLoadEnd= Sys.time()
            iteratingPermulationLoadTime = iteratingPermulationLoadEnd - iteratingPermulationStart
            message("Iterating permulation load end time:", Sys.time())
            message("Iterating permulation load time: ", iteratingPermulationLoadTime, attr(iteratingPermulationLoadTime, "units"))
            #Moving this function into the main code body to save on memory use 
            #combinedPermulationsData = combineCategoricalPermulationIntermediates(combinedPermulationsData, iteratingPermulationsData)
            
            CombinedIntermediates[[1]]$numMoreExtremePerms = CombinedIntermediates[[1]]$numMoreExtremePerms + iteratingPermulationsData[[1]]$numMoreExtremePerms  
            CombinedIntermediates[[1]]$numTotalPerms = CombinedIntermediates[[1]]$numTotalPerms + iteratingPermulationsData[[1]]$numTotalPerms  
            
            
            for (j in 1:length(CombinedIntermediates[[2]])) {
              CombinedIntermediates[[2]][[j]]$numMoreExtremePerms = CombinedIntermediates[[2]][[j]]$numMoreExtremePerms + iteratingPermulationsData[[2]][[j]]$numMoreExtremePerms
              CombinedIntermediates[[2]][[j]]$numTotalPerms = CombinedIntermediates[[2]][[j]]$numTotalPerms + iteratingPermulationsData[[2]][[j]]$numTotalPerms
            }
            
            
            iteratingPermulationCombineEnd = Sys.time()
            iteratingPermulationCombineTime = iteratingPermulationCombineEnd - iteratingPermulationLoadEnd
            message("Iterating permulation combination end time:", Sys.time())
            message("Iterating permulation combination time: ", iteratingPermulationCombineTime, attr(iteratingPermulationCombineTime, "units"))
            
            rm(iteratingPermulationsData)
            iteratingPermulationRemoveEnd = Sys.time()
            iteratingPermulationRemoveTime = iteratingPermulationRemoveEnd - iteratingPermulationCombineEnd
            message("Iterating permulation removal end time:", Sys.time())
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
  
    combinedDataFileName = paste(outputFolderName, filePrefix, "MetaCollected", permulationPrefix, "PermulationsIntermediatesExtemity", runInstanceValue, ".rds", sep="")
  
    fileSavingStart = Sys.time()
    saveRDS(CombinedIntermediates, file = combinedDataFileName)
    fileSavingEnd = Sys.time()
    fileSavingTime = fileSavingEnd - fileSavingStart
    message("Time to save combine permulations: ", fileSavingTime, attr(fileSavingTime, "units"))
  }
}  
  




# -- Calculate p values -- 

if(calulateValue){


  if(onlyCalulateValue){
    if(metacombineValue == F){
      combinedDataFileName = paste(outputFolderName, filePrefix, "Collected", permulationPrefix,"PermulationsIntermediatesExtemity", runInstanceValue, ".rds", sep="")
    }else{
      combinedDataFileName = paste(outputFolderName, filePrefix, "MetaCollected", permulationPrefix, "PermulationsIntermediatesExtemity", runInstanceValue, ".rds", sep="")
    }
    CombinedIntermediates = readRDS(combinedDataFileName)
  }
  permulationsPValuesOutput = CategoricalCalculatePValueFromCollectedIntermediates(CombinedIntermediates)
  
  
  
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






