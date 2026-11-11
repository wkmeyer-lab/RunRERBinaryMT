# -- Libraries 
source("Src/Reu/loadLibraries.R")
loadLibraries() # Note that this library load includes cluster run detection. 

source("Src/Reu/CategoricalPermulationsParallelFunctions.R")

# -- Usage:
# This script converts permulated trees into per-gene per-permulation results, called "Intermediates". 

#-------------------------------
# --- Command arguments list ---
# - Required - 
# r = filePrefix                            REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# m = mainTreeFilename.txt or .rds          This sets the location of the maintrees file

# - Required, Defaults likely acceptable -
# f = forceUpdate <T or F>                  This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# j = runInstanceValue                      This is used to generate unique filenames for each instance of the script. Used in parrallelization. Set as j to avoid overlap with alternates.
# - Optional - 
# h = <T>                                   This ARG causes the script to return its ARG list, instead of running the script. 
# b = <T or F>                              This is used to set if the script is targeting metacombined permulations. Called "metacombination". Used for parrallelization. 
# - Overrides - 
# y = targetFileInstance              This is used to manually set an instance value of the permulations data to be loaded (ie, instance "Dev"). Effects the filename to be loaded in. Defaults to this script's run instance value. Used in parallelization and debugging. 

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
# Defaults
mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"  #this is the standard location on Sol 
metacombineValue = FALSE
runInstanceValue = NULL
targetFileInstance = NULL

{ # Bracket used for collapsing purposes
  #MainTrees Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using default")
  }
  #Metacombination 
  if(!is.na(cmdArgImport('b'))){
    metacombineValue = cmdArgImport('b')
    metacombineValue = as.logical(metacombineValue)
    if(is.na(metacombineValue)){
      metacombineValue = FALSE
      message("Metacombination value not interpretable as logical. Did you remember to capitalize? Using FALSE.")
    }
  }else{
    message("Metacombination value not specified, using FALSE. If you aren't parrallelizing, don't worry about this.")
  }
  #Run instance value
  if(!is.na(cmdArgImport('j'))){
    runInstanceValue = cmdArgImport('j')
  }else{
    message("This script does not have a run instance value")
  }
  #target File Instance
  if(!is.na(cmdArgImport('y'))){
    targetFileInstance = cmdArgImport('y')
  }else{
    targetFileInstance = runInstanceValue
    message("No target file Instance specified, defaulting to this scripts run instance.")
  }
  
}


#                   ------- Code Body -------- 

# -- Import data -- 
#MainTrees
if(file_ext(mainTreesLocation) == "rds"){                                       #if the tree is an RDS file
  mainTrees = readRDS(mainTreesLocation)                                        #Read as RDS
}else{                                                                          #Otherwise
  mainTrees = readTrees(mainTreesLocation)                                      #read as text
}

#Phenotype Vector 
phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
phenotypeVector = readRDS(phenotypeVectorFilename)                              #Load the phenotype vector 

#RERs 
RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
RERObject = readRDS(RERFileName)                                                #Load the RERs 

#Paths
pathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "") #make a filename based on the prefix
pathsObject = readRDS(pathsFilename)                                            #Load the paths

#Correlations
correlationFileName = paste(outputFolderName, filePrefix, "CombinedCategoricalCorrelationFile.rds", sep= "") #Make a correlation filename based on the prefix
correlationsObject = readRDS(correlationFileName) 

#Permulations Trees
permulationsDataFileName = paste(outputFolderName, filePrefix, "PermulationsData", targetFileInstance, ".rds", sep= "")
permulationsData = readRDS(permulationsDataFileName)


# -- calculate p values -- 

permCorrelations = CategoricalPermulationGetCor(correlationsObject, permulationsData$trees, phenotypeVector, mainTrees, RERObject, report=T)
message("Correlations complete.")
permulationIntermediateFilename =  paste(outputFolderName, filePrefix, "CategoricalPermulationsIntermediates", runInstanceValue, ".rds", sep= "")
message("Attempting to save file.")
saveRDS(permCorrelations, permulationIntermediateFilename)
