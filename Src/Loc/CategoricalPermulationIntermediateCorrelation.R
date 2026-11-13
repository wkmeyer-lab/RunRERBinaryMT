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
# b = Metacombination <T or F>                This is used to set if the script is being run to combine previous combinations. Called "metacombination". Used for parrallelization. 
# - Overrides - 
# y = targetFileInstance                    This is used to manually set an instance value of the permulations data to be loaded (ie, instance "Dev"). Effects the filename to be loaded in. Defaults to this script's run instance value. Used in parallelization and debugging. 

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
  #Defaults
  mainTreesLocation = "Results/mainTreesFromAlignments.rds"
  runInstanceValue = NULL
  metacombineValue = FALSE
  targetFileInstance = NULL
  
  #Load from args 
  source("Src/Reu/ArgumentImport.R")
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
