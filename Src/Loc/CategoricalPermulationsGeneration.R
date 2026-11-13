# -- Usage:
# This script generates permulated species sets and trees for a given phenotype. 
# THIS IS THE ONLY PERMULATIONS SET WITH A NON-DETERMINISTIC OUTCOME. 
# If you wish to save your exact permulations inputs at a much smaller file size, save only the outputs of this script ('PermulationsData' files). 
# The rest of the results can be deterministicly re-generated, though at a non-trivial amount of compute time. 

# -- Libraries 
source("Src/Reu/loadLibraries.R")
loadLibraries() # Note that this library load includes cluster run detection. 
source("Src/Reu/RelaxedRejectionPermFuncs2.R")
#-------------------------------
# --- Command arguments list ---
# - Required - 
# r = filePrefix                            REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# m = mainTreeFilename.txt or .rds          This sets the location of the maintrees file
# n = numberOfPermulations                  This is the number of permulations to run in the script 

# - Required, Defaults likely acceptable -
# f = forceUpdate <T or F>                  This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# q = <ER or SYM or ARD>                    This sets the model type used to estimate ancestral branches; MUST BE THE SAME AS ONE USED IN MakeCategoricalPhenotypeTree.R
# j = runInstanceValuePerms                 This is used to generate unique filenames for each instance of the script. Used in parrallelization. Set as j to avoid overlap with alternates.
# u = <auto or stationary or flat>          This sets the probability of each state at the root of the tree.   
# - Optional - 
# h = <T>                                   This ARG causes the script to return its ARG list, instead of running the script. 
# x = relaxationValue <int, 0-1>            This is a value between 0 and 1, the percentage off of exact a permulation is allowed to be. Increases runspeed, but reduces permulation accuracy. Recommend 0.1 if more than 3 categories
# - Overrides - 
# --- End arguments list ---
#-------------------------------
# -- Section to include manually-specified args --
if(!exists("args", envir = environment(), inherits = FALSE)){
  #This is a space where you can manually write in arguments for the script. 
  #The above if statement ensures these lines will not override any existing args if the script is run whole-sale. 
  }
#-------------------------------



#----------------

# --- Standard start-up code; identical in all scripts---
source("Src/Reu/StandardizedStartupCode.R") #This runs the standardized code block at this location, ensuring identity across scripts. 
#This code reads the "report args" help ARG, and reutrns the args if it is active. 
#It then reads the filePrefix and sets the output directory to that file prefix, and reads the forceUpdate arg. 


# --- Argument Imports ---
{ # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
  #Defaults
  mainTreesLocation = "Results/mainTreesFromAlignments.rds"
  permulationAmount = 100
  modelType = "ER"
  runInstanceValue = NULL
  rootProbability = "auto"
  useRelaxation = TRUE
  relaxationValue = 0.1
  
  #Load from args 
  source("Src/Reu/ArgumentImport.R")
  
  if(relaxationValue ==0){useRelaxation = FALSE}
}

# ---------------------------
# ----- SCRIPT MAIN CODE ----
# ---------------------------
#                   ------- Code Body -------- 

# -- Read in files -- 
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
#RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
#RERObject = readRDS(RERFileName)                                                #Load the RERs 

#Paths
#pathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "") #make a filename based on the prefix
#pathsObject = readRDS(pathsFilename)                                            #Load the paths

#Correlations
#correlationFileName = paste(outputFolderName, filePrefix, "CorrelationFile.rds", sep= "") #Make a correlation filename based on the prefix
#correlationsObject = readRDS(correlationFileName)                               #Load the correlations

# -- Run Permulations --

permsStartTime = Sys.time()                                                     #get the time before start of permulations
if(useRelaxation){
  permulationData = categoricalPermulations(mainTrees, phenotypeVector, rm = modelType, rp = rootProbability, ntrees = permulationAmount, percent_relax = relaxationValue)
}else{
  permulationData = categoricalPermulations(mainTrees, phenotypeVector, rm = modelType, rp = rootProbability, ntrees = permulationAmount)
}
permsEndTime = Sys.time()                                                       #get time at end of permulations
message(paste("Time to run permulations: ", (permsEndTime - permsStartTime), attr(permsEndTime - permsStartTime, "units"))) #Print the time taken to calculate permulations
message(paste("Time per permulation: ", ((permsEndTime - permsStartTime)/permulationAmount), attr((permsEndTime - permsStartTime)/permulationAmount, "units")))

# -- Save Permulations -- 

saveStartTime = Sys.time()
permulationsDataFileName = paste(outputFolderName, filePrefix, "PermulationsData", runInstanceValue, ".rds", sep= "")
saveRDS(permulationData, permulationsDataFileName)
saveEndTime= Sys.time()
message(paste("Time to save permulations: ", (saveEndTime - saveStartTime), attr(saveEndTime - saveStartTime, "units")))
