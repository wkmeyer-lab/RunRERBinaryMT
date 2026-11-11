# -- Libraries 
source("Src/Reu/loadLibraries.R")
loadLibraries() # Note that this library load includes cluster run detection. 


# -- Usage:
# This script runs permulations for a continuous phenotype. 

#-------------------------------
# --- Command arguments list ---
# - Required - 
# r = filePrefix                            REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# m = mainTreeFilename.txt or .rds    This sets the location of the maintrees file
# n = numberOfPermulations            This is the number of permulations to run in the script 
# o = outgroupSpecies                 This is the outgroup species with the tip name as it appears on mainTrees$masterTree


# - Required, Defaults likely acceptable -
# f = forceUpdate <T or F>                    This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# j = runInstanceValue                      This is used to generate unique filenames for each instance of the script. Used in parrallelization. Set as j to avoid overlap with alternates.
# - Optional - 
# h = <T>                                   This ARG causes the script to return its ARG list, instead of running the script. 
# p = <T or F>                        This sets if the script should calculate permulation p vals. ONLY USE THIS IF NOT PARRALLELIZING, as it only considers permulations from this script.
# - Overrides - 
# --- End arguments list ---
#-------------------------------
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
permulationAmount = 100
runInstanceValue = NULL
outgroupVal = "Debug"
calculatePVal = F

{ # Bracket used for collapsing purposes
  #MainTrees Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using default")
  }
  
  #Number of permulations
  if(!is.na(cmdArgImport('n'))){
    permulationAmount = cmdArgImport('n')
    permulationAmount = as.numeric(permulationAmount)
  }else{
    paste("Number of permulations not specified, using 100")
  }
  #Run instance value
  if(!is.na(cmdArgImport('j'))){
    runInstanceValue = cmdArgImport('j')
  }else{
    paste("This script does not have a run instance value")
  }
  #Outgroup Species
  if(!is.na(cmdArgImport('o'))){
    outgroupVal = cmdArgImport('o')
  }else{
    stop("An outgroup is required, specify outgroup species")
  }
  #Calculate P 
  if(!is.na(cmdArgImport('p'))){
    calculatePVal = cmdArgImport('p')
  }else{
    paste("p caluclation not specified, not interally calculating p values.")
  }
}


#                   ------- Code Body -------- 

# -- Read in files -- 
#MainTrees
if(file_ext(mainTreesLocation) == "rds"){                                       #if the tree is an RDS file
  mainTrees = readRDS(mainTreesLocation)                                        #Read as RDS
}else{                                                                          #Otherwise
  mainTrees = readTrees(mainTreesLocation)                                      #read as text
}

#Phenotype Vector 
phenotypeVectorFilename = paste(outputFolderName, filePrefix, "ContinuousPhenotypeVector.rds",sep="") #make a filename based on the prefix
phenotypeVector = readRDS(phenotypeVectorFilename)                              #Load the phenotype vector 

#RERs 
RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
RERObject = readRDS(RERFileName)                                                #Load the RERs 

#Paths
#pathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "") #make a filename based on the prefix
#pathsObject = readRDS(pathsFilename)                                            #Load the paths

#Correlations
correlationFileName = paste(outputFolderName, filePrefix, "CorrelationFile.rds", sep= "") #Make a correlation filename based on the prefix
correlationsObject = readRDS(correlationFileName)                               #Load the correlations

# -- make rooted Mastertree --

rootedMaster = root.phylo(mainTrees$masterTree, outgroup = outgroupVal, resolve.root = T)
is.binary(rootedMaster)

# -- Run Permulations --

permsStartTime = Sys.time()                                                     #get the time before start of permulations
permulationData = getPermsContinuous(permulationAmount,phenotypeVector,RERObject, trees = mainTrees, mastertree = rootedMaster, calculateenrich = F)
permsEndTime = Sys.time()                                                       #get time at end of permulations
message(paste("Time to run permulations: ", (permsEndTime - permsStartTime), attr(permsEndTime - permsStartTime, "units"))) #Print the time taken to calculate permulations
message(paste("Time per permulation: ", ((permsEndTime - permsStartTime)/permulationAmount), attr((permsEndTime - permsStartTime)/permulationAmount, "units")))

# -- Save Permulations -- 

saveStartTime = Sys.time()
permulationsDataFileName = paste(outputFolderName, filePrefix, "PermulationsData", runInstanceValue, ".rds", sep= "")
saveRDS(permulationData, permulationsDataFileName)
saveEndTime= Sys.time()
message(paste("Time to save permulations: ", (saveEndTime - saveStartTime), attr(saveEndTime - saveStartTime, "units")))

# -- Calculate p values (optional) --

if(calculatePVal){
  permulatedPVals = permpvalcor(correlationsObject, permulationData)
  
  
  permulationsOutputFilename = paste(outputFolderName, filePrefix, "PermulationPValue", sep= "")
  write.csv(permulatedPVals, file= paste(permulationsOutputFilename, ".csv", sep=""), row.names = T, quote = F) #Save correlations as csv
  saveRDS(permulatedPVals, paste(permulationsOutputFilename, ".rds", sep=""))                          #and as an rds 
  
}
