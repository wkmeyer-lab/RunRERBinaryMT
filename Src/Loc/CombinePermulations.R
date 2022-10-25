#Library setup 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")

#-------
#Debug setup defaults
#permulationNumberValue = 3
#------

# --- Import prefix ----
args = commandArgs(trailingOnly = TRUE)
paste(args)
message(args)

#File Prefix
if(!is.na(cmdArgImport('r'))){
  filePrefix = cmdArgImport('r')
}else{
  stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
}

#------------

#-------
#Debug setup defaults
#permulationNumberValue = 3
#------

#------ Make Output directory -----

#Make output directory if it does not exist
if(!dir.exists("Output")){
  dir.create("Output")
}
#Make a specific subdirectory if it does not exist 
outputFolderNameNoSlash = paste("Output/",filePrefix, sep = "")
#create that directory if it does not exist
if(!dir.exists(outputFolderNameNoSlash)){
  dir.create(outputFolderNameNoSlash)
}
outputFolderName = paste("Output/",filePrefix,"/", sep = "")

#----------

#----- Default values -------

permulationNumberValue = 100

#-------


# -- Import number of permulations to combine ---
#Number of permulations
if(!is.na(cmdArgImport('n'))){
  permulationNumberValue = cmdArgImport('n')
  permulationNumberValue = as.numeric(permulationNumberValue)
}else{
  paste("Number of permulations not specified, using 100")
}





# -- Combine permulations data files ----

#Do the first combination:
basePermulationFilename = permualationsDataFileName = paste(outputFolderName, filePrefix, "PermulationsData")
firstPermulationFilename = paste(basePermulationFilename, 1, ".rds", sep="")
secondPermulationFilename = paste(basePermulationFilename, 2, ".rds", sep="")

combinedPermulationsData = combinePermData(firstPermulationFilename, secondPermulationFilename)

#Do all subsequent combinations
for(i in 3:permulationNumberValue){
  iteratingPermulationFilename = paste(basePermulationFilename, i, ".rds", sep="")
  if(file.exists(iteratingPermulationFilename)){
   combinedPermulationsData = combinePermData(combinedPermulationsData, interatingPermulationFilename)
   paste("Added file", i, "to combination.")
  }else{
    paste("Permulation file number", i, "does not exist. Combining other files.")
  }
}

#save output as file
combinedDataFileName = paste(outputFolderName, filePrefix, "combinedPermulationsData.rds")
saveRDS(combinedPermulationsData, file = combinedDataFileName)


??permulations
