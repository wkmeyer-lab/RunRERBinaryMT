.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/convertLogiToNumeric.R")


#---- USAGE -----
#used to combine permulations files made by RunPermulationsManual.R. 

# ARUGMENTS: 
#If an argument contains a '(' it is evaluated as code.
# 'r="filePrefix"'              This is the prefix attached to all files; a required argument. 
# 'e=F' OR 'e=T'                This is if the permulations being combined are enriched or not. Accepts 'T', 'F', 'TRUE', 'FALSE', '0', and '1'. 
# 'i=<number>'                  This is used to generate unique filenames for each instance of the script. Typically fed in by for loop used to run script in parallel.
# 'c=F' OR 'c=T'                This is used to set if the script is being run to combine previous combinations. Called "metacombination". Used for parrallelization. 
# 'm=mainTreeFilename.txt or .rds' This is the location of the maintree file. Accepts .txt or .rds. 
# 't = <s OR f OR p>            This sets which permulation filetype to look for. s is for slow, f is for fast, and p is for pruned-fast



#-------
#Debug setup defaults
#permulationNumberValue = 3
#Testing args:
#args = c('r=allInsectivory','n=5', 'e=F', 's=1', 't=p')
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

enrichValue = F
runInstanceValue = NULL
metacombineValue = FALSE


#-------





# -- Import if enriched or not --
if(!is.na(cmdArgImport('e'))){
  enrichValue = cmdArgImport('e')
  enrichValue = as.logical(enrichValue)
  if(is.na(enrichValue)){
    enrichValue = FALSE
    paste("Enrichment value not interpretable as logical. Did you remember to capitalize? Using FALSE.")
  }
}else{
  paste("enrichment not specified, using FLASE")
}



# -- Import the instance number of the script --- 
if(!is.na(cmdArgImport('i'))){
  runInstanceValue = cmdArgImport('i')
}else{
  paste("This script does not have a run instance value")
}

# -- Import if this is being run on combined correlations -- 
if(!is.na(cmdArgImport('c'))){
  metacombineValue = cmdArgImport('c')
  metacombineValue = as.logical(metacombineValue)
  if(is.na(metacombineValue)){
    metacombineValue = FALSE
    paste("Metacombination value not interpretable as logical. Did you remember to capitalize? Using FALSE.")
  }
}else{
  paste("Metacombination value not specified, using FLASE. If you aren't parrallelizing, don't worry about this.")
}

# -- Import which filetype to use  --- 
if(!is.na(cmdArgImport('t'))){
  fileType = cmdArgImport('t')
  fileType = as.character(fileType)
}else{
  paste("No fileType specified, defaulting to slow (PermulationsData)")
}


# ----- Calculation of p-values

# -- Get clades correlation --

cladesCorellationFileName = paste(outputFolderName, filePrefix, "CladesCorrelationFile", sep= "")
if(file.exists(paste(cladesCorellationFileName, ".rds", sep=""))){
  cladesCorrelation = readRDS(paste(cladesCorellationFileName, ".rds", sep=""))
}else{
  message("Clades Correlation file does not exist, p-values not calculated. (RunPermulationsManual, RunPermulationsFastScript")
  stop("Requires Clades Correlation file.")
}

# ---- Get Combined Permulations filename ---
#Do the first combination:
if(fileType == 's'){
  fileTypeString = ''
}else if(fileType == 'f'){
  fileTypeString = "UnprunedFast"
}else if(fileType == 'p'){
  fileTypeString = "PrunedFast"
}else{
  stop( "THIS IS AN ISSUE MESSAGE, IMPROPER FILETYPE ARGUMENT")
}

if(metacombineValue == F){
  combinedDataFileName = paste(outputFolderName, filePrefix, "Combined", fileTypeString,"PermulationsData", runInstanceValue, ".rds", sep="")
}else{
  combinedDataFileName = paste(outputFolderName, filePrefix, "MetaCombined", fileTypeString, "PermulationsData", runInstanceValue, ".rds", sep="")
}
combinedPermulationsData = readRDS(combinedDataFileName)

  


# -- run pValue calculation --
  permulationPValues = permpvalcor(cladesCorrelation, combinedPermulationsData)
  
  #save the permulations p values
  permulationPValueFileName = paste(outputFolderName, filePrefix, "Combined", fileTypeString, "PermulationsPValue", runInstanceValue, ".rds", sep= "")
  saveRDS(permulationPValues, file = permulationPValueFileName)


