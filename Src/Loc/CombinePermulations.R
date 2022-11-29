#Library setup 
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
# 'n=numberOfPermulations'      This is the number of permulation files the script will try to combine
# 'e=F' OR 'e=T'                This is if the permulations being combined are enriched or not. Accepts 'T', 'F', 'TRUE', 'FALSE', '0', and '1'. 
# 's=<number>'                  This is the permulation number to start at. Used for parrallelization. 
# 'i=<number>'                  This is used to generate unique filenames for each instance of the script. Typically fed in by for loop used to run script in parallel.
# 'c=F' OR 'c=T'                This is used to set if the script is being run to combine previous combinations. Called "metacombination". Used for parrallelization. 
# 't = <s OR f OR p>            This sets which permulation filetype to look for. s is for slow, f is for fast, and p is for pruned-fast

#-------
#Debug setup defaults
#permulationNumberValue = 3
#Testing args:
#args = c('r=allInsectivory','n=2', 'e=F', 's=1')
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

permulationNumberValue = 100
enrichValue = F
startValue = 1
runInstanceValue = NULL
metacombineValue = FALSE
fileType = 's'


#-------



# -- Import number of permulations to combine ---
if(!is.na(cmdArgImport('n'))){
  permulationNumberValue = cmdArgImport('n')
  permulationNumberValue = as.numeric(permulationNumberValue)
}else{
  paste("Number of permulations not specified, using 100")
}

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

# -- Import the permulation number to start at ---
if(!is.na(cmdArgImport('s'))){
  startValue = cmdArgImport('s')
  startValue = as.numeric(startValue)
}else{
  paste("Start value not specified, using 1")
}

# -- Import the instance number of the script --- 
if(!is.na(cmdArgImport('i'))){
  runInstanceValue = cmdArgImport('i')
}else{
  paste("This script does not have a run instance value")
}

# -- Import if this is being run to combine combinations -- 
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

# -- Combine permulations data files ----

#Do the first combination:
if(metacombineValue == F){
  if(fileType == 's'){
    basePermulationsFilename = paste(outputFolderName, filePrefix, "PermulationsData",  sep="")
  }else if(fileType == 'f'){
    basePermulationsFilename = paste(outputFolderName, filePrefix, "UnprunedFastPermulationsData",  sep="")
  }else if(fileType == 'p'){
    basePermulationsFilename = paste(outputFolderName, filePrefix, "PrunedFastPermulationsData",  sep="")
  }else{
    stop( "THIS IS AN ISSUE MESSAGE, IMPROPER FILETYPE ARGUMENT")
  }
}else{
  basePermulationsFilename = paste(outputFolderName, filePrefix, "CombinedPermulationsData",  sep="")
}

firstPermulationsFilename = paste(basePermulationsFilename, startValue, ".rds", sep="")
firstPermulationsData = readRDS(firstPermulationsFilename)
firstPermulationsData = convertLogiToNumericList(firstPermulationsData)

secondPermulationsFilename = paste(basePermulationsFilename, (startValue+1), ".rds", sep="")
secondPermulationsData = readRDS(secondPermulationsFilename)
secondPermulationsData = convertLogiToNumericList(secondPermulationsData)

combinedPermulationsData = combinePermData(firstPermulationsData, secondPermulationsData, enrich = enrichValue)

rm(firstPermulationsData)
rm(secondPermulationsData)


#Do all subsequent combinations
for(i in (startValue+2):(startValue+permulationNumberValue-1)){
  message(i)
  iteratingPermulationsFilename = paste(basePermulationsFilename, i, ".rds", sep="")
  
  if(file.exists(iteratingPermulationsFilename)){
    iteratingPermulationsData = readRDS(iteratingPermulationsFilename)
    iteratingPermulationsData = convertLogiToNumericList(iteratingPermulationsData)
    combinedPermulationsData = combinePermData(combinedPermulationsData, iteratingPermulationsData, enrich = enrichValue)
    rm(iteratingPermulationsData)
    message("Added file ", i, " to combination.")
  }else{
    message("Permulation file number ", i, " does not exist. Combining other files.")
  }
}

# ---- Calculate the pValues -----

# -- Get clades correlation --
#This relies on a correlationFile produced by the runPermulations initial scripts. If that file doesn't exist, it skips the step and goes directly to saving the combined permulations.

cladesCorellationFileName = paste(outputFolderName, filePrefix, "CladesCorrelationFile", sep= "")
if(file.exists(paste(cladesCorellationFileName, ".rds", sep=""))){
  cladesCorrelation = readRDS(paste(cladesCorellationFileName, ".rds", sep=""))

  # -- run pValue calculation --
  permulationPValues = permpvalcor(cladesCorrelation, combinedPermulationsData)
  
  #save the permulations p values
  permulationPValueFileName = paste(outputFolderName, filePrefix, "CombinedPermulationsPValue", runInstanceValue, ".rds", sep= "")
  saveRDS(permulationPValues, file = permulationPValueFileName)
}else{
  message("Clades Correlation file does not exist, p-values not calculated. (RunPermulationsManual, RunPermulationsFastScript")
}

# ---- Save Combined Permulations output as a file ---
if(metacombineValue == F){
  combinedDataFileName = paste(outputFolderName, filePrefix, "CombinedPermulationsData", runInstanceValue, ".rds", sep="")
}else{
  combinedDataFileName = paste(outputFolderName, filePrefix, "MetaCombinedPermulationsData", runInstanceValue, ".rds", sep="")
  
}
#TEMPORARY REDIRECT
#combinedDataFileName = paste("/home/mit221/wym219_123121/mit221/", filePrefix, "CombinedPermulationsData", runInstanceValue, ".rds", sep="")
saveRDS(combinedPermulationsData, file = combinedDataFileName)












# ---- Debug Code ----
#testCombinedPermsData = combinePermData(firstPermulationsData, iteratingPermulationsData, enrich = enrichValue)
#testCombinedPermsDataTwo = combinePermData(testCombinedPermsData, secondPermulationsData, enrich = enrichValue)
#i=3

# --- Permulations cleaning code; for use elsewhere ----

#filtering version
#cleanedfirstPermulationsData = firstPermulationsData
#cleanedfirstPermulationsData$corP = Filter(function(x)!all(is.na(x)), firstPermulationsData$corP)
#cleanedfirstPermulationsData$corRho = Filter(function(x)!all(is.na(x)), firstPermulationsData$corRho)
#cleanedfirstPermulationsData$corStat = Filter(function(x)!all(is.na(x)), firstPermulationsData$corStat)


#converting version
#convertLogiToNumericDataframe = function(inputDataframe){
#  collumnsLogical = sapply(inputDataframe, is.logical)
#  inputDataframe[,collumnsLogical] = lapply(inputDataframe[,collumnsLogical], as.numeric)
#  return(inputDataframe)
#}

#convertLogiToNumericList = function(inputList){
#  framesInList = sapply(inputList, is.data.frame)
#  inputList[framesInList] = lapply(inputList[framesInList], convertLogiToNumericDataframe)
#  return(inputList)
#}
