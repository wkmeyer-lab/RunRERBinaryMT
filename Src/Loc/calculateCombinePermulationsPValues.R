.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
Sys.setenv('R_MAX_VSIZE'=20000000000)
library(RERconverge) #load RERconverge package
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/convertLogiToNumeric.R")
source("Src/Loc/permPValCorReport.R")

#---- USAGE -----
#used to calculate p values on combined permulations files made by CombinePermualtions.R. 

# ARUGMENTS: 
#If an argument contains a '(' it is evaluated as code.
# 'r="filePrefix"'              This is the prefix attached to all files; a required argument. 
# 'i=<number>'                  This is used to generate unique filenames for each instance of the script. Typically fed in by for loop used to run script in parallel.
# 'c=F' OR 'c=T'                This is used to set if the script is being run to combine previous combinations. Called "metacombination". Used for parrallelization. 
# 't = <s OR f OR p>'           This sets which permulation filetype to look for. s is for slow, f is for fast, and p is for pruned-fast
# 's = <integer>'               This sets the gene number to start at for parallelization,                 
# 'n = <integer>                This is the number of genes to do, used to parallelization

#-------
#Debug setup defaults
#permulationNumberValue = 3
#Testing args:
args = c('r=demoinsectivory', 'n=3', 'e=F', 't=p')
args = c('r=allInsectivory', 'e=F', 's=1', 't=s', 'i=1')
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

runInstanceValue = NULL
metacombineValue = FALSE
startValue = 1
geneNumberValue = NA #This means that by defulat it does all of the genes

#-------





# -- Import if enriched or not --



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

# -- Import the permulation number to start at ---
if(!is.na(cmdArgImport('s'))){
  startValue = cmdArgImport('s')
  startValue = as.numeric(startValue)
}else{
  paste("Start value not specified, using 1")
}

# -- Import the number of genes to do 
if(!is.na(cmdArgImport('n'))){
  geneNumberValue = cmdArgImport('n')
  geneNumberValue = as.numeric(geneNumberValue)
}else{
  paste("Number ofgenes not specified, using all")
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

timeBeforeReadIn = Sys.time()
combinedPermulationsData = readRDS(combinedDataFileName)
timeAferReadIn = Sys.time()
readingTime = timeAferReadIn - timeBeforeReadIn
message("File reading time: ", readingTime, attr(readingTime, "units"))
  
#Explict order for garbage collection 
gc() 
# -- run pValue calculation --
  permulationPValues = permPValCorReport(cladesCorrelation, combinedPermulationsData, startNumber = startValue, geneNumber = geneNumberValue)
  
  #save the permulations p values
    #Make generate a marker saying which genes were calculated 
  if(geneNumberValue = NA){
    geneRange = "All"
  }else{
    geneRange = paste(startValue, "-", (startValue + geneNumberValue), sep = "" )
  }
  
  permulationPValueFileName = paste(outputFolderName, filePrefix, "Combined", fileTypeString, geneRange, "PermulationsPValue", runInstanceValue, ".rds", sep= "")
  saveRDS(permulationPValues, file = permulationPValueFileName)


