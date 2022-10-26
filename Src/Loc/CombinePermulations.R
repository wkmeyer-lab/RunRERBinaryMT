#Library setup 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")


#---- USAGE -----
#used to combine permulations files made by RunPermulationsManual.R. 

# ARUGMENTS: 
#If an argument contains a '(' it is evaluated as code.
# 'r="filePrefix"'              This is the prefix attached to all files a required argument. 
# 'n=numberOfPermulations'      This is the number of permulation files the script will try to combine
# 'e=F' OR 'e=T'                This is if the permulations being combined are enriched or not. Accepts 'T', 'F', 'TRUE', 'FALSE', '0', and '1'. 



#-------
#Debug setup defaults
#permulationNumberValue = 3
#Testing args:
#args = c('r=allInsectivory','n=4', 'e=F')
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


# -- Combine permulations data files ----

#Do the first combination:
basePermulationsFilename = paste(outputFolderName, filePrefix, "PermulationsData",  sep="")

firstPermulationsFilename = paste(basePermulationsFilename, 1, ".rds", sep="")
firstPermulationsData = readRDS(firstPermulationsFilename)

secondPermulationsFilename = paste(basePermulationsFilename, 2, ".rds", sep="")
secondPermulationsData = readRDS(secondPermulationsFilename)

combinedPermulationsData = combinePermData(firstPermulationsData, secondPermulationsData, enrich = enrichValue)

rm(firstPermulationsData)
rm(secondPermulationsData)


#Do all subsequent combinations
for(i in 3:permulationNumberValue){
  message(i)
  iteratingPermulationsFilename = paste(basePermulationsFilename, i, ".rds", sep="")
  
  if(file.exists(iteratingPermulationsFilename)){
    iteratingPermulationsData = readRDS(iteratingPermulationsFilename)
    combinedPermulationsData = combinePermData(combinedPermulationsData, iteratingPermulationsData, enrich = enrichValue)
    rm(iteratingPermulationsData)
    message("Added file ", i, " to combination.")
  }else{
    message("Permulation file number ", i, " does not exist. Combining other files.")
  }
}

# save output as file
combinedDataFileName = paste(outputFolderName, filePrefix, "combinedPermulationsData.rds")
saveRDS(combinedPermulationsData, file = combinedDataFileName)

# ---- Debug Code ----
#testCombinedPermsData = combinePermData(firstPermulationsData, iteratingPermulationsData, enrich = enrichValue)
#testCombinedPermsDataTwo = combinePermData(testCombinedPermsData, secondPermulationsData, enrich = enrichValue)
#i=3

# --- Permulations cleaning code; for use elsewhere ----
cleanedfirstPermulationsData = firstPermulationsData
cleanedfirstPermulationsData$corP = Filter(function(x)!all(is.na(x)), firstPermulationsData$corP)
cleanedfirstPermulationsData$corRho = Filter(function(x)!all(is.na(x)), firstPermulationsData$corRho)
cleanedfirstPermulationsData$corStat = Filter(function(x)!all(is.na(x)), firstPermulationsData$corStat)

firstPermulationsData$corP$X142
numericTest = as.numeric(firstPermulationsData$corP$X23)
numerifiedfirstPermulationsData = as.numeric(firstPermulationsData$corP)


