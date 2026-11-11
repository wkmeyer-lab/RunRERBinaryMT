# -- Usage:
# This script is used to generate gene set enrichment values for supplied genesets from the RERConverge Correlations. 
# This script requires that the .gmt files for the gene sets used to be supplied as a full path, because it needs to directly access the data in those files. 
# All files relating to a geneset (including those from later scripts) will have a standardized filename, with a gene set indicator created from the filename of the .gmt file. 
# Subsequent processing scripts only require this "gene set indicator" (the name of the .gmt file without path or .gmt extension), but are capable of processing a filepath to this indicator. 
# This script is also capable of handling "gene" arguments which tell later scripts to run on the gene results. It does this by ignoring any of the inputs which trigger gene processing in later scripts.

# -- Libraries 
source("Src/Reu/loadLibraries.R")
loadLibraries() # Note that this library load includes cluster run detection. 
#-------------------------------
# --- Command arguments list ---
# - Required - 
# r = filePrefix                              REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# g = c("gene", "sets", "file", "paths")      This sets the gene sets (file names) that the script should be run on, or "gene" for use on genes. 

# - Required, Defaults likely acceptable -
# f = forceUpdate <T or F>                    This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# p = permulations use <T or F or C or A>     This sets if the code should use permulated or unpermulated values. If set to "C" or "A", will use permulations for categorical values, which are stored in the main file. 
# s = "subdirectoryName"                      This is used to specify a subdirectory for the analysis to be run in. Primarily used for the components of categorical results. 
# - Optional - 
# h = <T>                                     This ARG causes the script to return its ARG list, instead of running the script. 
# - Overrides - 
# y = "correlationFileOverride.rds"           This can be used to manually set the correlation file. Used primarily to target pairwise comparisons of categorical phenotypes.
# z = "permulationPvalueFileLocation.rds"     This is a manual override to specify the script use a specific Permulation p-value file.
                                                  #If using permulations in a binary analysis, and using any file other than "CombinedBinaryPermulationsPValue" with no run instance number, it must be specified manually.

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

# ---------------------------
# ----- SCRIPT MAIN CODE ----
# ---------------------------

# --- Argument Imports ---


{ # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
  
  # Defaults
  gmtFileLocation = "Data/enrichmentGmtFile.gmt"
  usePermulations = FALSE
  usePermulationPValOverride = FALSE
  permulationPValOverride = NULL 
  useCorrelationOverride = FALSE
  correlationOverride = NULL
  useCategoricalPermulations = FALSE
  useSubdirectory = FALSE
  subdirectoryValueList = NULL
  validGeneArguments = c("NULL", "gene", "Gene")
  
  
  #Gmt file location
  if(!any(is.na(cmdArgImport('g')))){
    gmtFileLocation = cmdArgImport('g')
    
    #remove any arguments which tell later scripts to run a gene analysis 
    gmtFileLocation = gmtFileLocation[-which(gmtFileLocation %in% validGeneArguments)]
    
  }else{
    paste("No gmt location arugment, using Data/enrichmentGmtFile.gmt")                          #Report using default
    message("No gmt location arugment, using Data/enrichmentGmtFile.gmt")
  }
  
  #Permulation use
  if(!is.na(cmdArgImport('p'))){
    usePermulations = cmdArgImport('p')
    if(usePermulations %in% c("C","A")){
      usePermulations = TRUE
      useCategoricalPermulations = TRUE
    }
    usePermulations = as.logical(usePermulations)
    if(is.na(usePermulations)){
      usePermulations = TRUE
      message("Use Permualtions value not interpretable as logical. Did you remember to capitalize? Using TRUE.")
    }
  }else{
    message("Use Permulations value not specified, using TRUE.")
  }
  if(usePermulations){
    #Import permulation p-value Override 
    if(!is.na(cmdArgImport('z'))){
      usePermulationPValOverride = TRUE
      permulationPValOverride = cmdArgImport('z')
      message("Using Manually specified permulation p-value file.")
    }else{
      if(!useCategoricalPermulations){
        message("No Permulation pValue override, using standard CombinedBinaryPermulationsPValue.rds.")
      }
    }
  }
  #Import correlation file override
  if(!is.na(cmdArgImport('y'))){
    useCorrelationOverride = TRUE
    correlationOverride = cmdArgImport('y')
    message("Using Manually specified correlation file.")
  }else{
    message("No corelation override specified.")
  }
  
  #Import subdirectory
  if(!any(is.na(cmdArgImport('s')))){
    useSubdirectory = TRUE
    subdirectoryValueList = cmdArgImport('s')
    #message(paste("Using subdirectory", subdirectoryValue, "."))
    if(length(subdirectoryValueList ==1)){outputFolderName = paste(outputFolderName, subdirectoryValueList[1], "/", sep=""); subdirectoryValue = subdirectoryValueList[1]}
    
  }else{
    message("No subdirectory specified.")
    
    
    folderNames = basename(list.dirs(outputFolderName, full.names = TRUE, recursive = FALSE))
    if(length(folderNames) >0 ){
      useSubdirectory = TRUE
      subdirectoryValueList = folderNames
      subdirectoryValueList = subdirectoryValueList[!(subdirectoryValueList == "Alternates")]
      message("Folders not specified, using all folders but Alternates.")
    }else{
      message("Using main directory")
    }
    

  }
}


#                   ------- Code Body -------- 
for(i in 1:length(subdirectoryValueList)){
    outputFolderName = paste("Output/",filePrefix,"/", sep = "")
    message(paste("Using subdirectory", subdirectoryValueList[i], "."))
    if(useSubdirectory){
      outputFolderName = paste(outputFolderName, subdirectoryValueList[i], "/", sep="")
    }
    subdirectoryValue = subdirectoryValueList[i]
    
  #Load correlation file
  correlationFileLocation = paste(outputFolderName, filePrefix, subdirectoryValue, "CorrelationFile.rds", sep= "") #get the correlation file location based on prefix 
  if(useCategoricalPermulations){
    correlationFileLocation = paste(outputFolderName, filePrefix, subdirectoryValue, "PermulationsCorrelationFile", ".rds", sep= "")
  }
  if(useCorrelationOverride){                                                     #if a correlation override was specified, replace it with that
    correlationFileLocation = paste(outputFolderName, correlationOverride, sep="")                                      
  }
  correlationData = readRDS(correlationFileLocation)                              #Import the correlation data (non-permulated)
  
  if(usePermulations){                                                            #If permualtions are being used   
    if(useCategoricalPermulations){
        correlationData$P = correlationData$permP
      }else{
      if(usePermulationPValOverride){                                               #check for a location override
        permulationFileLocation = paste(outputFolderName, permulationPValOverride, sep="")                           #if so, use it 
      }else{                                                                        #if not, use the default 
        permulationFileLocation = paste(outputFolderName, filePrefix, "CombinedBinaryPermulationsPValue.rds", sep= "") #get the default location based on prefix 
      }
      permulationValues = readRDS(permulationFileLocation)                          #read the permulation file
      if(ncol(permulationValues)>1){                                              # add handling for continuous permulations having two columns
        correlationData$P = permulationValues$permpval                                         #replace the P column in the correlation data with the permulation values. This is the only column that the later function checks. 
        #correlationData$Rho = permulationValues$permstats
      }else{
        correlationData$P = permulationValues                                         #replace the P column in the correlation data with the permulation values. This is the only column that the later function checks. 
      }
      comparePermulationsTable = cbind(readRDS(correlationFileLocation), permulationValues)
    }
  }
  
  rerStats = getStat(correlationData)                                             #processes the RERs somewhat into stat values. only uses the P column, and the sign of the Rho column. 
  
  if(usePermulations == F){
    enrichmentCsvName = paste(outputFolderName, filePrefix, subdirectoryValue, "Enrichments.xlsx", sep= "") #make a filename based on the prefix and geneset
  }else{
    enrichmentCsvName = paste(outputFolderName, filePrefix, subdirectoryValue, "Permulation-Enrichments.xlsx", sep= "") #make a filename based on the prefix and geneset
  }
  if(!clusterRun){file.remove(enrichmentCsvName)}
  for(j in 1:length(gmtFileLocation)){
  #Load the gmt annotations 
  gmtAnnotations = read.gmt(gmtFileLocation[j])                                      #read the gmt file
  annotationsList = list(gmtAnnotations)                                          #reformat it into the format the next fuction expects
  enrichmentListName = sub("\\.gmt$", "", sub(".*/", "", gmtFileLocation[j])) #subset the file location to only the set's Name
  #enrichmentListName = substring(gmtFileLocation[j], 6, last = (nchar(gmtFileLocation[j]) - 4)) #Old code to make a set name based on the file name. 
  names(annotationsList) = enrichmentListName                                     #name geneset list with that name 
  
  enrichmentResult = fastwilcoxGMTall(rerStats, annotationsList, outputGeneVals = T, num.g =2) #run enrichment analysis 
  
  #This version run it on all of the genesets, regardless of the numebr f included genes; useful for order consistency. 
  #enrichmentResult = fastwilcoxGMTall(rerStats, annotationsList, outputGeneVals = T, num.g =0) #run enrichment analysis 
  
  #save the enrichment output
  if(usePermulations == F){
    enrichmentFileName = paste(outputFolderName, filePrefix, subdirectoryValue, "Enrichment-", enrichmentListName, ".rds", sep= "") #make a filename based on the prefix and geneset
  }else{
    enrichmentFileName = paste(outputFolderName, filePrefix, subdirectoryValue, "Enrichment-Permulation-", enrichmentListName, ".rds", sep= "") #make a filename based on the prefix and geneset
  }
  saveRDS(enrichmentResult, enrichmentFileName)                                   #Save the enrichment 
  gc()
  if(!i == 0){ #stops from writing two copies of the same sheet when no subdirectories. 
    if(!clusterRun){write.xlsx(enrichmentResult[[1]], file=enrichmentCsvName, sheetName=enrichmentListName, row.names=T, append = T)}
  }
  }
}

