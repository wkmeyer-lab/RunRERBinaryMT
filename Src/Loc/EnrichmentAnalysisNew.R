# -- Libraries 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")

# -- Usage:
# This script is used to generate enrichment values from a correlation file. Can include permulation p values if provided. 

# -- Command arguments list
# r = filePrefix                            This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                              This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# m = gmtFileLocation.gmt                   This is the location of the main gmt file
# p = <T or F>                              This sets if the code should use permulated or unpermualted values
# c = "correlationFileOverride.rds"         This can be used to manually set the correlation file. Used primarily to target pairwise comparisons of categroical phenotypes.
# f = "permulationPvalueFileLocation.rds"   This is a manual override to specify the script use a specific Permulation p-value file. 
    #If using any file other than "CombinedPrunedFastAll" with no run instance number, it must be specified manually.
#----------------
# --- Standard start-up code ---
args = commandArgs(trailingOnly = TRUE)
{  # Bracket used for collapsing purposes
  #File Prefix
  if(!is.na(cmdArgImport('r'))){
    filePrefix = cmdArgImport('r')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
  }
  
  #  Output Directory 
  if(!dir.exists("Output")){                                      #Make output directory if it does not exist
    dir.create("Output")
  }
  outputFolderNameNoSlash = paste("Output/",filePrefix, sep = "") #Set the prefix sub directory
  if(!dir.exists(outputFolderNameNoSlash)){                       #create that directory if it does not exist
    dir.create(outputFolderNameNoSlash)
  }
  outputFolderName = paste("Output/",filePrefix,"/", sep = "")
  
  #  Force update argument
  forceUpdate = FALSE
  if(!is.na(cmdArgImport('v'))){                                 #Import if update being forced with argument 
    forceUpdate = cmdArgImport('v')
    forceUpdate = as.logical(forceUpdate)
  }else{
    paste("Force update not specified, not forcing update")
  }
}

# --- Argument Imports ---
# Defaults
gmtFileLocation = "Data/enrichmentGmtFile.gmt"
usePermulations = TRUE
usePermulationPValOverride = FALSE
permulationPValOverride = NULL 
useCorrelationOverride = FALSE
correlationOverride = NULL

{ # Bracket used for collapsing purposes
  #Gmt file location
  if(!is.na(cmdArgImport('m'))){
    gmtFileLocation = cmdArgImport('m')
  }else{
    paste("No gmt location arugment, using Data/enrichmentGmtFile.gmt")                          #Report using default
    message("No gmt location arugment, using Data/enrichmentGmtFile.gmt")
  }
  
  #Permulation use
  if(!is.na(cmdArgImport('p'))){
    usePermulations = cmdArgImport('p')
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
    if(!is.na(cmdArgImport('f'))){
      usePermulationPValOverride = TRUE
      permulationPValOverride = cmdArgImport('f')
      message("Using Manually specified permulation p-value file.")
    }else{
      message("No Permulation pValue override, using standard CombinedPrunedFastAllPermulationsPValue.rds.")
    }
  }
  #Import correlation file override
  if(!is.na(cmdArgImport('f'))){
    useCorrelationOverride = TRUE
    correlationOverride = cmdArgImport('f')
    message("Using Manually specified correlation file.")
  }else{
    message("No corelation override specified.")
  }
}


#                   ------- Code Body -------- 

#Load correlation file
if(useCorrelationOverride){
  correlationFileLocation = correlationOverride
}
correlationFileLocation = paste(outputFolderName, filePrefix, "CorrelationFile.rds", sep= "")
correlationData = readRDS(correlationFileLocation)                            #Import the correlation data (non-permulated)

if(usePermulations){
  if(usePermulationPValOverride){
    permulationFileLocation = permulationPValOverride
  }else{
    permulationFileLocation = paste(outputFolderName, filePrefix, "CombinedPrunedFastAllPermulationsPValue.rds", sep= "")
  }
  permulationValues = readRDS(permulationFileLocation)
  correlationData$P = permulationValues
}

rerStats = getStat(correlationData)

#Load the gmt annotations 
gmtAnnotations = read.gmt(gmtFileLocation)
annotationsList = list(gmtAnnotations)
enrichmentListName = substring(gmtFileLocation, 6, last = (nchar(gmtFileLocation) - 4)) #determine geneset name from filename 
names(annotationsList) = enrichmentListName

enrichmentResult = fastwilcoxGMTall(rerStats, annotationsList, outputGeneVals = T, num.g =10)

#save the enrichment output
enrichmentFileName = paste(outputFolderName, filePrefix, "Enrichment-", enrichmentListName, ".rds", sep= "")
saveRDS(enrichmentResult, enrichmentFileName)



# --- Visualize the enrichment ----

{
  #This is manual only -- run-as-script does not accept a visualize output because no way to output result. 
  #For a script version, use PvQvGoVisualize.R 
  visualize = T
  visualize = F
  clean = T
  clean = F
  
  if(visualize){
    library(stringr)
    library(insight)
    enrichmentResult2 = enrichmentResult[[1]]
    makeGOTable = function(data, collumn){
      ValueHead = head(data[order(collumn, decreasing = T),], n=40)
      ValueHead$num.genes = as.character(ValueHead$num.genes)
      ValueHead$stat = round(ValueHead$stat, digits = 5)
      ValueHead$stat = as.character(ValueHead$stat)
      ValueHead = format_table(ValueHead, pretty_names = F, digits = "scientific5")
      ValueHead
    }
    enrichHead = makeGOTable(enrichmentResult2, abs(enrichmentResult2$stat))
    enrichHead
    textplot(enrichHead[1:4], mar = c(0,0,2,0), cmar = 1.5)
    if(usePermulations){
      title(main = paste("Top pathways by permulation"))
    }else{
      title(main = paste("Top pathways by non-permulation"))
    }
    
    if(clean){
      # ---- enrichment cleaning ----
      pathwaysToRemove = grep("CANCER", rownames(enrichHead))
      rowsToKeep = (!1:nrow(enrichHead) %in% pathwaysToRemove)
      cleanedHead = enrichHead[rowsToKeep,]
      textplot(cleanedHead[1:4], mar = c(0,0,2,0), cmar = 1.5)
      if(usePermulations){
        title(main = paste("Top pathways by permulation"))
      }else{
        title(main = paste("Top pathways by non-permulation"))
      }
      CleanheadFileName = paste(outputFolderName, filePrefix, "CleanedEnrichmentHead.csv", sep= "")
      write.csv(cleanedHead, CleanheadFileName)
    }
  }
}  
