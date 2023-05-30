# -- Libraries 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge)
library(tools)
source("Src/Reu/cmdArgImport.R")

# -- Usage:
# This text describes the purpose of the script 

# -- Command arguments list
# r = filePrefix    This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>      This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# m = mainTreeFilename.txt or .rds    This sets the location of the maintrees file
# c = <T or F>                        This is used to set if the script is targeting metacombined permulations. Called "metacombination". Used for parrallelization. 


#----------------
args = c('r=CategoricalDiet', 'm=data/RemadeTreesAllZoonomiaSpecies.rds') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.
permulationsDataDemo = readRDS("Output/CategoricalDiet3Phen/CategoricalDiet3PhenPermulationsData1.rds")
permulationsDataDemoToy = permulationsDataDemo[1:30]
# --- Standard start-up code ---
args = commandArgs(trailingOnly = TRUE)
{  # Bracket used for collapsing purposes
  #File Prefix
  if(!is.na(cmdArgImport('r'))){                                                #This cmdArgImport script is a way to import arguments from the command line. 
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
mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"  #this is the standard location on Sol 
metacombineValue = FALSE
{ # Bracket used for collapsing purposes
  #MainTrees Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using default")
  }
  #Metacombination 
  if(!is.na(cmdArgImport('c'))){
    metacombineValue = cmdArgImport('c')
    metacombineValue = as.logical(metacombineValue)
    if(is.na(metacombineValue)){
      metacombineValue = FALSE
      message("Metacombination value not interpretable as logical. Did you remember to capitalize? Using FALSE.")
    }
  }else{
    message("Metacombination value not specified, using FALSE. If you aren't parrallelizing, don't worry about this.")
  }
  
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
if(metacombineValue == F){
  combinedDataFileName = paste(outputFolderName, filePrefix, "Combined", permulationPrefix,"PermulationsData", runInstanceValue, ".rds", sep="")
}else{
  combinedDataFileName = paste(outputFolderName, filePrefix, "MetaCombined", permulationPrefix, "PermulationsData", runInstanceValue, ".rds", sep="")
}
permulationsData = readRDS(combinedDataFileName)


# -- calculate p values -- 

pValues = getPermPvalsCategoricalReport(correlationsObject, permulationsData$trees[1:30], phenotypeVector, mainTrees, RERObject)

permulationOutputFilename = correlationFileName = paste(outputFolderName, filePrefix, "CategroicalPermulationPValuesFile.rds", sep= "")
saveRDS(pValues, permulationOutputFilename)
