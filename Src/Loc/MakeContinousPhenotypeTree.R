# -- Libraries 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge)
library(tools)
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/ZoonomTreeNameToCommon.R")
source("Src/Reu/ZonomNameConvertVector.R")
# -- Usage:
# This script creates a categorical tree of a phenotype which has been annotated in the Manual Annotations spreadsheet of the meyer lab. 
# In theory, this script could be used on any spreadsheet, so long as the column containing the tip.labels is named "FaName", and the column with common names is named "Common.Name.Or.Group". 

# -- Command arguments list
# r = filePrefix                                         This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                                           This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# m = mainTreeFilename.txt or .rds                       This sets the location of the maintrees file
# d = spreadSheetFilename.csv                            This sets the spreadsheet to read the data from 
# a = "annotCollumn"                                     This is the column in the manual annotations spreadsheet to use
# s = "screenCollumn"                                    This is a collumn which must have a value of 1 for the species to be included. 
# c = < "Diff" or "mean" or "last" >                                           This is used for continuous traits, to determine if the metic should be the difference between the nodes (diff), the mean(mean of the two nodes), or last(the downstream value). Note that Mean and Last are not phylogenetically independent, and do not have downstream processing. 


#----------------
args = c('r=MaturityLifespanPercent', 'm=data/newHillerMainTrees.rds', 'd=Data/MaturityLifespanData.csv', 'a=MaturityPercentage','v=T')



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
{ # Bracket used for collapsing purposes
# Defaults
mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"
spreadSheetLocation = "Data/manualAnnotationsSheet.csv"
annotColumn = NULL
categoryList = NULL
useScreen = F
screenColumn = NULL

  #MainTrees Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using default")
  }
  #read in the tree based on filetype extension
  if(file_ext(mainTreesLocation) == "rds"){
    mainTrees = readRDS(mainTreesLocation)
  }else{
    mainTrees = readTrees(mainTreesLocation) 
  }

  #spreadsheet File
  if(!is.na(cmdArgImport('d'))){
    spreadSheetLocation = cmdArgImport('d')
  }else{
    message("Using Data/manualAnnotationsSheet.csv spreadsheet")
  }
  
  #Annots Column
  if(!is.na(cmdArgImport('a'))){
    annotColumn = cmdArgImport('a')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY ANNOTATION COLLUMN")
  }
  
  #Screen Column
  if(!is.na(cmdArgImport('s'))){
    useScreen = T
    screenColumn = cmdArgImport('s')
  }else{
    message("No screen column used.")
  }

  #Continuous Metric 
  if(!is.na(cmdArgImport('c'))){
    continousMetric = cmdArgImport('c')
    if(!any(continousMetric == validMetrics)){
      stop("Contious metric is not a valid option. Use 'Diff', 'mean', or 'last'")
    }
  }else{                                                                        #If a continuous phenotype, report using Diff
      message("No continuous metric specified, using Diff")
  }
}


#                   ------- Code Body --------  

manualAnnots = read.csv(spreadSheetLocation)                      #load the manual annotations file holding the phenotype data
manualAnnots[[annotColumn]] = trimws(manualAnnots[[annotColumn]])               #trim away whitespace to allow for better matching 

# - Species Filter - 
speciesFilterFilename = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 

if(!file.exists(speciesFilterFilename) | forceUpdate){                          #if no filter exists or update is forced, make a filter 
  # --- subset the manual annots to only those with data in the categories used, and optionally by the screen column
  relevantSpecies = manualAnnots[!is.na(manualAnnots[[annotColumn]]),]          #remove all species which have no value in the specified column
  if(useScreen){                                                                #if using a screening collumn 
    relevantSpecies = relevantSpecies[ relevantSpecies[screenColumn] %in% 1, ]  #remove all species not positive for that collumn 
  }
  relevantSpecies = relevantSpecies[!relevantSpecies$FaName %in% "", ]          #remove any species without an FA name (not on the master tree)
  speciesFilter = relevantSpecies$FaName                                        #make a list of the master tree tip labels of the included species

  saveRDS(speciesFilter, file = speciesFilterFilename)                          #save that as the species filter
  
  irrelevantSpecies = manualAnnots[! manualAnnots[["FaName"]] %in% speciesFilter,]
}else{ #if not, use the existing one 
  relevantSpecieslist = readRDS(speciesFilterFilename)                          #if not, use the existing list 
  speciesFilter = relevantSpecieslist                                           #make the speciesFilter object for later 
  relevantSpecies = manualAnnots[ manualAnnots[["FaName"]] %in% relevantSpecieslist,] #and select the manual annotations entries in that list (useful if the list is more restrictive than it would be by default) 
  irrelevantSpecies = manualAnnots[! manualAnnots[["FaName"]] %in% relevantSpecieslist,]
}

# - Phenotype Vector - 
speciesNames = relevantSpecies$FaName                                           #Exract the tip name of each species
speciesValues = relevantSpecies[[annotColumn]]                                  #extract the category of each species (in same order)

phenotypeVector = speciesValues                                                 #combine those into⌄
phenotypeVector = as.numeric(phenotypeVector)
names(phenotypeVector) = speciesNames                                           #the format the functions expect


phenotypeVectorFilename = paste(outputFolderName, filePrefix, "ContinuousPhenotypeVector.rds",sep="") #make a filename based on the prefix
saveRDS(phenotypeVector, file = phenotypeVectorFilename)                        #save the phenotype vector

# - Make common name versions of objects (used in visualization) - 
commonMainTrees = mainTrees
commonMainTrees$masterTree = ZoonomTreeNameToCommon(commonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation)
commonPhenotypeVector = phenotypeVector
names(commonPhenotypeVector) = ZonomNameConvertVectorCommon(names(commonPhenotypeVector), manualAnnotLocation = spreadSheetLocation)
commonSpeciesFilter = ZonomNameConvertVectorCommon(speciesFilter, manualAnnotLocation = spreadSheetLocation)

commonPhenotypeVectorFilename = paste(outputFolderName, filePrefix, "ContinuousCommonPhenotypeVector.rds",sep="") #make a filename based on the prefix
saveRDS(commonPhenotypeVector, commonPhenotypeVectorFilename)

# - Visualized Tree - 

#No visualized tree is created for continuous traits due to the negative branches lengths involved in continuous edge values. If desired to be made, this code could be referenced from MakeCategoricalPhenotypeTree.R

# - Paths - 
pathsFilename = paste(outputFolderName, filePrefix, "ContinuousPathsFile.rds", sep= "") #make a filename based on the prefix
paths = char2Paths(phenotypeVector, mainTrees, metric = continousMetric) #make a path based on the phenotype vector
saveRDS(paths, file = pathsFilename)                                            #save the path 


?char2Paths

