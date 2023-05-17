# -- Libraries 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge)
library(tools)
source("Src/Reu/cmdArgImport.R")

# -- Usage:
# This text describes the purpose of the script 

# -- Command arguments list
# 'r = filePrefix'                            This is a prefix used to organize and separate files by analysis run. Always required. 
# 'v = <T or F>'                              This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# 'm = mainTreeFilename.txt or .rds'          This sets the location of the maintrees file
# 'a = "annotCollumn"'                        This is the column in the manual annotations spreadsheet to use
# 'c = <c("nameOfCategory1,"nameOfCategory2")>   This is the list of category names 
# 's = "screenCollumn" '                         This is a collumn which must have a value of 1 for the species to be included. 
#----------------
args = c('r=categoryTest', 'a=Meyer.Lab.Classification', 'c=c("Carnivore","Generalist","Herbivore","Insectivore","Omnivore","Piscivore")', 'm=data/RemadeTreesAllZoonomiaSpecies.rds', 'v=T')
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
mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"
annotColumn = NULL
categoryList = NULL
useScreen = F
screenColumn = NULL

{ # Bracket used for collapsing purposes
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
  
  #Annots Column
  if(!is.na(cmdArgImport('a'))){
    annotColumn = cmdArgImport('a')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY ANNOTATION COLLUMN")
  }
  
  #Category list 
  if(!is.null(cmdArgImport('c'))){
    categoryList = cmdArgImport('c')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY CATEGORIES")
  }
  
  #Screen Column
  if(!is.na(cmdArgImport('s'))){
    useScreen = T
    screenColumn = cmdArgImport('s')
  }else{
    message("No screen column used.")
  }
}


# -- Code Body --- 

manualAnnots = read.csv("Data/manualAnnotationsSheet.csv")


# - Species Filter - 
speciesFilterFilename = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="")

if(!file.exists(speciesFilterFilename) | forceUpdate){ #if no existing filter or force update, make a filter
  # --- subset the manual annots to only those with data in the categories used, and optionally by the screen column
  relevantSpecies = manualAnnots[manualAnnots[[annotColumn]] %in% categoryList,]
  if(useScreen){
    relevantSpecies = relevantSpecies[ relevantSpecies[screenColumn] %in% 1, ]
  }
  relevantSpecies = relevantSpecies[!relevantSpecies$FaName %in% "", ]
  relevantSpeciesNames = relevantSpecies$FaName
  
  saveRDS(relevantSpeciesNames, file = speciesFilterFilename)
}else{ #if not, use the existing one 
  relevantSpecieslist = readRDS(speciesFilterFilename)
  relevantSpecies = manualAnnots[ manualAnnots[["FaName"]] %in% relevantSpecieslist,]
}

# - Phenotype Vector - 
speciesNames = relevantSpecies$FaName
speciesCategories = relevantSpecies[[annotColumn]]

phenotypeVector = speciesCategories
names(phenotypeVector) = speciesNames

phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="")
saveRDS(phenotypeVector, file = phenotypeVectorFilename)
