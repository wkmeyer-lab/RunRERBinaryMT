# -- Libraries 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge)
library(tools)
source("Src/Reu/cmdArgImport.R")

# -- Usage:
# This Script creates a binary tree from a csv file containing (at least) the following columns: 
# FaName: The name of the species as it appears in the Mutliphylo 
# Common.Name.Or.Group: The common name of the species
# Species.Name: The scientific name of the species
# [phenotypeCollumn]: A column (name set by argument) which contains the pehnotype data 
# (optional)[ScreenColumn]: A column which can be used to select a subset of the data. If a screen column is specified by argument, only species with a "1" in that column will be used to make the tree. 
# If this file is not named "manualAnnotationsSheet.csv", some functions may break. 

# -- Command arguments list
# r = filePrefix                           This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                             This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# m = mainTreeFilename.txt or rds          This is the location of the main trees file 
# a = annotationFileLocation.csv           This is the location of the annotation csv. If not the default "manualAnnotationsSheet.csv", some functions in other scripts may break. 
# p = phenotypeColumn                      This is the name of the column contianing the phenotype data 
# t = <"uni" OR "bi">                      This set unidirectional or bidirectional transition 
# c = <"ancestral" OR "all" OR "terminal"> This sets the clade type to be used in tree creation
# w = <T or F>                             This sets if the tree should be weighted 
# s = "screenColumn"                       This sets the column with the screening data


#----------------
args = c('m=data/RemadeTreesAllZoonomiaSpecies.rds', "r=EcholocationUpdate2", "t=bi", "p=Echolocation", "c=all", "v=T", "s=Laurasiatheria")
args = c("m=data/FirstExpressionTrees.rds", "r=LiverExpression", "p=Carnivory", "t=bi", "c=all", "w=F", "v=T", "a=Data/ExpressionAnnots.csv")
args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=CVHApplesToApples", "p=Carnivory", "t=bi", "c=all", "w=F", "v=T", "a=Data/ExpressionAnnots.csv")

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
mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"
annotationsLocation = "Data/manualAnnotationsSheet.csv"
phenotypeColumn = "ERRORDEFAULT"
transitionValue = "Default"
cladeValue = "Default"
weightValue = FALSE
useScreen = F
screenCollumn = NA

{ # Bracket used for collapsing purposes
  
  #Main Tree Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using default")
  }
  
  #Annotations Location
  if(!is.na(cmdArgImport('a'))){
    annotationsLocation = cmdArgImport('a')
  }else{
    message("No maintrees arg, using default 'manualAnnotationsSheet.csv' ")
  }
  
  #Phenotype Column
  if(!is.na(cmdArgImport('p'))){
    phenotypeColumn = cmdArgImport('p')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY PHENOTYPE COLLUMN")
  }
  
  #Transition value
  if(!is.na(cmdArgImport('t'))){
    transitionValue = cmdArgImport('t')
  }else{
    message("Using default bidirectional transistion")
  }
  
  #Clade value
  if(!is.na(cmdArgImport('c'))){
    cladeValue = cmdArgImport('c')
  }else{
    message("Using default all clade")
  }
  
  #Weight Value
  if(!is.na(cmdArgImport('w'))){
    weightValue = cmdArgImport('w')
  }else{
    message("Weight = false")
  }
  
  if(!is.na(cmdArgImport('s'))){
    useScreen = T
    screenCollumn = cmdArgImport('s')
  }else{
    message("Screen Column not specified, not using screen column.")
  }
}

#                   ------- Code Body -------- 

# - Import Files -
if(file_ext(mainTreesLocation) == "rds"){
  mainTrees = readRDS(mainTreesLocation)
}else{
  mainTrees = readTrees(mainTreesLocation) 
}
manualAnnots = read.csv(annotationsLocation)

# - Species Filter - 
speciesFilterFilename = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="")

if(!file.exists(speciesFilterFilename) | forceUpdate){
  relevantSpecies = manualAnnots[ manualAnnots[[phenotypeColumn]] %in% c(0,1),]
  if(useScreen){
    relevantSpecies = relevantSpecies[relevantSpecies[screenCollumn] == 1,]
  }
  relevantSpeciesNames = relevantSpecies$FaName
  saveRDS(relevantSpeciesNames, file = speciesFilterFilename)
  
}else{
  relevantSpecieslist = readRDS(speciesFilterFilename)
  relevantSpeciesNames = relevantSpecieslist
  relevantSpecies = manualAnnots[ manualAnnots[["FaName"]] %in% relevantSpecieslist,]
}

# - Setup foreground Species --

foregroundSpeciesAnnot = relevantSpecies[ relevantSpecies[[phenotypeColumn]] %in% 1,]

foregroundNames = foregroundSpeciesAnnot$FaName

foregroundFilename = paste(outputFolderName, filePrefix, "BinaryTreeForegroundSpecies.rds", sep="")
saveRDS(foregroundNames, file = foregroundFilename)

# -- set arguments for foreground2Trees --
f2tInputList = list(foregroundNames, mainTrees, useSpecies = relevantSpeciesNames)
#Transition
if(transitionValue == "uni"){
  f2tInputList[["transition"]]= "unidirectional"
  message("Unidirectional transition")
}else{
  f2tInputList[["transition"]]= "bidirectional"
  message("Bidirectional transition")
}

#clade
if(cladeValue == "ancestral"){
  f2tInputList[["clade"]] = "ancestral"
  message("ancestral clade")
}else if(cladeValue == "terminal"){
  f2tInputList[["clade"]] = "terminal"
  message("terminal clade")
}else{
  f2tInputList[["clade"]] = "all"
  message("all clade")
}

#Weight
if(weightValue ==TRUE || weightValue == 't'){
  f2tInputList[["weighted"]] = TRUE
  message("Weighted = TRUE")
}else{
  f2tInputList[["weighted"]] = FALSE
  message("Weighted = False")
}

# - Make the tree - 
binaryForegroundTreeOutput = do.call(foreground2Tree, f2tInputList)



# - Save the Tree - 
binaryTreeFilename = paste(outputFolderName, filePrefix, "BinaryForegroundTree.rds", sep="")
saveRDS(binaryForegroundTreeOutput, file = binaryTreeFilename)
binaryTreeFilename = paste(outputFolderName, filePrefix, "BinaryTree.rds", sep="")
saveRDS(binaryForegroundTreeOutput, file = binaryTreeFilename)

# - Read back in tree and print to pdf - 
readTest = readRDS(binaryTreeFilename)
testTreeDisplayable = readTest
testTreeDisplayable$edge.length = replace(testTreeDisplayable$edge.length, testTreeDisplayable$edge.length==0, 0.5)
testTreeDisplayable$edge.length = replace(testTreeDisplayable$edge.length, testTreeDisplayable$edge.length==1, 4)

source("Src/Reu/plotBinaryTree.R")

binaryTreePdfname = paste(outputFolderName, filePrefix, "BinaryForegroundTree.pdf", sep="")
pdf(binaryTreePdfname, width=8, height = 14)
plotBinaryTree(mainTrees, readTest, foregroundNames, mainTitle = paste(filePrefix, "Binary", "Foreground", "Tree"))
plotBinaryTree(mainTrees, readTest, foregroundNames, convertNames = F, mainTitle = paste(filePrefix, "Binary", "Foreground", "Tree"))
#plotTreeHighlightBranches(testTreeDisplayable, hlspecies=which(readTest$edge.length==1), hlcols="blue",)
dev.off()


