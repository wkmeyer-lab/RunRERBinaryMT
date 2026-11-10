# -- Libraries 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #DO NOT RUN LOCALLY. add Cluster path to custom libraries to searched locations

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


#-------------------------------
# --- Command arguments list ---
# - Required - 
# r = filePrefix                            REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# m = mainTreeFilename.txt or rds           This is the location of the main trees file 
# d = spreadSheetFilename.csv               This is the location of the annotation csv. If not the default "manualAnnotationsSheet.csv", some functions in other scripts may break. 
# a = annotCollumn                          This is the name of the column contianing the phenotype data 
# t = tipColumn                             This sets the column with the tip names as they appear in the maintrees file.


# - Required, Defaults likely acceptable -
# v = <T or F>                              This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# n = "nameColumn"                          This sets the column with common names of the species. 
# - Optional - 
# h = <T>                                   This ARG causes the script to return its ARG list, instead of running the script. 
# s = "screenColumn"                        This sets the column with the screening data
# n = <"uni" OR "bi">                       This set unidirectional or bidirectional transition 
# c = <"ancestral" OR "all" OR "terminal">  This sets the clade type to be used in tree creation
# w = <T or F>                              This sets if the tree should be weighted 
# - Pruning -
# z = <minimum branch length>               This sets the minimum branch length for terminal branches in the master tree. Branches shorter than this will be removed. 
# x = "pruningPrefrenceColumn"              This sets a column, where if the value is 1, the tip will be preferentially kept. If the value is TRUE, the tip will never be pruned.
# y = "c('unprunedtip1', 'unprunedtip2')"   This allows you to add a list of specific tips to not be dropped during pruning. Must use the tip name, not common name. 
# - Overrides - 
# --- End arguments list ---
#-------------------------------


# --- Standard start-up code ---
if(clusterRun){args = commandArgs(trailingOnly = TRUE)}
{  # Bracket used for collapsing purposes
  #Arg help 
  if(isTRUE(as.logical(cmdArgImport('h')))){
    message("Not Running Script, instead returning ARGs for script.")
    printArgDocumentation()
    stop()
  }
  
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
    message("Force update not specified, not forcing update")
  }
}

# --- Argument Imports ---
{ # Bracket used for collapsing purposes
  # Defaults
  mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"
  annotationsLocation = "Data/manualAnnotationsSheet.csv"
  phenotypeColumn = "ERRORDEFAULT"
  transitionValue = "Default"
  cladeValue = "Default"
  weightValue = FALSE
  useScreen = F
  screenCollumn = NA
  nameCollumn = NULL
  tipCol = NULL
  usingPruning = F
  pruningCutoff = NA
  pruningPrefrenceColumn = NA
  
  #Main Tree Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using default")
  }
  
  #Annotations Location
  if(!is.na(cmdArgImport('d'))){
    annotationsLocation = cmdArgImport('d')
  }else{
    message("No maintrees arg, using default 'manualAnnotationsSheet.csv' ")
  }
  
  #Phenotype Column
  if(!is.na(cmdArgImport('a'))){
    phenotypeColumn = cmdArgImport('a')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY PHENOTYPE COLLUMN")
  }
  
  #Transition value
  if(!is.na(cmdArgImport('n'))){
    transitionValue = cmdArgImport('n')
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
  
  #Screen Column 
  if(!is.na(cmdArgImport('s'))){
    useScreen = T
    screenCollumn = cmdArgImport('s')
  }else{
    message("Screen Column not specified, not using screen column.")
  }

  #Tip Column
  if(!is.na(cmdArgImport('t'))){
    tipCol = cmdArgImport('t')
  }else{
    message("Tip Column not specified, using 'tipName'.")
  }
  #Name Column
  if(!is.na(cmdArgImport('n'))){
    nameColumn = cmdArgImport('n')
  }else{
    message("Name Column not specified, using default list.")
  }
  
  #Pruning cutoff
  if(!is.na(cmdArgImport('z'))){
    usingPruning = T
    pruningCutoff = cmdArgImport('z')
  }else{
    message("Pruning Cutoff not specified, not pruning tree.")
  }
  
  #PruningPrefrenceColumn 
  if(!is.na(cmdArgImport('x'))){
    pruningPrefrenceColumn = cmdArgImport('x')
  }else{
    if(usingPruning){message("No pruning prefrence column specified")}
  }
  
  #ManualPruningProtections
  if(!is.na(cmdArgImport('y'))){
    manualPruningProtections = cmdArgImport('y')
  }else{
    if(usingPruning){message("No manually protected species specified")}
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

if(!file.exists(speciesFilterFilename) | forceUpdate){                          #If not using exisitng species filter
  relevantSpecies = manualAnnots[ manualAnnots[[phenotypeColumn]] %in% c(0,1),] #Set the relevant species to species with a 1 or a 0 in the annotation column 
  if(useScreen){                                                                #if using a screening column, 
    relevantSpecies = relevantSpecies[relevantSpecies[screenCollumn] == 1,]     #only include species with a 1 in the screen column
  }                                                                             
  relevantSpeciesNames = relevantSpecies[[tipCol]]                                 #use the species names 
  relevantNamesNA = which(is.na(relevantSpeciesNames))
  if(length(relevantNamesNA) > 0){  relevantSpeciesNames = relevantSpeciesNames[-relevantNamesNA]}

  
  if(usingPruning){
    source("Src/Reu/autoPruner.R")
    pruningProtectionSpecies = NA
    if(!is.na(pruningPrefrenceColumn)){
      pruningProtectionRows = manualAnnots[which(as.logical(manualAnnots[[pruningPrefrenceColumn]])),]
      pruningProtectionSpecies = pruningProtectionRows[[tipCol]]
      if(all(is.logical(manualAnnots[[pruningPrefrenceColumn]]))){
        pruningProtection = T
      }else{ 
        pruningProtection = F
      }
    }
    allProtectedSpecies = append(pruningProtectionSpecies, manualPruningProtections)
    
    
    workingTree = mainTrees$masterTree
    workingTree = drop.tip(workingTree, which(!workingTree$tip.label %in% relevantSpeciesNames))
    
    fewGeneSpecies = dropFewGeneSpecies(mainTrees, workingTree, nameConversionColumn = tipCol, nameConversionData = annotationsLocation)
    workingTree = drop.tip(workingTree, fewGeneSpecies)
    
    pruningFilename = paste(outputFolderName, filePrefix, "PruningTree.pdf", sep="")
    pdf(pruningFilename, width=16, height = 14)
    prunedTree = autopruner(workingTree, dropValue = pruningCutoff, tipsToKeep = pruningProtectionSpecies, nameConversionColumn = tipCol, nameConversionData = annotationsLocation, preDroppedTips = fewGeneSpecies)
    if(!pruningProtection){
      prunedTree = autopruner(prunedTree, dropValue = pruningCutoff, tipsToKeep = manualPruningProtections, nameConversionColumn = tipCol, nameConversionData = annotationsLocation, preDroppedTips = droppedTips, originalTree = workingTree)
    }
    dev.off()
    
    prunedSpecies = relevantSpeciesNames[!relevantSpeciesNames %in% prunedTree$tip.label]
    relevantSpeciesNames = relevantSpeciesNames[-which(relevantSpeciesNames %in% prunedSpecies)]
  }
  
  saveRDS(relevantSpeciesNames, file = speciesFilterFilename)
  
}else{
  relevantSpecieslist = readRDS(speciesFilterFilename)
  relevantSpeciesNames = relevantSpecieslist
  relevantSpecies = manualAnnots[ manualAnnots[[tipCol]] %in% relevantSpecieslist,]
}

# - Setup foreground Species --

foregroundSpeciesAnnot = relevantSpecies[ relevantSpecies[[phenotypeColumn]] %in% 1,] #set the foreground species to spcies with a 1 in the annotation column

foregroundNames = foregroundSpeciesAnnot[[tipCol]]
foregroundNamesNA = which(is.na(foregroundNames))
if(length(foregroundNamesNA) > 0){  foregroundNames = foregroundNames[-foregroundNamesNA]}


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


# ---- Saving and visualization ----

# - Save the Tree - 
binaryTreeFilename = paste(outputFolderName, filePrefix, "BinaryForegroundTree.rds", sep="")
saveRDS(binaryForegroundTreeOutput, file = binaryTreeFilename)
binaryTreeFilename = paste(outputFolderName, filePrefix, "BinaryTree.rds", sep="")
saveRDS(binaryForegroundTreeOutput, file = binaryTreeFilename)
binaryNewickFilename = paste(outputFolderName, filePrefix, "NewickTree.nwk", sep="")
write.tree(binaryForegroundTreeOutput, file = binaryNewickFilename)

# - Read back in tree and print to pdf - 
readTest = readRDS(binaryTreeFilename)
testTreeDisplayable = readTest
testTreeDisplayable$edge.length = replace(testTreeDisplayable$edge.length, testTreeDisplayable$edge.length==0, 0.5)
testTreeDisplayable$edge.length = replace(testTreeDisplayable$edge.length, testTreeDisplayable$edge.length==1, 4)

#source("Src/Reu/plotBinaryTree.R")

binaryTreePdfname = paste(outputFolderName, filePrefix, "BinaryForegroundTree.pdf", sep="")
#pdf(binaryTreePdfname, width=8, height = 14)
#plotBinaryTree(mainTrees, readTest, foregroundNames, mainTitle = paste(filePrefix, "Binary", "Foreground", "Tree"), tipColumn = tipCol)
#plotBinaryTree(mainTrees, readTest, foregroundNames, convertNames = F, mainTitle = paste(filePrefix, "Binary", "Foreground", "Tree"), tipColumn = tipCol)
plotTree(testTreeDisplayable)
#dev.off()


