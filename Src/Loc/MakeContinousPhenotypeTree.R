# -- Libraries 
runningLocally =identical(Sys.getenv("RSTUDIO"), "1")
clusterRun = F
if(!runningLocally){clusterRun = T}
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations

library(RERconverge)
library(tools)
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/ZonomNameConvertVectorCommon.R")
# -- Usage:
# This script creates a categorical tree of a phenotype which has been annotated in the Manual Annotations spreadsheet of the meyer lab. 
# In theory, this script could be used on any spreadsheet, so long as the column containing the tip.labels is named "FaName", and the column with common names is named "Common.Name.Or.Group". 

#-------------------------------
# --- Command arguments list ---
# - Required - 
# r = filePrefix                            REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# m = mainTreeFilename.txt or .rds          This sets the location of the maintrees file
# d = spreadSheetFilename.csv               This sets the spreadsheet to read the data from 
# a = "annotCollumn"                        This is the column in the manual annotations spreadsheet to use
# t = "tipColumn"                           This sets the column with the tip names as they appear in the maintrees file.
# n = "nameColumn"                          This sets the column with common names of the species. 
# - Required, Defaults likely acceptable -
# v = <T or F>                              This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# k = < "Diff" or "mean" or "last" >        This is used for continuous traits, to determine if the metic should be the difference between the nodes (diff), the mean(mean of the two nodes), or last(the downstream value). Note that Mean and Last are not phylogenetically independent, and do not have downstream processing. 
# - Optional - 
# h = <T>                                   This ARG causes the script to return its ARG list, instead of running the script. 
# s = "screenCollumn"                       This is a collumn which must have a value of 1 for the species to be included. 
# - Pruning - 
# z = <minimum branch length>                            This sets the minimum branch length for terminal branches in the master tree. Branches shorter than this will be removed. 
# x = "pruningPrefrenceColumn"                           This sets a column, where if the value is 1, the tip will be preferentially kept. If the value is TRUE, the tip will never be pruned.
# y = "c('unprunedtip1', 'unprunedtip2')"                This allows you to add a list of specific tips to not be dropped during pruning. Must use the tip name, not common name. 
# p = "c('prunedtip1', 'prunedtip2')"                    This allows you to manually specify additional branches to be pruned
# - Overrides - 
# --- End arguments list ---
#-------------------------------




# --- Standard start-up code ---
if(clusterRun)args = commandArgs(trailingOnly = TRUE)
{  # Bracket used for collapsing purposes
  #Arg help 
  if(isTRUE(as.logical(cmdArgImport('h')))){
    message("Not Running Script, instead returning ARGs for script.")
    printArgDocumentation()
    stop()
  }
  
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
    message("Force update not specified, not forcing update")
  }
}

# --- Argument Imports ---
{ # Bracket used for collapsing purposes
# Defaults
mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"
spreadSheetLocation = "Data/mergedData.csv"
annotColumn = NULL
categoryList = NULL
useScreen = F
screenColumn = NULL
nameColumn = NULL
tipCol = NULL
continousMetric = "diff"
validMetrics= c("diff", "mean", "last")
usingPruning = F
usingAutoPruning = F
manualPruningProtections = NULL
pruningPrefrenceColumn = NA
pruningProtection = F
manualPrunedSpecies = NULL

  #MainTrees Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using default")
  }
  #read in the tree based on filetype extension
  if(file_ext(mainTreesLocation) == "rds"){
    if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
  }else{
    if(!exists("mainTrees")){mainTrees = readTrees(mainTreesLocation)} 
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
  if(!is.na(cmdArgImport('k'))){
    continousMetric = cmdArgImport('k')
    if(!any(continousMetric == validMetrics)){
      stop("Contious metric is not a valid option. Use 'Diff', 'mean', or 'last'")
    }
  }else{                                                                        #If a continuous phenotype, report using Diff
      message("No continuous metric specified, using Diff")
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
    usingAutoPruning = T
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
  if(!all(is.na(cmdArgImport('y')))){
    manualPruningProtections = cmdArgImport('y')
  }else{
    if(usingPruning){message("No manually protected species specified")}
  }
  
  #ManualPruningSpecies
  if(!all(is.na(cmdArgImport('p')))){
    manualPruningSpecies = cmdArgImport('p')
    usingPruning = T
  }else{
    if(usingPruning){message("No manually pruned species specified")}
  }
}


#                   ------- Code Body --------  

manualAnnots = read.csv(spreadSheetLocation)                      #load the manual annotations file holding the phenotype data
manualAnnots[[annotColumn]] = trimws(manualAnnots[[annotColumn]])               #trim away whitespace to allow for better matching 

# - Species Filter - 
speciesFilterFilename = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 

if(!file.exists(speciesFilterFilename) | forceUpdate){                          #if no filter exists or update is forced, make a filter 
  # --- subset the manual annots to only those with data in the categories used, and optionally by the screen column
  relevantSpecies = manualAnnots[!is.na(manualAnnots[[annotColumn]]),]#remove all species which are not part of the specified categories
  if(useScreen){                                                                #if using a screening collumn 
    relevantSpecies = relevantSpecies[ relevantSpecies[[screenColumn]] %in% 1, ]  #remove all species not positive for that collumn 
  }
  relevantSpecies = relevantSpecies[!relevantSpecies[[tipCol]] %in% "", ]          #remove any species without an FA name (not on the master tree)
  speciesFilter = relevantSpecies[[tipCol]]                                       #make a list of the master tree tip labels of the included species
  if(usingPruning){
    if(usingAutoPruning){
      source("Src/Reu/autoPruner.R")
      pruningProtectionSpecies = NULL
      if(!is.na(pruningPrefrenceColumn)){
        if(all(is.logical(manualAnnots[[pruningPrefrenceColumn]]))){
          pruningProtection = T
        }else{ 
          pruningProtection = F
        }
        
        pruningProtectionRows = manualAnnots[which(as.logical(manualAnnots[[pruningPrefrenceColumn]])),]
        pruningProtectionSpecies = pruningProtectionRows[[tipCol]]
      }
      allProtectedSpecies = append(pruningProtectionSpecies, manualPruningProtections)
      
      workingTree = mainTrees$masterTree
      workingTree = drop.tip(workingTree, which(!workingTree$tip.label %in% speciesFilter))
      
      fewGeneSpecies = dropFewGeneSpecies(mainTrees, workingTree, nameConversionColumn = tipCol, nameConversionData = spreadSheetLocation)
      fewGeneSpecies = fewGeneSpecies[- which(fewGeneSpecies %in% allProtectedSpecies)]
      workingTree = drop.tip(workingTree, fewGeneSpecies)
      names(fewGeneSpecies)[1:length(fewGeneSpecies)] = "fewGenes"
      
      pruningFilename = paste(outputFolderName, filePrefix, "PruningTree.pdf", sep="")
      pdf(pruningFilename, width = 16, height = length(workingTree$tip.label)/8)
      prunedTree = autopruner(workingTree, dropValue = pruningCutoff, tipsToKeep = allProtectedSpecies, nameConversionColumn = tipCol, nameConversionData = spreadSheetLocation, preDroppedTips = fewGeneSpecies)
      if(!pruningProtection){
        prunedTree = autopruner(prunedTree, dropValue = pruningCutoff, tipsToKeep = manualPruningProtections, nameConversionColumn = tipCol, nameConversionData = spreadSheetLocation, preDroppedTips = droppedTips, originalTree = workingTree)
      }
      dev.off()
    }
    if(!is.null(all(manualPruningSpecies))){
      prunedTree = drop.tip(prunedTree, manualPruningSpecies)
      names(manualPruningSpecies)[1:length(manualPruningSpecies)] = "manualDrop"
      droppedTips = append(droppedTips, manualPruningSpecies)
    }
    
    
    prunedSpecies = speciesFilter[!speciesFilter %in% prunedTree$tip.label]
    speciesFilter = speciesFilter[-which(speciesFilter %in% prunedSpecies)]
    
    prunedSpeciesFilename = paste(outputFolderName, filePrefix, "prunedSpecies.rds",sep="")
    saveRDS(droppedTips, prunedSpeciesFilename)
    prunedSpeciesTextFilename = file(paste(outputFolderName, filePrefix, "prunedSpecies.txt",sep=""))
    writeLines(print(droppedTips),prunedSpeciesTextFilename)
    close(prunedSpeciesTextFilename)
  }
  
  
  
  saveRDS(speciesFilter, file = speciesFilterFilename)                          #save that as the species filter
  
  irrelevantSpecies = manualAnnots[! manualAnnots[[tipCol]] %in% speciesFilter,]
}else{ #if not, use the existing one 
  relevantSpecieslist = readRDS(speciesFilterFilename)                          #if not, use the existing list 
  speciesFilter = relevantSpecieslist                                           #make the speciesFilter object for later 
  relevantSpecies = manualAnnots[ manualAnnots[[tipCol]] %in% relevantSpecieslist,] #and select the manual annotations entries in that list (useful if the list is more restrictive than it would be by default) 
  irrelevantSpecies = manualAnnots[! manualAnnots[[tipCol]] %in% relevantSpecieslist,]
}

# - Phenotype Vector - 
relevantSpecies = relevantSpecies[relevantSpecies[[tipCol]] %in% speciesFilter, ] 
speciesValues = relevantSpecies[[annotColumn]]                                  #extract the category of each species (in same order)
speciesNames = relevantSpecies[[tipCol]]                                         #Exract the tip name of each species

phenotypeVector = speciesValues                                                 #combine those into⌄
phenotypeVector = as.numeric(phenotypeVector)
names(phenotypeVector) = speciesNames                                           #the format the functions expect

#remove species with an NA phenotype
speciesMissingPhenotypeData = which(is.na(phenotypeVector))
if(length(speciesMissingPhenotypeData) >0){
  phenotypeVector = phenotypeVector[-speciesMissingPhenotypeData]
}


phenotypeVectorFilename = paste(outputFolderName, filePrefix, "ContinuousPhenotypeVector.rds",sep="") #make a filename based on the prefix
saveRDS(phenotypeVector, file = phenotypeVectorFilename)                        #save the phenotype vector

# - Make common name versions of objects (used in visualization) - 
commonMainTrees = mainTrees
commonMainTrees$masterTree$tip.label = ZonomNameConvertVectorCommon(commonMainTrees$masterTree$tip.label, annotationLocation = spreadSheetLocation, tipColumn = tipCol, commonColumn = nameColumn)
commonPhenotypeVector = phenotypeVector
names(commonPhenotypeVector) = ZonomNameConvertVectorCommon(names(commonPhenotypeVector), annotationLocation = spreadSheetLocation, tipColumn= tipCol, commonColumn = nameColumn)
commonSpeciesFilter = ZonomNameConvertVectorCommon(speciesFilter, annotationLocation = spreadSheetLocation, tipColumn= tipCol, commonColumn = nameColumn)

commonPhenotypeVectorFilename = paste(outputFolderName, filePrefix, "ContinuousCommonPhenotypeVector.rds",sep="") #make a filename based on the prefix
saveRDS(commonPhenotypeVector, commonPhenotypeVectorFilename)

# - Paths - 
pathsFilename = paste(outputFolderName, filePrefix, "ContinuousPathsFile.rds", sep= "") #make a filename based on the prefix
paths = char2Paths(phenotypeVector, mainTrees, metric = continousMetric) #make a path based on the phenotype vector
saveRDS(paths, file = pathsFilename)                                            #save the path 

# - Visualized Tree - 

#No visualized tree is created for continuous traits due to the negative branches lengths involved in continuous edge values. If desired to be made, this code could be referenced from MakeCategoricalPhenotypeTree.R
  #This could be done by plotting the master tree, wit the continuous values as edge labels, but that requires having the phenotype values on each branch. This ould be extracted from the paths object, but I don't know how. 
