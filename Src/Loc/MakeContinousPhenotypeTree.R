# -- Usage:
# This script creates a phenotype tree of a phenotype with a continuous value. 
# If you wish to only include select species in the analysis, add a new column to your data as a "screen column". Set species to include to have a 1 in that column, and other species to 0. Use the 's' arg to specify the screen column's name.
# In practice, this function actually creates a PATHs file for a continuous phenotype (the next step of processing). 
# This is because in practice continuous values can be negative if using the 'diff' setting, which produces weird results on a tree, and all functions actually use the PATHs, not the tree. 
# The script is named this way to retain parity with the other Make-X-Tree functions; this is the continuous equivalent of such. 

# -- Libraries 
source("Src/Reu/loadLibraries.R")
loadLibraries() # Note that this library load includes cluster run detection. 
source("Src/Reu/ZonomNameConvertVectorCommon.R")
#-------------------------------
# --- Command arguments list ---
# - Required -
# r  = filePrefix                                         This is a prefix used to organize and separate files by analysis run. Always required. 
# m  = mainTreeFilename.txt or .rds                       This sets the location of the maintrees file
# td = spreadSheetFilename.csv                            This sets the spreadsheet to read the data from 
# tp = "phenotypeColumn"                                  This is the column with the phenotype information
# tt = "tipColumn"                                        This sets the column with the tip names as they appear in the maintrees file.

# - Required, Defaults likely acceptable -
# f  = forceUpdate <T or F>                               This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# tn = "nameColumn"                                       This sets the column with common names of the species. 
# k  = < "Diff" or "mean" or "last" >                     This is used for continuous traits, to determine if the metic should be the difference between the nodes (diff), the mean(mean of the two nodes), or last(the downstream value). Note that Mean and Last are not phylogenetically independent, and do not have downstream processing. 

# - Optional - 
# h  = <T>                                                This ARG causes the script to return its ARG list, instead of running the script. 
# ts = "screenCollumn"                                    This is a collumn which must have a value of 1 for the species to be included. 

# - Pruning - 
# pm = "c('prunedtip1', 'prunedtip2')"                    This allows you to manually specify additional branches to be pruned
# pc = <minimum branch length>                            This sets the minimum branch length for terminal branches in the master tree. Branches shorter than this will be removed. 
# pk = "pruningPrefrenceColumn"                           This sets a column, where if the value is 1, the tip will be preferentially kept. If the value is TRUE, the tip will never be pruned.
# ps = "c('unprunedtip1', 'unprunedtip2')"                This allows you to add a list of specific tips to not be dropped during pruning. Must use the tip name, not common name. 
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
  mainTreesLocation = "Results/mainTreesFromAlignments.rds"
  spreadSheetLocation = "Data/mergedData.csv"
  phenotypeColumn = NULL
  tipCol = "tipName"
  nameColumn = NULL
  continousMetric = "diff"
  useScreen = F
  screenColumn = NULL
  usingPruning = F
  manualPruningSpecies = NULL
  pruningCutoff = NULL
  
  #Load from args 
  source("Src/Reu/ArgumentImport.R")
  
  validMetrics= c("diff", "mean", "last")
  if(!any(continousMetric == validMetrics)){
    stop("Contious metric is not a valid option. Use 'diff', 'mean', or 'last'")
  }
}


#                   ------- Code Body --------  
if(file_ext(mainTreesLocation) == "rds"){
  if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
}else{
  if(!exists("mainTrees")){mainTrees = readTrees(mainTreesLocation)} 
}


manualAnnots = read.csv(spreadSheetLocation)                      #load the manual annotations file holding the phenotype data
manualAnnots[[phenotypeColumn]] = trimws(manualAnnots[[phenotypeColumn]])               #trim away whitespace to allow for better matching 

# - Species Filter - 
speciesFilterFilename = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 

if(!file.exists(speciesFilterFilename) | forceUpdate){                          #if no filter exists or update is forced, make a filter 
  # --- subset the manual annots to only those with data in the categories used, and optionally by the screen column
  relevantSpecies = manualAnnots[!is.na(manualAnnots[[phenotypeColumn]]),]#remove all species which are not part of the specified categories
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
speciesValues = relevantSpecies[[phenotypeColumn]]                                  #extract the category of each species (in same order)
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
