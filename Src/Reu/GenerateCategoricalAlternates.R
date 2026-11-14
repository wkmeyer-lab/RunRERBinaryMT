# -- Usage: 
# This script is a continuation of the MakeCategoricalPhenotypeTree.R Script, but partitioned. 
# It is done this way so that this portion of the script can call the earlier section to generate an all-species tree. 


useStartDrop = T
if(is.null(pruningCutoff)){ # This is a catch to allow alterate generation even without pruning (only useful for debugging)
  pruningCutoff = 0 
}

mainTrees = readRDS(mainTreesLocation) #refreshes the main tree from any changes made during pruning step 

#Load in the full tree with all species
fullTreePrefix = paste0(filePrefix, "AllSpecies")
fullTreeOutputFolder = paste0(outputFolderName, "AllSpecies", "/")


allSpeciesTreeFilename = paste(fullTreeOutputFolder, fullTreePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
if(file.exists(allSpeciesTreeFilename)){
  allSpeciesTree = readRDS(allSpeciesTreeFilename)
  message("All Species tree found. Generating Alternates. ")
  haveAllSpeciesTree = T
}else{
  message("No All Species tree found. Generating a tree with all species to subset from. ")
  haveAllSpeciesTree = T
  
  mainOutputFolder = outputFolderName
  mainFilePrefix = filePrefix
  
  
  outputFolderName = fullTreeOutputFolder
  if(!dir.exists(outputFolderName)){                       
    dir.create(outputFolderName)
  }
  filePrefix = fullTreePrefix
  generateAlternates = F
  usingPruning = F
  
  alreadyBlockingArgImports = exists("NoArgUpdate")
  if(!alreadyBlockingArgImports){NoArgUpdate = "Exists"}
  CurrentScript = paste("GenerateCategoricalAlternates.R -> MakeCategoricalPhenotypeTree.R")
  source("Src/Loc/MakeCategoricalPhenotypeTree.R")
  if(!alreadyBlockingArgImports){rm(NoArgUpdate)}
  
  
  CurrentScript = "GenerateCategoricalAlternates.R"
  outputFolderName = mainOutputFolder
  filePrefix = mainFilePrefix
  generateAlternates = T
  message("Created AllSpecies Tree")
}


#Figure out the number of short terminal branches in the main analysis 
masterTree = mainTrees$masterTree
mainAnalysisTree = mainTrees$masterTree
tipsDroppedForMain = mainAnalysisTree$tip.label[which(!mainAnalysisTree$tip.label %in% categoricalTree$tip.label)]
mainAnalysisTree = drop.tip(mainAnalysisTree, tipsDroppedForMain)
terminalEdges = which(mainAnalysisTree$edge[,2] <= length(mainAnalysisTree$tip.label))
mainNumShortBranches = length(which(mainAnalysisTree$edge.length[terminalEdges] <= pruningCutoff))

message(paste(mainNumShortBranches, "Branches below cutoff in primary analysis"))

masterTerminalEdges = which(masterTree$edge[,2] <= length(masterTree$tip.label))
masterShortBranches = which(masterTree$edge.length <= pruningCutoff)

masterShortTerminals = which((masterTree$edge[,2] <= length(masterTree$tip.label)) & (masterTree$edge.length <= pruningCutoff))
masterShortTips = masterTree$edge[masterShortTerminals,2]


phenotypeSizes = table(phenotypeVector)
print("Category sizes in main analysis ")
print(phenotypeSizes)

fullSpeciesNames = relevantSpecies[[tipCol]]                                         #Exract the tip name of each species
fullSpeciesCategories = relevantSpecies[[phenotypeColumn]]                              #extract the category of each species (in same order)
fullDataPhenotype = fullSpeciesCategories
names(fullDataPhenotype) = fullSpeciesNames 
if(!is.null(substitutions) & !all(is.na(substitutions))){
  for( i in 1:length(substitutions)){
    substitutePhenotypes = substitutions[[i]]
    message(paste("replacing", substitutePhenotypes[1], "with", substitutePhenotypes[2]))
    fullDataPhenotype = gsub(substitutePhenotypes[1], substitutePhenotypes[2], fullDataPhenotype)
  }
}
length(fullDataPhenotype)

fullPhenotypeSizes = table(fullDataPhenotype)
print("Category sizes in Phenotype Spreadsheet")
print(fullPhenotypeSizes)


mastertreeDataPhenotype = fullDataPhenotype[names(fullDataPhenotype) %in% mainTrees$masterTree$tip.label]
length(mastertreeDataPhenotype)
masterTreePhenotypeSizes = table(mastertreeDataPhenotype)

print("Category sizes in Master Tree")
print(masterTreePhenotypeSizes)




speciesSets = list()
for(j in 1:length(phenotypeSizes)){
  speciesSets[[j]] = names(mastertreeDataPhenotype)[which(mastertreeDataPhenotype == names(phenotypeSizes)[j])]
}
numTips = sum(phenotypeSizes) # This is the number of tips, calculated here so it is done once instead of in the loop


stableSpeciesSets = speciesSets
stablePhenotypeSizes = phenotypeSizes



message("Beginning alternate search")
alternateSets = list()
i=1
minShort = length(mainAnalysisTree$edge.length) 

#Add in a check for alternate generation without pruning 
if(sum(phenotypeSizes) == sum(masterTreePhenotypeSizes)){
  message("------ WARNING!!! ------")
  message("It seems that you are generating alterates without pruning the tree.")
  message("Because all species are included in the main analysis, any alternative selection of species sets will produce the same thing -- all of the species.")
  message("This makes alternates functionally pointless when not pruning.")
  message("The ability to generate alternates when not pruning is only present for demonstration and debugging.")
  message("If this is an actual analysis, you should probably either turn this off, or prune your tree.")
  message("------")
  for(i in 1:numberOfAlternates)
  alternateSets[[i]] = names(phenotypeVector) 
}

while(length(alternateSets) < numberOfAlternates){
  
  randomizedSpeciesSet = character()
  
  if(useStartDrop){
    #message(i)
    keptFraction = 2
    validRemoval = F
    regectedRemovals = 0
    
    while(validRemoval == F){
      speciesSets = stableSpeciesSets
      startRemoval = sample(masterShortTips, length(masterShortTips)/keptFraction)
      startRemoveTips = masterTree$tip.label[startRemoval]
      validSize = logical()
      
      #remove from selectable list 
      for(j in 1:length(speciesSets)){
        speciesSets[[j]] = speciesSets[[j]][!(speciesSets[[j]] %in% startRemoveTips)]
        #length(speciesSets[[j]])
        if(length(speciesSets[[j]]) > phenotypeSizes[j]){
          validSize[j] = T
        }else{
          validSize[j] = F
          #message(names(phenotypeSizes[j]))
          #message(paste(length(speciesSets[[j]]), "<", phenotypeSizes[j]))
        }
      }
      removePhens = fullDataPhenotype[names(fullDataPhenotype) %in% startRemoveTips] 
      #print(masterTreePhenotypeSizes)
      #print(phenotypeSizes)
      #print(table(removePhens))
      #length(speciesSets[[j]])
      if(all(validSize)){
        validRemoval = T
        #message(paste("Removal Accepted after rejections:", regectedRemovals))
        regectedRemovals = 0
      }else{
        regectedRemovals = regectedRemovals +1
        #message("Removal Rejected")
      }
      
    }
    
    
  }
  
  
  for(j in 1:length(phenotypeSizes)){
    #numberOfSpecies = phenotypeSizes[j]
    #speciesSet = mastertreeDataPhenotype[which(mastertreeDataPhenotype == names(phenotypeSizes)[j])]
    chosenSpecies = sample(speciesSets[[j]], phenotypeSizes[j])
    randomizedSpeciesSet = append(randomizedSpeciesSet, chosenSpecies)
  } 
  
  
  
  
  
  
  tipsToDrop = mainTrees$masterTree$tip.label[!mainTrees$masterTree$tip.label %in% randomizedSpeciesSet]
  testTree = drop.tip(mainTrees$masterTree, tipsToDrop)
  i=i+1
  if(i %% 10000 == 0){message(i)}
  terminalEdges = which(testTree$edge[,2] <= numTips)
  
  numberOfShortBranches = length(which(testTree$edge.length[terminalEdges] <= pruningCutoff))
  if(numberOfShortBranches < minShort){
    minShort = numberOfShortBranches
    message(i)
    message(numberOfShortBranches)
  }
  
  
  if(numberOfShortBranches <= mainNumShortBranches*1.5){
    message("Found Valid Alternate")
    message(i)
    alternateTips = testTree$tip.label
    alternateTips = list(alternateTips)
    alternateSets = append(alternateSets, alternateTips)
  }
}

saveRDS(alternateSets, paste0(outputFolderName, filePrefix, "AlternatePruningSpecies.rds"))


if(!dir.exists(paste0("Output/", filePrefix, "/Alternates"))){                                      #Make output directory if it does not exist
  dir.create(paste0("Output/", filePrefix, "/Alternates"))
}


if(haveAllSpeciesTree){
  
  
  commonMainTrees = mainTrees
  commonMainTrees$masterTree$tip.label = ZonomNameConvertVectorCommon(commonMainTrees$masterTree$tip.label, annotationLocation = spreadSheetLocation, tipColumn = tipCol, commonColumn = nameColumn)
  
  
  for(i in 1:length(alternateSets)){
    message(i)
    currentSet = alternateSets[[i]]
    alternateFilePrefix = paste0("/Alternates/Alternate", i)
    
    
    currentPhenotypeVector = fullDataPhenotype[which(names(fullDataPhenotype) %in% currentSet)]
    currentCommonPhenotypeVector = currentPhenotypeVector
    names(currentCommonPhenotypeVector) = ZonomNameConvertVectorCommon(names(currentCommonPhenotypeVector), annotationLocation = spreadSheetLocation, tipColumn = tipCol, commonColumn = nameColumn)
    
    
    currentPhenotypeVectorSaving = data.frame(names(currentPhenotypeVector), names(currentCommonPhenotypeVector), currentPhenotypeVector)
    currentDirectReadablePhenotypeVectorFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalPhenotypeVector.csv",sep="")
    #write.csv(currentPhenotypeVectorSaving, currentDirectReadablePhenotypeVectorFilename)  
    
    currentPhenotypeVectorFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
    saveRDS(currentPhenotypeVector, file = currentPhenotypeVectorFilename)                        #save the phenotype vector
    
    currentSpeciesFilterFilename =  paste(outputFolderName, alternateFilePrefix, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 
    currentSpeciesFilter = currentSet
    saveRDS(currentSpeciesFilter, currentSpeciesFilterFilename)
    
    
    tipsToRemove = allSpeciesTree$tip.label[!allSpeciesTree$tip.label %in% currentSet]
    currentTree = categoricalDropTip(allSpeciesTree, tipsToRemove)
    currentCommonSpeciesFilter = ZonomNameConvertVectorCommon(currentSet, annotationLocation = spreadSheetLocation, tipColumn = tipCol, commonColumn = nameColumn)
    
    # - Categorical Tree - 
    currentCategoricalTree = currentTree
    currentCommonCategoricalTree = currentTree
    currentCommonCategoricalTree$tip.label = ZonomNameConvertVectorCommon(currentCommonCategoricalTree$tip.label, annotationLocation = spreadSheetLocation, tipColumn = tipCol, commonColumn = nameColumn)
    
    
    currentCategoricalTreeFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
    saveRDS(currentCategoricalTree, currentCategoricalTreeFilename)                               #save the tree
    currentCategoricalCommonTreeFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalCommonTree.rds", sep="") #make a filename based on the prefix
    saveRDS(currentCommonCategoricalTree, currentCategoricalCommonTreeFilename)
    
  }
}


#Try to test for broken paths 
alternateSets = readRDS(paste0(outputFolderName, filePrefix, "AlternatePruningSpecies.rds"))
brokenPaths = integer()
for(i in 1:length(alternateSets)){
  message(i)
  currentSet = alternateSets[[i]]
  alternateFilePrefix = paste0("/Alternates/Alternate", i)
  currentCategoricalTreeFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
  currentCommonCategoricalTreeFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalCommonTree.rds", sep="") #make a filename based on the prefix
  currentPhenotypeVectorFilename = paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
  currentSpeciesFilterFilename =  paste(outputFolderName, alternateFilePrefix, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 
  
  currentPhenotypeTree = readRDS(currentCategoricalTreeFilename)
  currentSpeciesFilter = readRDS(currentSpeciesFilterFilename)
  
  currentPathsObject = tree2Paths(currentPhenotypeTree, mainTrees, useSpecies = currentSpeciesFilter, categorical = TRUE)
  currentPathsFilename =  paste(outputFolderName, alternateFilePrefix, filePrefix, "CategoricalPathsFile.rds",sep="") #set a filename for the species filter based on the prefix 
  saveRDS(currentPathsObject, currentPathsFilename)
  
  #print(table(currentPathsObject))
  
  if(length(table(currentPathsObject)) != length(phenotypeSizes)){
    message("Alternate is broken")
    file.rename(currentCategoricalTreeFilename, paste(outputFolderName, "Alternates/Broken/Alternate", i, filePrefix, "CategoricalTree.rds", sep=""))
    file.rename(currentCommonCategoricalTreeFilename, paste(outputFolderName, "Alternates/Broken/Alternate", i, filePrefix, "CategoricalCommonTree.rds", sep=""))
    file.rename(currentPhenotypeVectorFilename, paste(outputFolderName, "Alternates/Broken/Alternate", i, filePrefix, "CategoricalPhenotypeVector.rds", sep=""))
    file.rename(currentSpeciesFilterFilename, paste(outputFolderName, "Alternates/Broken/Alternate", i, filePrefix, "SpeciesFilter.rds", sep=""))
    
    brokenPaths = append(brokenPaths, i)
    
    
  }else{
    message("Paths fine, maintaining")
  }
  
}
if(length(brokenPaths) > 0){
  if(!dir.exists(paste0(outputFolderName, "/Alternates/Broken"))){                                      #Make output directory if it does not exist
    dir.create(paste0(outputFolderName, "/Alternates/Broken"))
  }
  
  brokenAlternateSet = alternateSets[brokenPaths]
  names(brokenAlternateSet) = brokenPaths
  brokenAlternatesFilename = paste(outputFolderName, "Alternates/Broken/", filePrefix, "BrokenAlternateSets.rds", sep="")
  saveRDS(brokenAlternateSet, brokenAlternatesFilename)
  
  backupAlternatesFilename = paste(outputFolderName, "Alternates/Broken/", filePrefix, "OriginalAlternateSets.rds", sep="")
  saveRDS(alternateSets, backupAlternatesFilename)
}



#Exit the script if being run automatically before making any changes 
if(!length(brokenPaths) ==0){
  stop("SOME OF THE ALTERNATES ARE PRODUCING BROKEN PATHS! Check the Broken folder to investigate further.")
  alternateSets = alternateSets[-brokenPaths]
}