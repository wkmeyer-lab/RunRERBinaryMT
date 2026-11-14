# -- Usage: 
# This script is a continuation of the MakeCategoricalPhenotypeTree.R Script, but partitioned. 
# It is done this way so that this portion of the script can call the earlier section to generate an all-species tree. 

useCladeLimit = F
useStartDrop = T

#Load in necessary script
source("Src/Reu/CategoricalDropTip.R")
mainTrees = readRDS(mainTreesLocation) #refreshes the main tree from any changes made during pruning step 

#Load in the full tree with all species
fullTreePrefix = paste0(filePrefix, "AllSpecies")
fullTreeOutputFolder = paste0(outputFolderNameNoSlash, "AllSpecies", "/")





#Used to convert v2 back to v1 to reuse all species
#fullTrePrefixStub = substr(filePrefix, 1, (nchar(filePrefix)-1))
#fullTreeOutputFolderStub = substr(outputFolderNameNoSlash, 1, (nchar(outputFolderNameNoSlash)-1))
#fullTreePrefix = paste0(fullTrePrefixStub, "AllSpecies")
#fullTreeOutputFolder = paste0(fullTreeOutputFolderStub, "AllSpecies", "/")


tryCatch({
  allSpeciesTreeFilename = paste(fullTreeOutputFolder, fullTreePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
  allSpeciesTree = readRDS(allSpeciesTreeFilename)
  haveAllSpeciesTree = T
}, 
error=function(i){
  haveAllSpeciesTree = F
  message("No All Species tree found. Generating alternates, but not creating trees for them. ")
}
)


#Figure out the number of short terminal branches in the main analysis 
masterTree = mainTrees$masterTree
mainAnalysisTree = mainTrees$masterTree
tipsDroppedForMain = mainAnalysisTree$tip.label[which(!mainAnalysisTree$tip.label %in% categoricalTree$tip.label)]
mainAnalysisTree = drop.tip(mainAnalysisTree, tipsDroppedForMain)
terminalEdges = which(mainAnalysisTree$edge[,2] <= length(mainAnalysisTree$tip.label))
mainNumShortBranches = length(which(mainAnalysisTree$edge.length[terminalEdges] <= pruningCutoff))


masterTerminalEdges = which(masterTree$edge[,2] <= length(masterTree$tip.label))
masterShortBranches = which(masterTree$edge.length <= pruningCutoff)

masterShortTerminals = which((masterTree$edge[,2] <= length(masterTree$tip.label)) & (masterTree$edge.length <= pruningCutoff))
masterShortTips = masterTree$edge[masterShortTerminals,2]



phenotypeSizes = table(phenotypeVector)

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

mastertreeDataPhenotype = fullDataPhenotype[names(fullDataPhenotype) %in% mainTrees$masterTree$tip.label]
length(mastertreeDataPhenotype)
masterTreePhenotypeSizes = table(mastertreeDataPhenotype)





speciesSets = list()
for(j in 1:length(phenotypeSizes)){
  speciesSets[[j]] = names(mastertreeDataPhenotype)[which(mastertreeDataPhenotype == names(phenotypeSizes)[j])]
}
numTips = sum(phenotypeSizes) # This is the number of tips, calculated here so it is done once instead of in the loop


stableSpeciesSets = speciesSets
stablePhenotypeSizes = phenotypeSizes




alternateSets = list()
i=1
minShort = length(mainAnalysisTree$edge.length) 
while(length(alternateSets) < 100){
  
  randomizedSpeciesSet = character()
  if(useCladeLimit){
    #reset the selection data to undo previous cuts
    phenotypeSizes = stablePhenotypeSizes
    speciesSets = stableSpeciesSets
    
    {
      # add preselected
      numPerClade = 3
      
      Bovidae = manualAnnots[manualAnnots$MSWC_Family %in% "Bovidae",][[tipCol]]
      Pteropodidae = manualAnnots[manualAnnots$MSWC_Family %in% "Pteropodidae",][[tipCol]]
      Cervidae = manualAnnots[manualAnnots$MSWC_Family %in% "Cervidae",][[tipCol]]
      Cricetidae = manualAnnots[manualAnnots$MSWC_Family %in% "Cricetidae",][[tipCol]]
      HystricognathiFamilies =  c("Caviidae", "Chinchillidae", "Ctenomyidae", "Dasyproctidae", "Dinomyidae", "Caviidae", "Bathyergidae", "Caviidae", "Hystricidae", "Myocastoridae", "Octodontidae", "Petromuridae", "Thryonomyidae", "Erethizontidae")
      Hystricognathi = manualAnnots[manualAnnots$MSWC_Family %in% HystricognathiFamilies,][[tipCol]]
      
      Vespertilionidae = manualAnnots[manualAnnots$MSWC_Family %in% "Vespertilionidae",][[tipCol]]
      speciesInLargeFamilies = c(Bovidae, Pteropodidae, Cervidae, Cricetidae, Hystricognathi, Vespertilionidae)
      
      selectedBovids = sample(Bovidae, numPerClade)
      selectedPteropodidae = sample(Pteropodidae, numPerClade)
      selectedCervidae = sample(Cervidae, numPerClade)
      selectedCricetidae = sample(Cricetidae, numPerClade)
      selectedHystricognathi = sample(Hystricognathi, numPerClade)
      
      selectedVespertilionidae = sample(Vespertilionidae, numPerClade)
      
      selectedHerbs = c(selectedBovids, selectedPteropodidae, selectedCervidae, selectedCricetidae, selectedHystricognathi)
      selectedInvert = c(selectedVespertilionidae)
      
      
      
      
      
      preselectedSpecies = c(selectedHerbs, selectedInvert)
      speciesToNotAutoselect = speciesInLargeFamilies
      
      
    }
    
    #Handle the preselcted species 
    {
      #Add to sample
      randomizedSpeciesSet = append(randomizedSpeciesSet, preselectedSpecies)
      
      #reduce category sizes 
      for(j in 1:length(speciesSets)){
        numPreselectedInCategory = length(which(preselectedSpecies %in% speciesSets[[j]]))
        phenotypeSizes[j] = phenotypeSizes[j] - numPreselectedInCategory
      }
      #remove from selectable list 
      for(j in 1:length(speciesSets)){
        speciesSets[[j]] = speciesSets[[j]][-!(speciesSets[[j]] %in% speciesToNotAutoselect)]
      }
    }
  }
  
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


if(!haveAllSpeciesTree){
  for(i in 1:length(alternateSets)){
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
  }
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
    
    if(useCladeLimit){
      
      treeImageFilename = paste(outputFolderName, alternateFilePrefix, filePrefix,"CategoricalTree.pdf", sep="") #make a filename based on the prefix
      pdf(treeImageFilename, height = length(currentPhenotypeVector)/18, width = 10)                     #make a pdf to store the plot, sized based on tree size
      palette(c( "darkgreen", "darkblue","black", "red"))
      
      plotTreeCategorical(currentCommonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMainTrees$masterTree)
      plotTreeCategorical(currentCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)
      
      dev.off()  
    }
    
    
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
  
  print(table(currentPathsObject))
  
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