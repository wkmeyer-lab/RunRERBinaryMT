# -- Usage:
# This script creates a phenotype tree of a phenotype with a categorical value.
# If you wish to only include select species in the analysis, add a new column to your data as a "screen column". Set species to include to have a 1 in that column, and other species to 0. Use the 's' arg to specify the screen column's name.
# Similarly, if you wish to only include species with specific category values, specify used categories with the 'c' argument. 
# Note the various options to combine categories -- see the 'u' and 'o' arguments.
# If you want to perform an analysis with a category composed of several other categories combined, see "RunMergedAnalysis.R". Include them as separate categories in this script. 

# -- Libraries 
source("Src/Reu/loadLibraries.R")
loadLibraries() # Note that this library load includes cluster run detection. 
source("Src/Reu/ZonomNameConvertVectorCommon.R")
#-------------------------------
# --- Command arguments list ---
# - Required -
# r  = filePrefix                                         This is a prefix used to organize and separate files by analysis run. Always required. 
# m  = mainTreeFilename.txt or .rds                       This sets the location of the mainTrees file
# td = spreadSheetFilename.csv                            This sets the spreadsheet to read the data from 
# tp = "phenotypeColumn"                                  This is the column with the phenotype information
# tt = "tipColumn"                                        This sets the column with the tip names as they appear in the mainTrees file.

# - Required, Defaults likely acceptable -
# f  = forceUpdate <T or F>                               This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# tn = "nameColumn"                                       This sets the column with common names of the species. 
# tc = <c("nameOfCategory1","nameOfCategory2")>           This sets a list of category names to be used in the analysis. Species with values not in this list will not be included in the analysis. Include ALL categories to be included (even those merged later).  
# - Optional - 
# h = <T>                                                 This ARG causes the script to return its ARG list, instead of running the script. 
# ts = "screenColumn"                                     This is a column which must have a value of 1 for the species to be included. 
# tr = list(c("replace1", "with1"),c("replace2, with2"))  This combined multiple entries to a new value. The new value can be an existing category, or a new label. 
# tm = list(c("phenotype1", "intophen1"), c("2", "i2"))   This causes combination phenotypes to be merged into the second phenotype, but does not replace standalone phenotypes
# ta = "ancestralTrait"                                   This can be used to set all non-terminal branches to this category. Use be one of the categories in the list. 
# tb = c( "colors", "in ","the", "palette")               This sets the color palette used in the tree graphs. 
# tl = <T or F>                                           This sets if alternate species sets should be created. 
# ta = <number of alternates>                             This sets the number of alternate species sets to generate, if alternate species sets are generated. 
# tz = <T or F>                                           DO NOT USE UNLESS SPECIFIED -- This determines if the liam infrence nodes should be added 
# q = <ER or SYM or ARD>                                  This sets the model type used to estimate ancestral branches.  
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
  
  args = c('r=ComplexDietHighSugarAllSpeciesTest', 'm=data/zoonomiaAllMammalsTrees.rds', 'td=Data/mergedData.csv', 'tp=ComplexHerbDiet', 'f=T', 'q=ER', 'tt=ZoonomiaTip', 'tz=T', 'tl=F',
           'tc=c(
              "Vertivore", "Omnivore", "Invertivore", "Frugivore", "Generalist Herbivore", "Granivore", "Mixed Herbivore", "Nectivore", "Strict Frugivore", "Strict Nectivore"
            )',
           'tr=list(
            c("Mixed Herbivore", "HerbivoreGeneralist"), c("Generalist Herbivore", "HerbivoreGeneralist"), c("Granivore", "HerbivoreGeneralist"), 
            c("Strict Frugivore", "SugarHerbviore"), c("Strict Nectivore", "SugarHerbviore"),c("Frugivore", "SugarHerbviore"), c("Nectivore", "SugarHerbviore")
           )',
           'tb=c("darkgreen", "darkblue","black", "yellow", "red")'
  )
  
    
}
#-------------------------------


# --- Standard start-up code; identical in all scripts---
source("Src/Reu/StandardizedStartupCode.R") #This runs the standardized code block at this location, ensuring identity across scripts. 
#This code reads the "report args" help ARG, and reutrns the args if it is active. 
#It then reads the filePrefix and sets the output directory to that file prefix, and reads the forceUpdate arg. 


# --- Argument Imports ---
{ # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
  #usedArgs = c('r', 'm', 'td', 'tp', 'tt', 
  #             'f', 'tn','tc',
  #             'h', 'ts', 'tr', 'tm', 'ta', 'tb', 'tl', 'ta', 'tz', 'q',
  #             'pm', 'pk', 'ps', 'pd')
  
  # Defaults
  mainTreesLocation = "Results/mainTreesFromAlignments.rds"
  spreadSheetLocation = "Data/mergedData.csv"
  phenotypeColumn = NULL
  tipCol = "tipName"
  nameColumn = NULL
  categoryList = NULL
  screenColumn = NULL
  useScreen = F
  substitutions = NULL
  mergeOnlys = NULL
  ancestralTrait = NULL
  colorSet = c( "darkgreen", "darkblue","black", "red", "yellow")
  generateAlternates = F
  useLiam = F
  modelType = "ER"
  usingPruning = F
  manualPruningSpecies = NULL
  pruningCutoff = NULL
  
  #Load from args 
  source("Src/Reu/ArgumentImport.R")

}
palette(colorSet)

# ---------------------------
# ----- SCRIPT MAIN CODE ----
# ---------------------------


#                   ------- Code Body --------  
#read in the tree based on filetype extension
if(file_ext(mainTreesLocation) == "rds"){
  if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
}else{
  if(!exists("mainTrees")){mainTrees = readTrees(mainTreesLocation)} 
}


manualAnnots = read.csv(spreadSheetLocation)                      #load the manual annotations file holding the phenotype data
manualAnnots[[phenotypeColumn]] = trimws(manualAnnots[[phenotypeColumn]])               #trim away whitespace to allow for better matching 

#Use all categories, if no categories specified. 
if(is.null(categoryList)){
  categoryList = unique(manualAnnots[[phenotypeColumn]])
  message("Using all categories in phenotype column.")
  print(categoryList)
}


# - Merge hyrbid of either substituted phenotypes or merge-only phenotypes - 
if(!is.null(substitutions) & !all(is.na(substitutions))){                                                    #Consider species with multiple combined categories as the merged category
  for( i in 1:length(substitutions)){                                           #Eg if [X] is replaced with [Y], [X/Y] becomes [Y]
    substitutePhenotypes = substitutions[[i]]
    message(paste("Combining", substitutePhenotypes[1], "/", substitutePhenotypes[2]))
    entriesWithPhen1 = grep(substitutePhenotypes[1], manualAnnots[[phenotypeColumn]])
    entriesWithPhen2 = grep(substitutePhenotypes[2], manualAnnots[[phenotypeColumn]])
    combineEntries = which(entriesWithPhen1 %in% entriesWithPhen2)
    combineIndexes = entriesWithPhen1[combineEntries]
    manualAnnots[[phenotypeColumn]][combineIndexes] = substitutePhenotypes[2]
  }
}

if(!is.null(mergeOnlys) & ! all(is.na(mergeOnlys))){                                                    #Consider species with multiple combined categories as the merged category
  for( i in 1:length(mergeOnlys)){                                           #Eg if [X] is replaced with [Y], [X/Y] becomes [Y]
    substitutePhenotypes = mergeOnlys[[i]]
    message(paste("Merging Hybrids of", substitutePhenotypes[1], "/", substitutePhenotypes[2], "to", substitutePhenotypes[2]))
    entriesWithPhen1 = grep(substitutePhenotypes[1], manualAnnots[[phenotypeColumn]])
    entriesWithPhen2 = grep(substitutePhenotypes[2], manualAnnots[[phenotypeColumn]])
    combineEntries = which(entriesWithPhen1 %in% entriesWithPhen2)
    combineIndexes = entriesWithPhen1[combineEntries]
    manualAnnots[[phenotypeColumn]][combineIndexes] = substitutePhenotypes[2]
  }
}


# - Species Filter - 
speciesFilterFilename = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #set a filename for the species filter based on the prefix 

if(!file.exists(speciesFilterFilename) | forceUpdate){                          #if no filter exists or update is forced, make a filter 
  # --- subset the manual annots to only those with data in the categories used, and optionally by the screen column
  relevantSpecies = manualAnnots[manualAnnots[[phenotypeColumn]] %in% categoryList,]#remove all species which are not part of the specified categories
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
      if(all(is.null(allProtectedSpecies))){allProtectedSpecies = ""}
      
      workingTree = mainTrees$masterTree
      workingTree = drop.tip(workingTree, which(!workingTree$tip.label %in% speciesFilter))
      
      fewGeneSpecies = dropFewGeneSpecies(mainTrees, workingTree, nameConversionColumn = tipCol, nameConversionData = spreadSheetLocation)
      if(length(which(fewGeneSpecies %in% allProtectedSpecies)) > 0){
        fewGeneSpecies = fewGeneSpecies[-which(fewGeneSpecies %in% allProtectedSpecies)]
      }
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
    if(!all(is.null(manualPruningSpecies))){
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
  relevantSpecies = manualAnnots[manualAnnots[[phenotypeColumn]] %in% categoryList,]#remove all species which are not part of the specified categories
  if(useScreen){                                                                #if using a screening collumn 
    relevantSpecies = relevantSpecies[ relevantSpecies[[screenColumn]] %in% 1, ]  #remove all species not positive for that collumn 
  }
  relevantSpecies = relevantSpecies[!relevantSpecies[[tipCol]] %in% "", ]   
}

# - Phenotype Vector - 
filteredSpecies = relevantSpecies[relevantSpecies[[tipCol]] %in% speciesFilter, ] 
speciesNames = filteredSpecies[[tipCol]]                                         #Exract the tip name of each species
speciesCategories = filteredSpecies[[phenotypeColumn]]                              #extract the category of each species (in same order)

phenotypeVector = speciesCategories                                             #combine those into⌄
names(phenotypeVector) = speciesNames                                           #the format the functions expect
if(!is.null(substitutions) & !all(is.na(substitutions))){
  for( i in 1:length(substitutions)){
    substitutePhenotypes = substitutions[[i]]
    message(paste("replacing", substitutePhenotypes[1], "with", substitutePhenotypes[2]))
    phenotypeVector = gsub(substitutePhenotypes[1], substitutePhenotypes[2], phenotypeVector)
  }
}

phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
saveRDS(phenotypeVector, file = phenotypeVectorFilename)                        #save the phenotype vector


if(!useLiam){
  # - Make common name versions of objects (used in visualization) - 
  commonMainTrees = mainTrees
  commonMainTrees$masterTree$tip.label = ZonomNameConvertVectorCommon(commonMainTrees$masterTree$tip.label, annotationLocation = spreadSheetLocation, tipColumn = tipCol, commonColumn = nameColumn)
  commonPhenotypeVector = phenotypeVector
  names(commonPhenotypeVector) = ZonomNameConvertVectorCommon(names(commonPhenotypeVector), annotationLocation = spreadSheetLocation, tipColumn = tipCol, commonColumn = nameColumn)
  commonSpeciesFilter = ZonomNameConvertVectorCommon(speciesFilter, annotationLocation = spreadSheetLocation, tipColumn = tipCol, commonColumn = nameColumn)
  
  # - Categorical Tree - 
  treeImageFilename = paste(outputFolderName, filePrefix, "CategoricalTree.pdf", sep="") #make a filename based on the prefix
  pdf(treeImageFilename, height = length(phenotypeVector)/18, width = 10)                     #make a pdf to store the plot, sized based on tree size
    commonCategoricalTree = char2TreeCategorical(commonPhenotypeVector, commonMainTrees, commonSpeciesFilter, model = modelType, anctrait = ancestralTrait, plot = T)
  
    categoricalTree = char2TreeCategorical(phenotypeVector, mainTrees, speciesFilter, model = modelType, anctrait = ancestralTrait, plot = T) #use the phenotype vector to make a tree
  dev.off()                                                                       #save the plot to the pdf
}else{
  
  masterTree = mainTrees$masterTree
  
  nodesToAdd = c(455, 457, 471, 650, 492)
  names(nodesToAdd) = c("Mammalia", "Marsupalia", "Placentalia", "Chiroptera", "Primates")
  #phenToAdd = c("Insectivore", "Insectivore", "Insectivore", "Insectivore", "Omnivore")
  phenToAdd = c("Invertivore", "Invertivore", "Invertivore", "Invertivore", "Omnivore")
  #nodesToAdd = c(785)
  #names(nodesToAdd) = c("Placentalia")
  #phenToAdd = c("Carnivore")
  
  masterTreeAdded = masterTree
  for(i in 1:length(nodesToAdd)){
    M<-matchNodes(masterTree,masterTreeAdded,method="distances")
    masterTreeAdded<-bind.tip(masterTreeAdded,names(nodesToAdd)[i],edge.length=0,
                              where=M[which(M[,1]==as.numeric(nodesToAdd[i])),2])
  }
  mainTrees$masterTree = masterTreeAdded
  
  
  names(phenToAdd) = names(nodesToAdd)
  
  phenotypeVector = append(phenotypeVector, phenToAdd)
  speciesFilter = append(speciesFilter, names(phenToAdd))
  
  
  # - Make common name versions of objects (used in visualization) - 
  commonMainTrees = mainTrees
  commonMainTrees$masterTree$tip.label = ZonomNameConvertVectorCommon(commonMainTrees$masterTree$tip.label, annotationLocation = spreadSheetLocation, tipCol = tipCol, commonColumn = nameColumn)
  commonPhenotypeVector = phenotypeVector
  names(commonPhenotypeVector) = ZonomNameConvertVectorCommon(names(commonPhenotypeVector), annotationLocation = spreadSheetLocation, tipColumn = tipCol, commonColumn = nameColumn)
  commonSpeciesFilter = ZonomNameConvertVectorCommon(speciesFilter, annotationLocation = spreadSheetLocation, tipColumn = tipCol, commonColumn = nameColumn)
  
  # - Categorical Tree - 
  treeImageFilename = paste(outputFolderName, filePrefix, "CategoricalTree.pdf", sep="") #make a filename based on the prefix

  pdf(treeImageFilename, height = length(phenotypeVector)/18, width = 10)                     #make a pdf to store the plot, sized based on tree size
  commonCategoricalTree = char2TreeCategorical(commonPhenotypeVector, commonMainTrees, commonSpeciesFilter, model = modelType, anctrait = ancestralTrait, plot = F)
  categoricalTree = char2TreeCategorical(phenotypeVector, mainTrees, speciesFilter, model = modelType, anctrait = ancestralTrait, plot = F) #use the phenotype vector to make a tree
  
  commonCategoricalTreeExtraTip = commonCategoricalTree
  categoricalTreeExtraTip = categoricalTree
  
  commonCategoricalTree = drop.tip(commonCategoricalTree, names(nodesToAdd))
  categoricalTree = drop.tip(categoricalTree, names(nodesToAdd))
  mainTrees$masterTree = drop.tip(mainTrees$masterTree, names(nodesToAdd))
  commonMasterAdded = commonMainTrees$masterTree
  commonMainTrees$masterTree = drop.tip(commonMainTrees$masterTree, names(nodesToAdd))
  
  #plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMainTrees$masterTree)
  #plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)
  
  #plotTreeCategorical(commonCategoricalTreeExtraTip, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMasterAdded)
  #plotTreeCategorical(categoricalTreeExtraTip, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = masterTreeAdded)
  
  categoryLegend = names(table(phenotypeVector))
  
  plotTreeCategorical(commonCategoricalTree, categoryLegend, master = commonMainTrees$masterTree)
  plotTreeCategorical(categoricalTree, categoryLegend, master = mainTrees$masterTree)
  
  plotTreeCategorical(commonCategoricalTreeExtraTip, categoryLegend, master = commonMasterAdded)
  plotTreeCategorical(categoricalTreeExtraTip, categoryLegend, master = masterTreeAdded)
  dev.off()  
  rm(mainTrees) #removes the altered mainTrees object, so that it doesn't affect any other scripts. This is a potential issue unique to the mainTrees object, because it is not reloaded if detected, due to long load times.)
}


scientificCategoricalTree = categoricalTree
scientificCategoricalTree$tip.label = ZonomNameConvertVectorCommon(commonCategoricalTree$tip.label, annotationLocation = spreadSheetLocation, tipColumn = tipCol, toScientific = T)

categoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
saveRDS(categoricalTree, categoricalTreeFilename)                               #save the tree
categoricalCommonTreeFilename = paste(outputFolderName, filePrefix, "CategoricalCommonTree.rds", sep="") #make a filename based on the prefix
saveRDS(commonCategoricalTree, categoricalCommonTreeFilename)
scientificCategoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalScientificTree.rds", sep="") #make a filename based on the prefix
saveRDS(scientificCategoricalTree, scientificCategoricalTreeFilename)

phenotypeVectorSaving = data.frame(names(phenotypeVector), names(commonPhenotypeVector), phenotypeVector)
directReadablePhenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.csv",sep="")
write.csv(phenotypeVectorSaving, directReadablePhenotypeVectorFilename)



if(generateAlternates){
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
}

