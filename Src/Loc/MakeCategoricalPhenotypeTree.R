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
source("Src/Reu/CategoricalDropTip.R")
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
# tn = "ancestralTrait"                                   This can be used to set all non-terminal branches to this category. Use be one of the categories in the list. 
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
if(!exists("NoArgUpdate")){ #This line allows other scripts to run this code without regenerating args, if that script creates a "NoArgUpdate" variable.
    
  source("Src/Reu/StandardizedStartupCode.R") #This runs the standardized code block at this location, ensuring identity across scripts. 
  #This code reads the "report args" help ARG, and reutrns the args if it is active. 
  #It then reads the filePrefix and sets the output directory to that file prefix, and reads the forceUpdate arg. 
  
  
  # --- Argument Imports ---
  { # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
    #usedArgs = c('r', 'm', 'td', 'tp', 'tt', 
    #             'f', 'tn','tc',
    #             'h', 'ts', 'tr', 'tm', 'ta', 'tb', 'tl', 'ta', 'tz', 'q',
    #             'pm', 'pk', 'ps', 'pd')
    
    CurrentScript = "MakeCategoricalPhenotypeTree.R"
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
}
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

  
  source("Src/Reu/GenerateCategoricalAlternates.R")
  
}

