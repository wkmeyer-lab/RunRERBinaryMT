#---------------------------------------------------------------------
# --- RunAlternateEnrichments.R  --- 
# --------------------------------------------------------------------
# Defaults
gmtFileLocation = "Data/enrichmentGmtFile.gmt"
usePermulations = FALSE
usePermulationPValOverride = FALSE
permulationPValOverride = NULL 
useCorrelationOverride = FALSE
correlationOverride = NULL
useCategoricalPermulations = FALSE
useSubdirectory = FALSE
subdirectoryValueList = NULL
validGeneArguments = c("NULL", "gene", "Gene")
alternatesSpecified = F
usedAlternates = NULL
validLocalRunAllAlternates = c("local", "sure", "all")



#Gmt file location
if(!any(is.na(cmdArgImport('g', T)))){
  gmtFileLocation = cmdArgImport('g')
  
  gmtFileLocation = cmdArgImport('g')
  
  #remove any arguments which tell later scripts to run a gene analysis 
  gmtFileLocation = gmtFileLocation[-which(gmtFileLocation %in% validGeneArguments)]
  
}else{
  paste("No gmt location arugment, using Data/enrichmentGmtFile.gmt")                          #Report using default
  message("No gmt location arugment, using Data/enrichmentGmtFile.gmt")
}

#Permulation use
if(!is.na(cmdArgImport('p', T))){
  usePermulations = cmdArgImport('p')
  if(usePermulations %in% c("C", "A")){
    usePermulations = TRUE
    useCategoricalPermulations = TRUE
  }
  usePermulations = as.logical(usePermulations)
  if(is.na(usePermulations)){
    usePermulations = TRUE
    message("Use Permualtions value not interpretable as logical. Did you remember to capitalize? Using TRUE.")
  }
}else{
  message("Use Permulations value not specified, using TRUE.")
}

#Import subdirectory
if(!any(is.na(cmdArgImport('s', T)))){
  useSubdirectory = TRUE
  subdirectoryValueList = cmdArgImport('s')
  
}else{
  message("No subdirectory specified.")
  alternatesFolderName = paste0(outputFolderName,"Alternates")
  folderNames = basename(list.dirs(alternatesFolderName, full.names = TRUE, recursive = FALSE))
  if(length(folderNames) >0 ){
    useSubdirectory = TRUE
    subdirectoryValueList = folderNames
    message("Folders not specified, using all folders.")
  }else{
    message("Using main directory")
  }
}

#alternate set
if(!is.na(cmdArgImport('i', T))){
  usedAlternates = cmdArgImport('i')
  if(tolower(usedAlternates) %in% tolower(validLocalRunAllAlternates)){
    message("Run all alternates on local machine CONFIRMED.")
    message("Caution: May have long run times.")
    usedAlternates = NULL
    alternatesSpecified = F
  }else{
    usedAlternates = as.integer(usedAlternates)
    alternatesSpecified = T
  }
}else{
  message("Instance not specified, running on all alternates.")
  if(!clusterRun){
    stop("Attempting to process all alternates on a local machine! This may have run times in excess of 24 hours. Aborting. \n If you wish to avoid this abort, use i=all")
  }
}

#Import permulation p-value Override  
if(usePermulations){
  if(!is.na(cmdArgImport('z', T))){
    usePermulationPValOverride = TRUE
    permulationPValOverride = cmdArgImport('z')
    message("Using Manually specified permulation p-value file.")
  }else{
    message("No Permulation pValue override, using standard CombinedPrunedFastAllPermulationsPValue.rds.")
  }
}

#Import correlation file override
if(!is.na(cmdArgImport('y', T))){
  useCorrelationOverride = TRUE
  correlationOverride = cmdArgImport('y')
  message("Using Manually specified correlation file.")
}else{
  message("No corelation override specified.")
}



#---------------------------------------------------------------------
# --- ProcessAlternates.R  --- 
# --------------------------------------------------------------------

# Defaults
mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"  #this is the standard location on Sol 
binaryPhenotypeTreeLocation = NULL
speciesFilter = NULL
phenotypeStyle = "Categorical"
continousMetric = "diff"
validMetrics = c("diff", "mean", "last")
minSpValue = 10
alternatesSpecified = F
usedAlternates = NULL
validLocalRunAllAlternates = c("local", "sure", "all")

toMerge = NULL
mergeValue = NULL
runMerge= F


#MainTrees Location
if(!is.na(cmdArgImport('m', T))){
  mainTreesLocation = cmdArgImport('m')
}else{
  message("No maintrees arg, using default")
}

#Phenotype Style
if(!is.na(cmdArgImport('a', T))){
  phenotypeStyle = cmdArgImport('a')
  #Convert the various input options; so that the style can be specified with several names 
  if(phenotypeStyle == "b" | phenotypeStyle == "B" | phenotypeStyle == "binary"){phenotypeStyle = "Binary"}
  if(phenotypeStyle == "c" | phenotypeStyle == "C" | phenotypeStyle == "continuous"){phenotypeStyle = "Continuous"}
  if(phenotypeStyle == "g" | phenotypeStyle == "G" | phenotypeStyle == "categorical"){phenotypeStyle = "Categorical"}
}else{
  message("PhenotypeStyle not specified, using continuous")
}

#Continuous Metric 
if(!is.na(cmdArgImport('k', T))){
  continousMetric = cmdArgImport('k')
  if(phenotypeStyle == "Continuous"){ 
    if(!any(continousMetric == validMetrics)){
      stop("Contious metric is not a valid option. Use 'diff', 'mean', or 'last'")
    }
  }
}else{                                                                        #If a continuous phenotype, report using Diff
  if(phenotypeStyle == "Continuous"){
    message("No continuous metric specified, using diff")
  }
}

#min.sp Value
if(!is.na(cmdArgImport('l', T))){
  minSpValue = cmdArgImport('l')
}else{
  message("min.sp not specified, using 10")
}

#alternate set
if(!is.na(cmdArgImport('i', T))){
  usedAlternates = cmdArgImport('i')
  if(tolower(usedAlternates) %in% tolower(validLocalRunAllAlternates)){
    message("Run all alternates on local machine CONFIRMED.")
    message("Caution: May have long run times.")
    usedAlternates = NULL
    alternatesSpecified = F
  }else{
    usedAlternates = as.integer(usedAlternates)
    alternatesSpecified = T
  }
}else{
  message("Instance not specified, running on all alternates.")
  if(!clusterRun){
    stop("Attempting to process all alternates on a local machine! This may have run times in excess of 24 hours. Aborting. \n If you wish to avoid this abort, use i=all")
  }
}

#toMerge
if(!all(is.na(cmdArgImport('e', T)))){
  toMerge = cmdArgImport('e')
  runMerge = T
}else{
  message("No phenotypes set to merge. No merge will occur; if you are not running a merged phenotype, this is fine.")
  runMerge = F
}

#mergeValue
if(!is.na(cmdArgImport('t', T))){
  mergeValue = cmdArgImport('t')
}else{
  if(!runMerge){
    message("No phenotypes set to merge. Script will effectively do nothing; if you are not running a merged phenotype, this is fine.")  
  }else{
    stop("ERROR: Attempting to run merge, but the merged phenotype is unnamed. Exiting.")
  }
  
}

#---------------------------------------------------------------------
# --- CategoricalCombinePairwiseResults.R  --- 
# --------------------------------------------------------------------
# Defaults
significanceCutoff = 0.05
pairwiseSets = NULL
geneSet = NULL
usingGo = F
usingGene = T
validGeneArguments = c("NULL", "gene", "Gene")
writeData = T
usePermulations = F
addPermulations = F

#geneset
if(!all(is.na(cmdArgImport('g', T)))){
  geneSets = cmdArgImport('g', F)
  if(is.null(geneSets) | any(geneSets %in% validGeneArguments)){
    usingGene = T
    message("Null geneset specified or gene specified, performing gene correlation.")
  }else{
    usingGene = F
  }
  
  if(any(!geneSet %in% validGeneArguments)){
    usingGo = T
    geneSet = trimGmtFilenameToSetname(geneSet)
  }else{
    usingGo = F
  }
  
}else{
  message("No geneset specified or gene specified, performing gene correlation.")
  usingGo = F
}

#pairwise Sets
if(!all(is.na(cmdArgImport('s', T)))){
  if(all(cmdArgImport('s', F) == "NULL")){
    folderNames = basename(list.dirs(outputFolderName, full.names = TRUE, recursive = FALSE))
    pairwiseSets = folderNames[grep("-", folderNames)]
  }else{
    pairwiseSets = cmdArgImport('s', F) 
  }
}else{
  folderNames = basename(list.dirs(outputFolderName, full.names = TRUE, recursive = FALSE))
  pairwiseSets = folderNames[grep("-", folderNames)]
  message("Folders not specified, using all folders. Note that if driver analysis is active, this may cause issues.")
}

#Use Permulations 
if(!is.na(cmdArgImport('p', T))){
  if(cmdArgImport('p', F) == "A"){
    usePermulations = T
    addPermulations = T
  }else if(cmdArgImport('p', F) == "C"){
    usePermulations = T
  }else{
    usePermulations = as.logical(cmdArgImport('p', F))
  }
  
}else{
  message("use Permulations not specified, using FLASE")
}

#writeData 
if(!is.na(cmdArgImport('w', T))){
  writeData = as.logical(cmdArgImport('w', F))
}else{
  message("writeData not specified, using TRUE")
} 

#cutoff 
if(!is.na(cmdArgImport('c', T))){
  significanceCutoff = as.numeric(cmdArgImport('c', F))
}else{
  message("pvalue cuttoff not specified, using 0.05")
}  


#---------------------------------------------------------------------
# --- RunEnrichmentAnalysis.R  --- 
# --------------------------------------------------------------------

# Defaults
gmtFileLocation = "Data/enrichmentGmtFile.gmt"
usePermulations = FALSE
usePermulationPValOverride = FALSE
permulationPValOverride = NULL 
useCorrelationOverride = FALSE
correlationOverride = NULL
useCategoricalPermulations = FALSE
useSubdirectory = FALSE
subdirectoryValueList = NULL
validGeneArguments = c("NULL", "gene", "Gene")


#Gmt file location
if(!any(is.na(cmdArgImport('g', T)))){
  gmtFileLocation = cmdArgImport('g', F)
  
  #remove any arguments which tell later scripts to run a gene analysis 
  gmtFileLocation = gmtFileLocation[-which(gmtFileLocation %in% validGeneArguments)]
  
}else{
  message("No gmt location arugment, using Data/enrichmentGmtFile.gmt")
}

#Permulation use
if(!is.na(cmdArgImport('p', T))){
  usePermulations = cmdArgImport('p', F)
  if(usePermulations %in% c("C","A")){
    usePermulations = TRUE
    useCategoricalPermulations = TRUE
  }
  usePermulations = as.logical(usePermulations)
  if(is.na(usePermulations)){
    usePermulations = TRUE
    message("Use Permualtions value not interpretable as logical. Did you remember to capitalize? Using TRUE.")
  }
}else{
  message("Use Permulations value not specified, using TRUE.")
}

#Import subdirectory
if(!any(is.na(cmdArgImport('s', T)))){
  useSubdirectory = TRUE
  subdirectoryValueList = cmdArgImport('s', F)
  #message(paste("Using subdirectory", subdirectoryValue, "."))
  if(length(subdirectoryValueList ==1)){outputFolderName = paste(outputFolderName, subdirectoryValueList[1], "/", sep=""); subdirectoryValue = subdirectoryValueList[1]}
  
}else{
  message("No subdirectory specified.")
  
  
  folderNames = basename(list.dirs(outputFolderName, full.names = TRUE, recursive = FALSE))
  if(length(folderNames) >0 ){
    useSubdirectory = TRUE
    subdirectoryValueList = folderNames
    subdirectoryValueList = subdirectoryValueList[!(subdirectoryValueList == "Alternates")]
    message("Folders not specified, using all folders but Alternates.")
  }else{
    message("Using main directory")
  }
  
  
}


#Import permulation p-value Override 
if(usePermulations){
  if(!is.na(cmdArgImport('z', T))){
    usePermulationPValOverride = TRUE
    permulationPValOverride = cmdArgImport('z', F)
    message("Using Manually specified permulation p-value file.")
  }else{
    if(!useCategoricalPermulations){
      message("No Permulation pValue override, using standard CombinedBinaryPermulationsPValue.rds.")
    }
  }
}

#Import correlation file override
if(!is.na(cmdArgImport('y', T))){
  useCorrelationOverride = TRUE
  correlationOverride = cmdArgImport('y', F)
  message("Using Manually specified correlation file.")
}else{
  message("No corelation override specified.")
}




#---------------------------------------------------------------------
# --- RunMergedAnalysis.R  --- 
# --------------------------------------------------------------------

# Defaults
mainTreesLocation = "data/zoonomiaAllMammalsTrees.rds"  #this is the standard location on Sol 
minSpValue = 10
toMerge = NULL
mergeValue = NULL
runMerge= F

#MainTrees Location
if(!is.na(cmdArgImport('m', T))){
  mainTreesLocation = cmdArgImport('m', F)
}else{
  message("No maintrees arg, using default")
}
#min.sp Value
if(!is.na(cmdArgImport('l', T))){
  minSpValue = cmdArgImport('l', F)
}else{
  message("min.sp not specified, using 10")
}

#toMerge
if(!all(is.na(cmdArgImport('e', T)))){
  toMerge = cmdArgImport('e', F)
  runMerge = T
}else{
  message("No phenotypes set to merge. Script will effectively do nothing; if you are not running a merged phenotype, this is fine.")
  runMerge = F
}

#mergeValue
if(!is.na(cmdArgImport('t', T))){
  mergeValue = cmdArgImport('t', F)
}else{
  if(!runMerge){
    message("No phenotypes set to merge. Script will effectively do nothing; if you are not running a merged phenotype, this is fine.")  
  }else{
    stop("ERROR: Attempting to run merge, but the merged phenotype is unnamed. Exiting.")
  }
}


#---------------------------------------------------------------------
# --- RunRERAndCorrelation.R  --- 
# --------------------------------------------------------------------

# Defaults
mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"  #this is the standard location on Sol 
binaryPhenotypeTreeLocation = NULL
speciesFilter = NULL
phenotypeStyle = "continuous"
continousMetric = "diff"
validMetrics = c("diff", "mean", "last")
minSpValue = 10


#MainTrees Location
if(!is.na(cmdArgImport('m', T))){
  mainTreesLocation = cmdArgImport('m', F)
}else{
  message("No maintrees arg, using default")
}

#Phenotype Style
if(!is.na(cmdArgImport('a', T))){
  phenotypeStyle = cmdArgImport('a', F)
  #Convert the various input options; so that the style can be specified with several names 
  if(phenotypeStyle == "b" | phenotypeStyle == "B" | phenotypeStyle == "binary"){phenotypeStyle = "Binary"}
  if(phenotypeStyle == "c" | phenotypeStyle == "C" | phenotypeStyle == "continuous"){phenotypeStyle = "Continuous"}
  if(phenotypeStyle == "g" | phenotypeStyle == "G" | phenotypeStyle == "categorical"){phenotypeStyle = "Categorical"}
}else{
  message("No phenotype style specified, attempting to detect automatically.")
  
  #Check for categorical 
  phenotypeTreeFilename = paste(outputFolderName, filePrefix, "Categorical", "Tree.rds", sep="")
  if(file.exists(phenotypeTreeFilename)){
    message("Found information for Categorical phenotype, using Categorical phenotype")
    phenotypeStyle = "Categorical"
  }else{
    phenotypeTreeFilename = paste(outputFolderName, filePrefix, "Continuous", "PhenotypeVector.rds", sep="")
    if(file.exists(phenotypeTreeFilename)){
      message("Found information for continuous phenotype, using continuous phenotype")
      phenotypeStyle = "Continuous"
    }else{
      phenotypeTreeFilename = paste(outputFolderName, filePrefix, "Binary", "Tree.rds", sep="")
      if(file.exists(phenotypeTreeFilename)){
        message("Found information for binary phenotype, using binary phenotype")
        phenotypeStyle = "binary"
      }else{
        message("No phenotyep specified and no correct tree found.")
        stop()
      }
    } 
  }
  message("PhenotypeStyle not specified, using continuous")
}

#phenotype tree location
if(!is.na(cmdArgImport('y', T))){
  phenotypeTreeLocation = cmdArgImport('y', F)
}else{                                                                        #See if a pre-made tree for this prefix and style exists 
  if(phenotypeStyle == "Continuous"){
    phenotypeTreeFilename = paste(outputFolderName, filePrefix, phenotypeStyle, "PhenotypeVector.rds", sep="")
  }else{
    phenotypeTreeFilename = paste(outputFolderName, filePrefix, phenotypeStyle, "Tree.rds", sep="")
  }
  
  if(file.exists(paste(phenotypeTreeFilename))){                              #if so, use it                
    phenotypeTreeLocation = phenotypeTreeFilename                     
    message("Pre-made Phenotype tree found, using pre-made tree.")
  }else{
    #paste("THIS IS AN ISSUE MESSAGE; SPECIFY PHENOTYPE TREE")
    stop("THIS IS AN ISSUE MESSAGE; NO PHENOTYPE TREE/VECTOR FOUND. ENSURE PHENOTYPE TREE/VECTOR FOR SELECTED STYLE AVAILABLE AND PREFIX IS CORRECT.")
  }
}

#Species Filter   
if(!is.na(cmdArgImport('z', T))){
  speciesFilter = cmdArgImport('z', F)
}else{ #See if a pre-made filter for this prefix exists 
  speciesFilterFileName = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #Make the name of the location a pre-made filter would have to test for it
  if (file.exists(paste(speciesFilterFileName))){                  
    speciesFilter = readRDS(speciesFilterFileName)                       #if so, use it 
    message("Pre-made filter found, using pre-made filter.")
  }else{                                                    
    message("No speciesFilter arg, using all species")                           #if not, use no filter
  }
}

#Continuous Metric 
if(!is.na(cmdArgImport('k', T))){
  continousMetric = cmdArgImport('k', F)
  if(phenotypeStyle == "Continuous"){ 
    if(!any(continousMetric == validMetrics)){
      stop("Contious metric is not a valid option. Use 'diff', 'mean', or 'last'")
    }
  }
}else{                                                                        #If a continuous phenotype, report using Diff
  if(phenotypeStyle == "Continuous"){
    message("No continuous metric specified, using diff")
  }
}

#min.sp Value
if(!is.na(cmdArgImport('l', T))){
  minSpValue = cmdArgImport('l', F)
}else{
  message("min.sp not specified, using 10")
}



#---------------------------------------------------------------------
# --- MakeContinousPhenotypeTree.R  --- 
# --------------------------------------------------------------------

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
if(!is.na(cmdArgImport('m', T))){
  mainTreesLocation = cmdArgImport('m', F)
}else{
  message("No maintrees arg, using default")
}
#read in the tree based on filetype extension


#spreadsheet File
if(!is.na(cmdArgImport('td', T))){
  spreadSheetLocation = cmdArgImport('td', F)
}else{
  message("Using Data/mergedData.csv spreadsheet")
}

#Phenotype Column
if(!is.na(cmdArgImport('tp', T))){
  annotColumn = cmdArgImport('tp', F)
}else{
  stop("THIS IS AN ISSUE MESSAGE; SPECIFY PHENOTYPE COLUMN")
}

#Tip Column
if(!is.na(cmdArgImport('tt', T))){
  tipCol = cmdArgImport('tt', F)
}else{
  message("Tip Column not specified, using 'tipName'.")
}

#Name Column
if(!is.na(cmdArgImport('tn', T))){
  nameColumn = cmdArgImport('tn', F)
}else{
  message("Name Column not specified, using default list.")
}

#Continuous Metric 
if(!is.na(cmdArgImport('k', T))){
  continousMetric = cmdArgImport('k', F)
  if(!any(continousMetric == validMetrics)){
    stop("Contious metric is not a valid option. Use 'Diff', 'mean', or 'last'")
  }
}else{                                                                        #If a continuous phenotype, report using Diff
  message("No continuous metric specified, using Diff")
}

#Screen Column
if(!is.na(cmdArgImport('ts', T))){
  useScreen = T
  screenColumn = cmdArgImport('ts', F)
}else{
  message("No screen column used.")
}



#Pruning cutoff
if(!is.na(cmdArgImport('pm', T))){
  usingPruning = T
  usingAutoPruning = T
  pruningCutoff = cmdArgImport('pm', F)
}else{
  message("Pruning Cutoff not specified, not pruning tree.")
}

#PruningPrefrenceColumn 
if(!is.na(cmdArgImport('pk', T))){
  pruningPrefrenceColumn = cmdArgImport('pk', F)
}else{
  if(usingPruning){message("No pruning prefrence column specified")}
}

#ManualPruningProtections
if(!all(is.na(cmdArgImport('ps', T)))){
  manualPruningProtections = cmdArgImport('ps', F)
}else{
  if(usingPruning){message("No manually protected species specified")}
}

#ManualPruningSpecies
if(!all(is.na(cmdArgImport('pd', T)))){
  manualPruningSpecies = cmdArgImport('pd', F)
  usingPruning = T
}else{
  if(usingPruning){message("No manually pruned species specified")}
}




#---------------------------------------------------------------------
# --- MakeCategoricalPhenotypeTree.R (new) --- 
# --------------------------------------------------------------------


# Defaults
mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"
spreadSheetLocation = "Data/manualAnnotationsSheet.csv"
annotColumn = NULL
categoryList = NULL
useScreen = F
screenColumn = NULL
modelType = "ER"
ancestralTrait = NULL
substitutions = NULL
tipCol = "tipName"
nameColumn = NULL
usingPruning = F
usingAutoPruning = F
manualPruningProtections = NULL
pruningPrefrenceColumn = NA
pruningProtection = F
manualPruningSpecies = NULL
useLiam = F
generateAlternates = F
numberOfAlternates = 100
colorSet = c( "darkgreen", "darkblue","black", "red", "yellow")

#MainTrees Location
if(!is.na(cmdArgImport('m', T))){
  mainTreesLocation = cmdArgImport('m', F)
}else{
  message("No maintrees arg, using default")
}


#spreadsheet File
if(!is.na(cmdArgImport('td', T))){
  spreadSheetLocation = cmdArgImport('td', F)
}else{
  message("Using Data/manualAnnotationsSheet.csv spreadsheet")
}

#Phenotype Column
if(!is.na(cmdArgImport('tp', T))){
  annotColumn = cmdArgImport('tp', F)
}else{
  stop("THIS IS AN ISSUE MESSAGE; SPECIFY PHENOTYPE COLUMN")
}

#Tip Column
if(!is.na(cmdArgImport('tt', T))){
  tipCol = cmdArgImport('tt', F)
}else{
  message("Tip Column not specified, using 'tipName'.")
}



#Name Column
if(!is.na(cmdArgImport('tn', T))){
  nameColumn = cmdArgImport('tn', F)
}else{
  message("Name Column not specified, using default list.")
}

#Category list 
if(!all(is.na(cmdArgImport('tc', T)))){
  categoryList = cmdArgImport('tc', F)
}else{
  message("Using all available categories in phenotype column.")
  categoryList = NULL
}

#Screen Column
if(!is.na(cmdArgImport('ts', T))){
  useScreen = T
  screenColumn = cmdArgImport('ts', F)
}else{
  message("No screen column used.")
}

#Substitution list 
if(!is.null(cmdArgImport('tr', T))){
  substitutions = cmdArgImport('tr', F)
}else{
  message("No substitutions provided")
}

#Merge list 
if(!is.null(cmdArgImport('tm', T))){
  mergeOnlys = cmdArgImport('tm', F)
}else{
  message("No merges provided")
}

#Ancestral Trait
if(!is.na(cmdArgImport('ta', T))){
  ancestralTrait = cmdArgImport('ta', F)
}else{
  message("No ancestral trait specified, using NULL")
}

#Generate alternates
if(!is.na(cmdArgImport('tl', T))){
  generateAlternates = as.logical(cmdArgImport('tl', F))
  if(generateAlternates){message("Generating Alternates")}
}else{
  message("Not Generating Alternates.")
}

#Number of Alternates
if(generateAlternates){
  if(!is.na(cmdArgImport('ta', T))){
    numberOfAlternates = cmdArgImport('ta', F)
    numberOfAlternates = as.numeric(numberOfAlternates)
  }else{
    message("NUmber of alternates not specified, using 100.")
  }
}

#Color pallete
if(!all(is.na(cmdArgImport('tb', T)))){
  colorSet = cmdArgImport('tb', F)
  message("using specified colorset")
}else{
  message("Using Default palette Generating Alternates.")
}

#use Liam Infrence
if(!is.na(cmdArgImport('tz', T))){
  useLiam = as.logical(cmdArgImport('tz', F))
  message("Using Liam infrence -- DO NOT DO THIS UNLESS USING A SPECIFIC LIAM PHENOTYPE")
}else{
  message("Liam infrence not used.")
}

#Model Type
if(!is.na(cmdArgImport('q', T))){
  modelType = cmdArgImport('q', F)
}else{
  message("No model specified, using ER.")
}



#Pruning cutoff
if(!is.na(cmdArgImport('pm', T))){
  usingPruning = T
  usingAutoPruning = T
  pruningCutoff = cmdArgImport('pm', F)
  pruningCutoff = as.numeric(pruningCutoff)
}else{
  message("Pruning Cutoff not specified, not pruning tree.")
}

#PruningPrefrenceColumn 
if(!is.na(cmdArgImport('pk', T))){
  pruningPrefrenceColumn = cmdArgImport('pk', F)
}else{
  if(usingPruning){message("No pruning prefrence column specified")}
}

#ManualPruningProtections
if(!all(is.na(cmdArgImport('ps', T)))){
  manualPruningProtections = cmdArgImport('ps', F)
}else{
  if(usingPruning){message("No manually protected species specified")}
}

#ManualPruningSpecies
if(!all(is.na(cmdArgImport('pd', T)))){
  manualPruningSpecies = cmdArgImport('pd', F)
  usingPruning = T
}else{
  if(usingPruning){message("No manually pruned species specified")}
}
