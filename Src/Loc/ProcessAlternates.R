# -- Usage:
# This script is used to run RERs and Correlations, as well as the Merged Analysis versions of them, for all of the alternate species sets generated.
# It is essentially a combined version of RunRERAndCorrelation and RunMergedAnalysis, but slightly altered to target alternate species sets. 
# Because this can take a very long time, the script is set up to automatically decline running on a local machine if it has more than 10 alternates to process. 
# To avoid this automatic shutdown use one of the 'i' values specified below (e.g. 'all'), which disables this behavior. 

# -- Libraries 
source("Src/Reu/loadLibraries.R")
loadLibraries() # Note that this library load includes cluster run detection. 
#-------------------------------
# --- Command arguments list ---
# - Required - 
# r = filePrefix                                                                REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# m = mainTreeFilename.txt or .rds                                              This sets the location of the maintrees file

# - Required, Defaults likely acceptable -
# f = forceUpdate <T or F>                                                      This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# a = < ["b" or "binary"] or ["c" or "continuous"] or ["g" or "categorical"]>   This prefix is used to set the type of phenotype being supplied
# k = < "diff" or "mean" or "last" >                                            This is used for continuous traits, to determine if the metic should be the difference between the nodes (diff), the mean(mean of the two nodes), or last(the downstream value). Note that Mean and Last are not phylogenetically independent, and do not have downstream processing. 
# - Optional - 
# h = <T>                                                                       This ARG causes the script to return its ARG list, instead of running the script. 
# l = <min.sp value>                                                            This sets the min.sp value to be used in the correlation. 
# i = run instance                                                              This sets which alternate(s) to run on for parallelization. if left blank, will run on all alternates. Set to "all" to allow a local script to process all alternates.
# e = c("phenotypes", "to", "merge)                                             This sets which phenotype should be combined into a single phenotype in a merged analysis. 
# t = c("nameOfMergedPhenotype")                                                This sets the name of the now merged phenotype
# - Overrides - 
# --- End arguments list ---
#-------------------------------
# -- Section to include manually-specified args --
if(!exists("args", envir = environment(), inherits = FALSE)){
  #This is a space where you can manually write in arguments for the script. 
  #The above if statement ensures these lines will not override any existing args if the script is run whole-sale. 
  
  args = c('r=ComplexDietCentralAnalysis', 's=g', 'm=Data/zoonomiaAllMammalsTrees.rds', 'l=170', 'i=1')
  args = c('r=ComplexDietCentralAnalysisSimplify', 's=g', 'm=Data/zoonomiaAllMammalsTrees.rds', 'l=170', 'i=11')
  args = c('r=ComplexDietCentralAnalysisSimplifyStrictPred', 's=g', 'm=Data/zoonomiaAllMammalsTrees.rds', 'l=170', 'i=7')
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
{  # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
  
  # Defaults
  validLocalRunAllAlternates = c("local", "sure", "all")
  validMetrics= c("diff", "mean", "last")
  
  mainTreesLocation = "Results/mainTreesFromAlignments.rds"
  continousMetric = "diff"
  phenotypeStyle = "categorical"
  minSpValue = 10
  usedAlternates = NULL
  alternatesSpecified = F
  toMerge = NULL
  mergeValue = NULL
  runMerge= F
  
  
  #Load from args 
  source("Src/Reu/ArgumentImport.R")
  
  
  if(phenotypeStyle == "Continuous" && !any(continousMetric == validMetrics)){
    stop("Contious metric is not a valid option. Use 'diff', 'mean', or 'last'")
  }
}


# -- Read in the trees --
#MainTrees
if(file_ext(mainTreesLocation) == "rds"){
  if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
}else{
  if(!exists("mainTrees")){mainTrees = readTrees(mainTreesLocation)} 
}

# ---- Make required adjustments to run on alternates properly ---- 

alternateSets = readRDS(paste0(outputFolderName, filePrefix, "AlternatePruningSpecies.rds"))

outputFolderName = paste0(outputFolderName, "Alternates/")
primaryFilePrefix = filePrefix


#                   ------- Code Body -------- 

if(!alternatesSpecified){
  usedAlternates = 1:length(alternateSets)
}

usedAlternates = as.integer(usedAlternates)
for(i in usedAlternates){
  message(i)
  currentSet = alternateSets[[i]]
  alternateFilePrefix = paste0("/Alternates/Alternate", i)
  filePrefix = paste0("Alternate", i, primaryFilePrefix)
  phenotypeTreeLocation = paste(outputFolderName, filePrefix, phenotypeStyle, "Tree.rds", sep="")
  speciesFilterFileName = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #Make the name of the location a pre-made filter would have to test for it
  speciesFilter = readRDS(speciesFilterFileName) 
  correlationFileNameValue = paste(outputFolderName, filePrefix, "CorrelationFile", sep= "")
  phenotypeVectorFilenameValue = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #select the phenotype vector based on prefix
  pairwiseCorrelationFileNameValue = paste(outputFolderName, filePrefix, "PairwiseCorrelationFile", sep= "")
  
  # --- ensure that the species filter and the alternate set match 
  alternateInFilter = all(currentSet %in% speciesFilter)
  filterInAlternate = all(speciesFilter %in% currentSet)
  if(!alternateInFilter & filterInAlternate){
    stop("Species filter and alterate set do not match. Aborting.")
  }
  
  
  
  #Run the RER and Correlation 
  alreadyBlockingArgImports = exists("NoArgUpdate")
  if(!alreadyBlockingArgImports){NoArgUpdate = "Exists"}
  
  source("Src/Loc/RunRERAndCorrelation.R")
  
  if(!alreadyBlockingArgImports){rm(NoArgUpdate)}
  
  # ------
  # -- Run the merged Analysis for the alternate --
  # -----
  if(runMerge){
    alreadyBlockingArgImports = exists("NoArgUpdate")
    if(!alreadyBlockingArgImports){NoArgUpdate = "Exists"}
    
    source("Src/Loc/RunMergedAnalysis.R")
    
    if(!alreadyBlockingArgImports){rm(NoArgUpdate)}
  }
}

