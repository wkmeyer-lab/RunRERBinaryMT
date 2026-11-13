# -- Usage:
# This script is used to run gene set enrichment, potentially on multiple gene sets, for all of the alternate species sets generated.
# It is essentially a version of EnrichmentAnalysisNew, but slightly altered to target alternate species sets. 
# Because this can take a very long time, the script is set up to automatically decline running on a local machine if it has more than 10 alternates to process. 
# To avoid this automatic shutdown use one of the 'i' values specified below (e.g. 'all'), which disables this behavior. 

# -- Libraries 
source("Src/Reu/loadLibraries.R")
loadLibraries() # Note that this library load includes cluster run detection. 
#-------------------------------
# --- Command arguments list ---
# - Required - 
# r = filePrefix                              REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# g = c("gene", "sets", "file", "paths")      This sets the gene sets (file names) that the script should be run on, or "gene" for use on genes. 

# - Required, Defaults likely acceptable -
# f = forceUpdate <T or F>                    This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# p = <T or F or C or A>                      This sets if the code should use permulated or unpermulated values. If set to "C", will use permulations for categorical values, which are stored in the main file.
# s = "subdirectoryName"                      This is used to specify a subdirectory for the analysis to be run in. Primarily used for the components of categorical results. 
# - Optional - 
# h = <T>                                     This ARG causes the script to return its ARG list, instead of running the script. 
# i = usedAlternates                          This sets which alternate(s) to run on for parallelization. if left blank, will run on all alternates. Set to "all" to allow a local script to process all alternates.
# - Overrides - 
# y = "correlationFileOverride.rds"           This can be used to manually set the correlation file. Used primarily to target pairwise comparisons of categorical phenotypes.
# z = "permulationPvalueFileLocation.rds"     This is a manual override to specify the script use a specific Permulation p-value file.
                                                #If using permulations in a binary analysis, and using any file other than "CombinedBinaryPermulationsPValue" with no run instance number, it must be specified manually.

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
{  # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
  #Defaults
  geneSets = NULL
  validGeneArguments = c("NULL", "gene", "Gene")
  usePermulations = FALSE
  useCategoricalPermulations = FALSE
  useSubdirectory = FALSE
  subdirectoryValueList = NULL
  usePermulationPValOverride = FALSE
  permulationPValOverride = NULL 
  useCorrelationOverride = FALSE
  correlationOverride = NULL
  usingGo = T
  
  alternatesSpecified = F
  usedAlternates = NULL
  validLocalRunAllAlternates = c("local", "sure", "all")
  
  #Load from args 
  source("Src/Reu/ArgumentImport.R")
  
}

# ---- Make required adjustments to run on alternates properly ---- 

alternateSets = readRDS(paste0("Output/",filePrefix,"/", filePrefix, "AlternatePruningSpecies.rds"))

mainOutputFolderName = outputFolderName
primaryFilePrefix = filePrefix
outputFolderName = paste0(outputFolderName, "Alternates/")



#                   ------- Code Body -------- 

if(!alternatesSpecified){
  usedAlternates = 1:length(alternateSets)
}

for(i in usedAlternates){
  currentSet = alternateSets[[i]]
  alternateFilePrefix = paste0("/Alternates/Alternate", i)
  filePrefix = paste0("Alternate", i, primaryFilePrefix)


  
  alreadyBlockingArgImports = exists("NoArgUpdate")
  if(!alreadyBlockingArgImports){NoArgUpdate = "Exists"}
  
  source("Src/Loc/RunEnrichmentAnalysis.R")
  
  if(!alreadyBlockingArgImports){rm(NoArgUpdate)}
    
}