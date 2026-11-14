# -- Usage:
# This script runs the pipeline in this package, by sourcing the various included .R files. 
# In theory, it can be the only file used in an analysis. 
# It includes places for the user to supply the args required to run the scripts. 
# -- Libraries 
source("Src/Reu/loadLibraries.R")
loadLibraries(useVisualizePackages = T, loadAllScripts = T) # Note that this library load includes cluster run detection. 
globalHasLoadedLibrariesAlready = T #because we just ran the most extensive version of library loading, we use a global variable to allow all other scripts to skip running library loading. 

# Here included is an ABRIDGED version of the args used to control the pipeline. 
# It includes all of the args a used is likely to want to change, when running the analysis via this script. 
# This notably lacks arguments for permulations (not run via this script).
# A full list of of the ARGs used, as well as which script those args are used in, can be found at "Src/Reu/ArgImports.R
#-------------------------------
# --- Command arguments list ---
# --- Combined Pipeline Args ---
# - Required -
# r  = filePrefix                                         This is a prefix used to organize and separate files by analysis run. Always required. 
# m  = mainTreeFilename.txt or .rds                       This sets the location of the maintrees file
#Analysis
# a = <["b" or "binary"] or ["c" or "continuous"] or ["g" or "categorical"]>   This prefix is used to set the Analysis Style type of phenotype being supplied; it will detect automatically if not specified.
# g = c("gene", "sets", "file", "paths")                  This sets the gene sets (file names) that the script should be run on, or "gene" for use on genes. 
#Tree Creation
# td = spreadSheetFilename.csv                            This sets the spreadsheet to read the data from 
# tp = "phenotypeColumn"                                  This is the column with the phenotype information
# tt = "tipColumn"                                        This sets the column with the tip names as they appear in the maintrees file.

# - Required, Defaults likely acceptable -
# f  = forceUpdate <T or F>                               This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# q = <ER or SYM or ARD>                                  This is used for categorical traits, and sets the model type used to estimate ancestral branches.
# k  = < "Diff" or "mean" or "last" >                     This is used for continuous traits, to determine if the metic should be the difference between the nodes (diff), the mean(mean of the two nodes), or last(the downstream value). Note that Mean and Last are not phylogenetically independent, and do not have downstream processing. 
# s = "subdirectoryName"                                  This is used to specify a subdirectory for the analysis to be run in. Primarily used for the components of categorical results. 

# - Optional - 
# h  = <T>                                                This ARG causes the script to return its ARG list, instead of running the script. 
# l = <min.sp value>                                      This sets the min.sp value to be used in the correlation. 
# c = significance cutoff                                 This sets the cutoff for inclusion in the significant columns
# ts = "screenCollumn"                                    This is a collumn which must have a value of 1 for the species to be included. 
#Categorical
# tc = <c("nameOfCategory1","nameOfCategory2")>           This is used for categorical traits, and sets a list of category names to be used in the analysis. Species with values not in this list will not be included in the analysis. Include ALL categories to be included (even those merged later).  
# tr = list(c("replace1", "with1"),c("replace2, with2"))  This is used for categorical traits, and combined multiple entries to a new value. The new value can be an existing category, or a new label. 
# e = c("phenotypes", "to", "merge)                       This sets which phenotype should be combined into a single phenotype in a merged analysis. 
# t = c("nameOfMergedPhenotype")                          This sets the name of the now merged phenotype
# tl = <T or F>                                           This sets if alternate species sets should be created. Currently supports only categorical phenotypes. 
# ta = <number of alternates>                             This sets the number of alternate species sets to generate, if alternate species sets are generated. 
# i = usedAlternates                                      This sets which alternate(s) to run on for parallelization. if left blank, will run on all alternates. Set to "all" to allow a local script to process all alternates.
# o = alternateRatio <decimal>                            This sets the percentage of alternates which must be significant for a result to be considered robust to alternate species selection. 
# - Pruning - 
# pm = <minimum branch length>                            This sets the minimum branch length for terminal branches in the master tree. Branches shorter than this will be removed. 
# ps = "c('unprunedtip1', 'unprunedtip2')"                This allows you to add a list of specific tips to not be dropped during pruning. Must use the tip name, not common name. 
# --- End arguments list ---
#-------------------------------
# -- Section to include manually-specified args --
if(!exists("args", envir = environment(), inherits = FALSE)){
  #This is a space where you can manually write in arguments for the script. 
  #The above if statement ensures these lines will not override any existing args if the script is run whole-sale. 
  #Additionally, if no args exist and the script is being run locally, it will load args from the control file.
  
}
# --- Standard start-up code; identical in all scripts---
source("Src/Reu/StandardizedStartupCode.R") #This runs the standardized code block at this location, ensuring identity across scripts. 
#This code reads the "report args" help ARG, and returns the args if it is active. 
#It then reads the filePrefix and sets the output directory to that file prefix, and reads the forceUpdate arg

# --- Argument Imports ---
{ # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
  CurrentScript = "ExecutePackage.R"
  
  # Defaults
  phenotypeStyle = NULL
  generateAlternates = NULL
  usedAlternates = NULL
  
  
  #Load from args 
  source("Src/Reu/ArgumentImport.R")
  
  phenotypeStyleExecute = phenotypeStyle
  generateAlternatesExecute = generateAlternates
  usedAlternatesExecute = usedAlternates
  rm(phenotypeStyle)
  rm(generateAlternates)
  rm(usedAlternates)

}


# ---------------------------
# ----- SCRIPT MAIN CODE ----
# ---------------------------

#---------------------------------------------------------------------
# --- Phenotype Tree Creation  --- 
# --------------------------------------------------------------------
# Here included is an ABRIDGED version of the args used to control the phenotype tree portion of the pipeline. 
# It includes all of the args a used is likely to want to change, when running the analysis via this script. 
# A full list of of the ARGs used, as well as which script those args are used in, can be found at "Src/Reu/ArgImports.R


if(phenotypeStyleExecute == "Categorical"){
  source("Src/Loc/MakeCategoricalPhenotypeTree.R")
}else if(phenotypeStyleExecute == "Continuous"){
  source("Src/Loc/MakeContinousPhenotypeTree.R")
}else if(phenotypeStyleExecute == "Binary"){
  source("Src/Loc/MakeBinaryPhenotypeTree.R")
}else{
  stop("Invalid Phenotype Style supplied.")
}



#---------------------------------------------------------------------
# --- RER Analysis Section  --- 
# --------------------------------------------------------------------

#Run the RER generation, and gene correlation 
source("Src/Loc/RunRERAndCorrelation.R")

#Run the Merged Analysis Script. Only applies to categorical phenotypes, but if not applicable it safely does nothing. 
source("Src/Loc/RunMergedAnalysis.R")

#Run gene set enrichment analysis. If no genesets are provided, it safely does nothing. 
source("Src/Loc/RunEnrichmentAnalysis.R")

#If running a categorical phenotype, combine the pairwise results into a single table. 
if(phenotypeStyleExecute == "Categorical"){
  source("Src/Loc/CategoricalCombinePairwiseResults.R")
}

# --- Alternate Species Sets ---
if(generateAlternatesExecute){
  #Run RER and Correlation on the alternates
  source("Src/Loc/ProcessAlternates.R")
  
  #Run Enrichments for the alternates
  source("Src/Loc/RunAlternateEnrichments.R")
  
  #Combine the alternates
  source("Src/Loc/CombineAlternates.R")
  
  #If running a categorical phenotype, add the alternates to the combined files 
  if(phenotypeStyleExecute == "Categorical"){
    source("Src/Loc/AddAlternatesToCombinations.R")
  }
  
}

