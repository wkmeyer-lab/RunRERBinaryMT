# -- Usage:
# This script runs the pipeline in this package, by sourcing the various included .R files. 
# In theory, it can be the only file used in an analysis. 
# It includes places for the user to supply the args required to run the scripts. 
# -- Libraries 
source("Src/Reu/loadLibraries.R")
loadLibraries(useVisualizePackages = T, loadAllScripts = T) # Note that this library load includes cluster run detection. 
globalHasLoadedLibrariesAlready = T #because we just ran the most extensive version of library loading, we use a global variable to allow all other scripts to skip running library loading. 


# Here included is an ABRIDGED version of the args used to control the analysis portion of the pipeline. 
# It includes all of the args a used is likely to want to change, when running the analysis via this script. 
# A full list of of the ARGs used, as well as which script those args are used in, can be found at "Src/Reu/ArgIndexes.R

#-------------------------------
# --- Combined Analysis Args --- 
# --- Command arguments list ---
# - Required - 
# r = filePrefix                                                                REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# m = mainTreeFilename.txt or .rds                                              REQUIRED: This sets the location of the maintrees file
# a = < ["b" or "binary"] or ["c" or "continuous"] or ["g" or "categorical"]>   This prefix is used to set the Analysis Style type of phenotype being supplied; it will detect automatically if not specified.
# g = c("gene", "sets", "file", "paths")                                        This sets the gene sets (file names) that the script should be run on, or "gene" for use on genes. 

# - Required, Defaults likely acceptable -
# f = forceUpdate <T or F>                                                      This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# k = < "diff" or "mean" or "last" >                                            This is used for continuous traits, to determine if the metic should be the difference between the nodes (diff), the mean(mean of the two nodes), or last(the downstream value). Note that Mean and Last are not phylogenetically independent, and do not have downstream processing. 
# p = <T or F or C or A>                                                        This sets if the code should use permulated or unpermulated values. If set to "C" or "A", will use permulations for categorical values, which are stored in the main file. 
# s = "subdirectoryName"                                                        This is used to specify a subdirectory for the analysis to be run in. Primarily used for the components of categorical results. 
# o = alternateRatio <decimal>                                                  This sets the percentage of alternates which must be significant for a result to be considered robust to alternate species selection. 
# Permulations: q = <ER or SYM or ARD>                                          This sets the model type used to estimate ancestral branches; MUST BE THE SAME AS ONE USED IN MakeCategoricalPhenotypeTree.R

# - Optional -
# h = <T>                                                                       This ARG causes the script to return its ARG list, instead of running the script. 
# l = <min.sp value>                                                            This sets the min.sp value to be used in the correlation. 
# e = c("phenotypes", "to", "merge)                                             This sets which phenotype should be combined into a single phenotype in a merged analysis. 
# t = c("nameOfMergedPhenotype")                                                This sets the name of the now merged phenotype
# c = significance cutoff                                                       This sets the cutoff for inclusion in the significant columns
# Permulations: x = relaxationValue <int, 0-1>                                  This is a value between 0 and 1, the percentage off of exact a permulation is allowed to be. Increases runspeed, but reduces permulation accuracy. Recommend 0.1 if more than 3 categories
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





