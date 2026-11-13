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
# A full list of of the ARGs used, as well as which script those args are used in, can be found at "Src/Reu/ArgIndexes.R
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
  
  
}
# ---------------------------
# ----- SCRIPT MAIN CODE ----
# ---------------------------

# --- Argument Imports ---
{ # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
  
  # Defaults
  phenotypeStyleExecute = "categorical"
  generateAlternatesExecute = F
  
  
  #Phenotype Style (For Execute script)
  if(!is.na(cmdArgImport('a', T))){
    phenotypeStyleExecute = cmdArgImport('a', F)
    #Convert the various input options; so that the style can be specified with several names 
    if(phenotypeStyleExecute == "b" | phenotypeStyleExecute == "B" | phenotypeStyleExecute == "binary"){phenotypeStyleExecute = "Binary"}
    if(phenotypeStyleExecute == "c" | phenotypeStyleExecute == "C" | phenotypeStyleExecute == "continuous"){phenotypeStyleExecute = "Continuous"}
    if(phenotypeStyleExecute == "g" | phenotypeStyleExecute == "G" | phenotypeStyleExecute == "categorical"){phenotypeStyleExecute = "Categorical"}
  }else{
    message("No phenotype style specified, attempting to detect automatically.")
    
    #Check for categorical 
    phenotypeTreeFilename = paste(outputFolderName, filePrefix, "Categorical", "Tree.rds", sep="")
    if(file.exists(phenotypeTreeFilename)){
      message("Found information for Categorical phenotype, using Categorical phenotype")
      phenotypeStyleExecute = "Categorical"
    }else{
      phenotypeTreeFilename = paste(outputFolderName, filePrefix, "Continuous", "PhenotypeVector.rds", sep="")
      if(file.exists(phenotypeTreeFilename)){
        message("Found information for continuous phenotype, using continuous phenotype")
        phenotypeStyleExecute = "Continuous"
      }else{
        phenotypeTreeFilename = paste(outputFolderName, filePrefix, "Binary", "Tree.rds", sep="")
        if(file.exists(phenotypeTreeFilename)){
          message("Found information for binary phenotype, using binary phenotype")
          phenotypeStyleExecute = "binary"
        }else{
          message("No phenotype specified and no correct tree found.")
          message("PhenotypeStyle not specified, using categorical")
        }
      } 
    }
  }
  
  #Generate alternates (for Execute Script)
  if(!is.na(cmdArgImport('tl', T))){
    generateAlternatesExecute = as.logical(cmdArgImport('tl', F))
    if(generateAlternatesExecute){message("Generating Alternates")}
  }else{
    message("Not Generating Alternates.")
  }
  
  
  

}




#---------------------------------------------------------------------
# --- Phenotype Tree Creation  --- 
# --------------------------------------------------------------------
# Here included is an ABRIDGED version of the args used to control the phenotype tree portion of the pipeline. 
# It includes all of the args a used is likely to want to change, when running the analysis via this script. 
# A full list of of the ARGs used, as well as which script those args are used in, can be found at "Src/Reu/ArgIndexes.R


if(phenotypeStyleExecute == "Categorical"){
  source("Src/Loc/MakeCategoricalPhenotypeTree.R")
}elseif(phenotypeStyleExecute == "Continuous"){
  source("Src/Loc/MakeContinousPhenotypeTree.R")
}elseif(phenotypeStyleExecute == "Continuous"){
  source("Src/Loc/MakeBinaryPhenotypeTree.R")
}else{
  stop("Invalid Phenotype Style supplied.")
}



#---------------------------------------------------------------------
# --- RER Analysis Section  --- 
# --------------------------------------------------------------------



# --- Standard start-up code; identical in all scripts---
source("Src/Reu/StandardizedStartupCode.R") #This runs the standardized code block at this location, ensuring identity across scripts. 
#This code reads the "report args" help ARG, and reutrns the args if it is active. 
#It then reads the filePrefix and sets the output directory to that file prefix, and reads the forceUpdate arg. 





