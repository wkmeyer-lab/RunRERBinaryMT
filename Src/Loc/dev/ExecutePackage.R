# -- Usage:
# This script runs the pipeline in this package, by sourcing the various included .R files. 
# In theory, it can be the only file used in an analysis. 
# It includes places for the user to supply the args required to run the scripts. 
# -- Libraries 
source("Src/Reu/loadLibraries.R")
loadLibraries(useVisualizePackages = T, loadAllScripts = T) # Note that this library load includes cluster run detection. 
globalHasLoadedLibrariesAlready = T #because we just ran the most extensive version of library loading, we use a global variable to allow all other scripts to skip running library loading. 


#-------------------------------
# --- Command arguments list ---
# - Required - 
# r = filePrefix                            REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# m = mainTreeFilename.txt or .rds          This sets the location of the maintrees file
# d = spreadSheetFilename.csv               This sets the spreadsheet to read the data from 
# a = "annotCollumn"                        This is the column in the manual annotations spreadsheet to use
# t = "tipColumn"                           This sets the column with the tip names as they appear in the maintrees file.
# n = "nameColumn"                          This sets the column with common names of the species. 
# - Required, Defaults likely acceptable -
# v = <T or F>                              This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# k = < "Diff" or "mean" or "last" >        This is used for continuous traits, to determine if the metic should be the difference between the nodes (diff), the mean(mean of the two nodes), or last(the downstream value). Note that Mean and Last are not phylogenetically independent, and do not have downstream processing. 
# - Optional - 
# h = <T>                                   This ARG causes the script to return its ARG list, instead of running the script. 
# s = "screenCollumn"                       This is a collumn which must have a value of 1 for the species to be included. 
# - Pruning - 
# z = <minimum branch length>                            This sets the minimum branch length for terminal branches in the master tree. Branches shorter than this will be removed. 
# x = "pruningPrefrenceColumn"                           This sets a column, where if the value is 1, the tip will be preferentially kept. If the value is TRUE, the tip will never be pruned.
# y = "c('unprunedtip1', 'unprunedtip2')"                This allows you to add a list of specific tips to not be dropped during pruning. Must use the tip name, not common name. 
# p = "c('prunedtip1', 'prunedtip2')"                    This allows you to manually specify additional branches to be pruned
# - Overrides - 
# --- End arguments list ---
#-------------------------------
# -- Section to include manually-specified args --
if(!exists("args", envir = environment(), inherits = FALSE)){
  #This is a space where you can manually write in arguments for the script. 
  #The above if statement ensures these lines will not override any existing args if the script is run whole-sale. 
  
}
#-------------------------------



# --- Standard start-up code; identical in all scripts---
{  # Bracket used for collapsing purposes
  if(clusterRun)args = commandArgs(trailingOnly = TRUE)
  #Arg help 
  if(isTRUE(as.logical(cmdArgImport('h')))){
    message("Not Running Script, instead returning ARGs for script.")
    printArgDocumentation()
    stop()
  }
  
  #File Prefix
  if(!is.na(cmdArgImport('r'))){
    filePrefix = cmdArgImport('r')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
  }
  
  #  Output Directory 
  if(!dir.exists("Output")){                                      #Make output directory if it does not exist
    dir.create("Output")
  }
  outputFolderNameNoSlash = paste("Output/",filePrefix, sep = "") #Set the prefix sub directory
  if(!dir.exists(outputFolderNameNoSlash)){                       #create that directory if it does not exist
    dir.create(outputFolderNameNoSlash)
  }
  outputFolderName = paste("Output/",filePrefix,"/", sep = "")
  
  #  Force update argument
  forceUpdate = FALSE
  if(!is.na(cmdArgImport('v'))){                                 #Import if update being forced with argument 
    forceUpdate = cmdArgImport('v')
    forceUpdate = as.logical(forceUpdate)
  }else{
    message("Force update not specified, not forcing update")
  }
}

