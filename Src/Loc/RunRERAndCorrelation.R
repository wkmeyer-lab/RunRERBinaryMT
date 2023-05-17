# -- Libraries 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge)
library(tools)
source("Src/Reu/cmdArgImport.R")

# -- Usage:
# This text describes the purpose of the script 

# -- Command arguments list
# r = filePrefix                                                               This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                                                                 This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# m = mainTreeFilename.txt or .rds                                             This sets the location of the maintrees file
# p = phenotypeTreeFilename.txt or .rds                                        This can be used to manually override the phenotype tree being used 
# f = speciesFilterText                                                        This can be used to manually specify a species filer; leave blank for automatic
# s = < ["b" or "binary"] or ["c" or "continuous"] or ["g" or "categorical"]>  This prefix is used to set the type of phenotype being supplied

#----------------
# --- Standard start-up code ---
args = commandArgs(trailingOnly = TRUE)
{  # Bracket used for collapsing purposes
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
    paste("Force update not specified, not forcing update")
  }
}

# --- Argument Imports ---
# Defaults
mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"  #this is the standard location on Sol 
binaryPhenotypeTreeLocation = NULL
speciesFilter = NULL
phenotypeStyle = "continuous"

{ # Bracket used for collapsing purposes

  #MainTrees Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using default")
  }
  
  #Phenotype Style
  if(!is.na(cmdArgImport('t'))){
    phenotypeStyle = cmdArgImport('t')
    #Convert the various input options 
    if(phenotypeStyle == "b" | phenotypeStyle == "B" | phenotypeStyle = "binary"){phenotypeStyle = "Binary"}
    if(phenotypeStyle == "c" | phenotypeStyle == "C" | phenotypeStyle = "continuous"){phenotypeStyle = "Continuous"}
    if(phenotypeStyle == "g" | phenotypeStyle == "G" | phenotypeStyle = "categorical"){phenotypeStyle = "Categorical"}
  }else{
    message("PhenotypeStyle not specified, using continuous")
  }
  
  #phenotype tree location
  if(!is.na(cmdArgImport('p'))){
    phenotypeTreeLocation = cmdArgImport('p')
  }else{ #See if a pre-made tree for this prefix and style exists 
    phenotypeTreeFilename = paste(outputFolderName, filePrefix, phenotypeStyle, "Tree.rds", sep="")
    if(file.exists(paste(phenotypeTreeFilename))){   #if so, use it                
      phenotypeTreeLocation = phenotypeTreeFilename                     
      paste("Pre-made Phenotype tree found, using pre-made tree.")
    }else{
      #paste("THIS IS AN ISSUE MESSAGE; SPECIFY PHENOTYPE TREE")
      stop("THIS IS AN ISSUE MESSAGE; SPECIFY PHENOTYPE TREE")
    }
  }
  
  #Species Filter   
  if(!is.na(cmdArgImport('f'))){
    speciesFilter = cmdArgImport('f')
  }else{ #See if a pre-made filter for this prefix exists 
    speciesFilterFileName = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #Make the name of the location a pre-made filter would have to test for it
    if (file.exists(paste(speciesFilterFileName))){                  
      speciesFilter = readRDS(speciesFilterFileName)                       #if so, use it 
      paste("Pre-made filter found, using pre-made filter.")
  }else{                                                    
    message("No speciesFilter arg, using NULL")                           #if not, use no filter
  }
  }
  

}

#                   ------- Code Body -------- 

# -- Read in the trees --
#MainTrees
if(file_ext(mainTreesLocation) == "rds"){
  mainTrees = readRDS(mainTreesLocation)
}else{
  mainTrees = readTrees(mainTreesLocation) 
}
#Phenotype tree
if(file_ext(phenotypeTreeLocation) == "rds"){
  phenotypeTree = readRDS(phenotypeTreeLocation)
}else{
  phenotypeTree = readTrees(phenotypeTreeLocation) 
}

# --- RERs ---

RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")

if(!file.exists(paste(RERFileName)) | forceUpdate){
  RERObject = getAllResiduals(mainTrees, useSpecies = speciesFilter, plot = F)
  saveRDS(RERObject, file = RERFileName)
}else{
  RERObject = readRDS(RERFileName)
}

# --- PATHS ---

pathsFileName = paste(outputFolderName, filePrefix, phenotypeStyle, "PathsFile.rds", sep= "")
if(!file.exists(paste(pathsFileName)) | forceUpdate){
  if(phenotypeStyle == "Binary"){
    pathsObject = tree2Paths(phenotypeTree, mainTrees, binarize=T, useSpecies = speciesFilter)
  }else if(phenotypeStyle == "Continous"){
    stop("This function hasn't been completed")
  } else if(phenotypeStyle == "Categorical"){ 
      pathsObject = tree2Paths(phenotypeTree, mainTrees, useSpecies = speciesFilter)
  }
  saveRDS(pathsObject, file = pathsFileName)
}else{
  pathsObject = readRDS(pathsFileName)
}


# --- CORRELATION ---
correlationFileName = paste(outputFolderName, filePrefix, "CorrelationFile", sep= "")

if(phenotypeStyle == "Binary"){
  correlation = correlateWithBinaryPhenotype(RERObject, pathsObject, min.sp =10)
}else if(phenotypeStyle == "Continous"){
  stop("This function hasn't been completed")
} else if(phenotypeStyle == "Categorical"){ 
  categoricalCorrelation = correlateWithCategoricalPhenotype(RERObject, pathsObject, min.sp = 10, min.pos = 2, )
  overalCategorical = categoricalCorrelation[[1]]
  correlation = overalCategorical
  
  #process the pairwise outputs
  pairwiseCategorical = categoricalCorrelation[[2]]
  
  phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="")
  phenotypeVector = readRDS(phenotypeVectorFilename)
  categories = map_to_state_space(phenotypeVector)
  categoryNames = categories$name2index
  
  pairwiseTableNames = names(pairwiseCategorical)
  for(i in 1:length(categoryNames)){
    gsub(i, names(categoryNames)[i], pairwiseTableNames)
  }
  names(pairwiseCategorical) = pairwiseTableNames
  
  pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "PairwiseCorrelationFile", sep= "")
  write.csv(pairwiseCategorical, file= paste(pairwiseCorrelationFileName, ".csv", sep=""), row.names = T, quote = F)
  saveRDS(pairwiseCategorical, paste(pairwiseCorrelationFileName, ".rds", sep=""))
}

write.csv(correl, file= paste(outputFileName, ".csv", sep=""), row.names = T, quote = F)
saveRDS(correl, paste(outputFileName, ".rds", sep=""))

