# -- Usage:
# This script can be used to run RER calculations and phenotype correlations with Binary, continuous, or categorical phenotypes.
# It inputs a phenotype tree made by the appropriate script for the style, and outputs an RER file, a Paths file, and a Correlation file. 
# For continuous phenotypes, running both "diff" and "mean" using different file prefixes is recommended. 

# -- Libraries 
source("Src/Reu/loadLibraries.R")
loadLibraries() # Note that this library load includes cluster run detection. 
#-------------------------------
# --- Command arguments list ---
# - Required - 
# r = filePrefix                                                               REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# m = mainTreeFilename.txt or .rds                                             REQUIRED: This sets the location of the maintrees file

# - Required, Defaults likely acceptable -
# f = forceUpdate <T or F>                                                     This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# k = < "diff" or "mean" or "last" >                                           This is used for continuous traits, to determine if the metic should be the difference between the nodes (diff), the mean(mean of the two nodes), or last(the downstream value). Note that Mean and Last are not phylogenetically independent, and do not have downstream processing. 
# a = < ["b" or "binary"] or ["c" or "continuous"] or ["g" or "categorical"]>  This prefix is used to set the Analysis Style type of phenotype being supplied; it will detect automatically if not specified.
# - Optional -
# h = <T>                                                                      This ARG causes the script to return its ARG list, instead of running the script. 
# l = <min.sp value>                                                           This sets the min.sp value to be used in the correlation. 
# - Overrides - 
# y = speciesFilterText                                                        OPTIONAL OVERRIDE: This can be used to manually specify a species filer; leave blank for automatic
# z = phenotypeTreeFilename.txt or .rds                                        OPTIONAL OVERRIDE: This can be used to manually override the phenotype tree being used. For continuous analyses, this is the location of the trait vector.
# --- End arguments list ---
#-------------------------------
# -- Section to include manually-specified args --
if(!exists("args", envir = environment(), inherits = FALSE)){
  #This is a space where you can manually write in arguments for the script. 
  #The above if statement ensures these lines will not override any existing args if the script is run whole-sale. 
  
  args = c('r=ComplexDietHighSugarAllSpeciesTest', 'f=T', 'm=data/zoonomiaAllMammalsTrees.rds', 
           'a=g', 'l=10'
           
  )
  
}
#-------------------------------

if(!exists("NoArgUpdate")){ #This line allows other scripts to run this code without regenerating args, if that script creates a "NoArgUpdate" variable.
  # --- Standard start-up code; identical in all scripts---
  source("Src/Reu/StandardizedStartupCode.R") #This runs the standardized code block at this location, ensuring identity across scripts. 
  #This code reads the "report args" help ARG, and returns the args if it is active. 
  #It then reads the filePrefix and sets the output directory to that file prefix, and reads the forceUpdate arg. 

  
  # --- Argument Imports ---
  { # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
    #Defaults
    mainTreesLocation = "Results/mainTreesFromAlignments.rds"
    continousMetric = "diff"
    phenotypeStyle = "categorical"
    speciesFilter = NULL
    phenotypeTreeLocation = NULL
    minSpValue = 10
    
    #Load from args 
    source("Src/Reu/ArgumentImport.R")
    
    validMetrics= c("diff", "mean", "last")
    if(phenotypeStyle == "Continuous" && !any(continousMetric == validMetrics)){
      stop("Contious metric is not a valid option. Use 'diff', 'mean', or 'last'")
    }
  }
  
  
  #Generate the output filenames here, so that it can be overwritten by other scripts sourcing this one. 
  correlationFileNameValue = paste(outputFolderName, filePrefix, "CorrelationFile", sep= "")
  
  if(phenotypeStyle == "Binary"){
    phenotypeVectorFilenameValue = paste(outputFolderName, filePrefix, "phenotypeVector.rds", sep="")
  }else if(phenotypeStyle == "Categorical"){
    phenotypeVectorFilenameValue = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #select the phenotype vector based on prefix
    pairwiseCorrelationFileNameValue = paste(outputFolderName, filePrefix, "PairwiseCorrelationFile", sep= "")
    combinedCategoricalCorrelationFilenameValue = paste(outputFolderName, filePrefix, "CombinedCategoricalCorrelationFile", sep= "") # make this file for later functions that want it in combo
    correlationsOverallFilenameValue = paste(outputSubdirectory, filePrefix, "OverallCorrelationFile.rds", sep= "")
  }
}
#                   ------- Code Body -------- 

# ---------------------------
# ----- SCRIPT MAIN CODE ----
# ---------------------------


# -- Read in the trees --
#MainTrees
if(file_ext(mainTreesLocation) == "rds"){
  if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
}else{
  if(!exists("mainTrees")){mainTrees = readTrees(mainTreesLocation)} 
}

#Phenotype tree
if(file_ext(phenotypeTreeLocation) == "rds"){                                   #if the tree is an RDS file
  phenotypeTree = readRDS(phenotypeTreeLocation)                                #Read as RDS
}else{                                                                          #Otherwise
  phenotypeTree = readTrees(phenotypeTreeLocation)                              #read as text
}

#Species filter
if(all(is.null(speciesFilter))){ # If the species filter is meant to be empty, that is, all of the species should be used
  speciesFilter = mainTrees$masterTree$tip.label #include all of the species on the tree in the filter 
}



# --- RERs ---

RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix

if(!file.exists(paste(RERFileName)) | forceUpdate){                             #if it does not exist, or update is forced 
  gc()
  RERObject = getAllResiduals(mainTrees, useSpecies = speciesFilter, plot = F)  #Calculate the RERs
  saveRDS(RERObject, file = RERFileName)                                        #Save them
}else{                                                                          #Otherwise
  RERObject = readRDS(RERFileName)                                              #Use the existing ones
}

# --- PATHS ---

pathsFileName = paste(outputFolderName, filePrefix, phenotypeStyle, "PathsFile.rds", sep= "") #Set a filename for the pathss based on the prefix and style
if(!file.exists(paste(pathsFileName)) | forceUpdate){                           #if it does not exist, or update is forced 
  if(phenotypeStyle == "Binary"){                                               #If binary 
    pathsObject = tree2Paths(phenotypeTree, mainTrees, binarize=T, useSpecies = speciesFilter)  #make path with binary function
  }else if(phenotypeStyle == "Continuous"){                                      #If continous
    continousTraitVector = readRDS(phenotypeTreeLocation)                         #repurposing the phenotype tree argument for this, as it is the equivalent
    pathsObject = char2Paths(phenotypeTree, mainTrees, metric = continousMetric)
  } else if(phenotypeStyle == "Categorical"){                                   #if categorical
      pathsObject = tree2Paths(phenotypeTree, mainTrees, useSpecies = speciesFilter, categorical = TRUE) #do not binarize; the categorical data is already contained in the phenotype tree.
  }
  saveRDS(pathsObject, file = pathsFileName)                                    #Save the paths
}else{
  pathsObject = readRDS(pathsFileName)                                          #If the file already exists, use the existing one.
}


# --- CORRELATION ---
correlationFileName = correlationFileNameValue #Make a correlation filename based on the prefix, unless overwritten.

if(phenotypeStyle == "Binary"){                                                 #If binary
  correlation = correlateWithBinaryPhenotype(RERObject, pathsObject, min.sp =minSpValue)#Correlate with binary phenotype
  
  #Generate a phenotype vector 
  fgEdgeObjects = phenotypeTree$edge[which(phenotypeTree$edge.length>=1) ,]                                        #Make an object of the edges in the foreground. This is used as opposed to just referencing the tree directly to allow for "walking" in the final loop of the code
  foregroundNodes = which(1:length(phenotypeTree$tip.label) %in% as.vector(fgEdgeObjects))
  foregroundSpecies = phenotypeTree$tip.label[foregroundNodes]
  phenotypeVector = c(0,0);length(phenotypeVector) = length(phenotypeTree$tip.label);phenotypeVector[] = 0 
  names(phenotypeVector) = phenotypeTree$tip.label
  phenotypeVector[(names(phenotypeVector) %in% foregroundSpecies)] = 1
  phenotypeVector = phenotypeVector[names(phenotypeVector) %in% speciesFilter]
  phenotypeVectorFilename = phenotypeVectorFilenameValue
  saveRDS(phenotypeVector, file = phenotypeVectorFilename)
  
}else if(phenotypeStyle == "Continuous"){                                        #if continuous
  correlation = correlateWithContinuousPhenotype(RERObject, pathsObject, min.sp = minSpValue)
  
} else if(phenotypeStyle == "Categorical"){                                     #if categorical
  categoricalCorrelation = correlateWithCategoricalPhenotype(RERObject, pathsObject, min.sp = minSpValue, min.pos = 2) #Calculate with categorical, min 2 species per category 
  overalCategorical = categoricalCorrelation[[1]]                               #select the results relating to overall difference between all categories
  correlation = overalCategorical                                               # and classify it as the main correlation file
  
  #process the pairwise outputs
  pairwiseCategorical = categoricalCorrelation[[2]]                             #select the group of pairwise comparisons
  
  phenotypeVectorFilename = phenotypeVectorFilenameValue
  phenotypeVector = readRDS(phenotypeVectorFilename)                            #load in the phenotype vector 
  categories = map_to_state_space(phenotypeVector)                              #and use it to connect branch lengths to phenotype name
  categoryNames = categories$name2index                                         #store the length-phenotype connection
  
  pairwiseTableNames = names(pairwiseCategorical)                               #Prepare to repalce the number-number titles with phenotype-phenotype titles
  for(i in 1:length(categoryNames)){                                            #for each phenotype
    pairwiseTableNames= gsub(i, names(categoryNames)[i], pairwiseTableNames)                        #replace the number with the phenotype name  
  }
  names(pairwiseCategorical) = pairwiseTableNames                               #update the dataframe titles
  
  pairwiseCorrelationFileName = pairwiseCorrelationFileNameValue
  write.csv(pairwiseCategorical, file= paste(pairwiseCorrelationFileName, ".csv", sep=""), row.names = T, quote = F) #save the correlations as a csv
  saveRDS(pairwiseCategorical, paste(pairwiseCorrelationFileName, ".rds", sep="")) #and as an rds 
  
  combinedCategoricalCorrelationFilename = combinedCategoricalCorrelationFilenameValue # make this file for later functions that want it in combo
  saveRDS(categoricalCorrelation, paste(combinedCategoricalCorrelationFilename, ".rds", sep="")) #and as an rds 
  
  #save the outputs to subdirectories 
  outputSubdirectoryNoslash = paste(outputFolderName, "Overall", sep = "")
  if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
    dir.create(outputSubdirectoryNoslash)
  }
  outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
  
  correlationsOverallFilename = correlationsOverallFilenameValue
  saveRDS(categoricalCorrelation[[1]], correlationsOverallFilename)
  
  for(i in 1:length(pairwiseTableNames)){
    pairwiseTableNames= gsub(" ", "", pairwiseTableNames)
    
    outputSubdirectoryNoslash = paste(outputFolderName, pairwiseTableNames[i], sep = "")
    if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
      dir.create(outputSubdirectoryNoslash)
    }
    outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
    
    correlationsPairFilename = paste(outputSubdirectory, filePrefix, pairwiseTableNames[i], "CorrelationFile",".rds", sep= "")
    saveRDS(categoricalCorrelation[[2]][[i]], correlationsPairFilename)
  }
  
}

write.csv(correlation, file= paste(correlationFileName, ".csv", sep=""), row.names = T, quote = F) #Save correlations as csv
saveRDS(correlation, paste(correlationFileName, ".rds", sep=""))                          #and as an rds 

