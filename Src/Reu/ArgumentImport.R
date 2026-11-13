# -- Usage:
# This script performs all of the argument imports for all scripts within this package. 
# By having all such imports in a centralized location, changes are shared across the package. 
# Additionally, argument overlap prevention is natively handled via having a shared import script. 
# This also reduces the dead space at the header of each script. 
source("Src/Reu/loadLibraries.R")
loadLibraries() # Note that this library load includes cluster run detection. 

#The first section Of this script is a combined agruments list, similar to those found at the beginging of each script. 
#Arguments are separated by analysis step. 




#---------------------------------------------------------------------
# --- Combined Tree Generation Control Args  --- 
# --------------------------------------------------------------------
{ # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
#-------------------------------
# --- Command arguments list ---
# - Required -
# r  = filePrefix                                         This is a prefix used to organize and separate files by analysis run. Always required. 
# m  = mainTreeFilename.txt or .rds                       This sets the location of the maintrees file
# td = spreadSheetFilename.csv                            This sets the spreadsheet to read the data from 
# tp = "phenotypeColumn"                                  This is the column with the phenotype information
# tt = "tipColumn"                                        This sets the column with the tip names as they appear in the maintrees file.

# - Required, Defaults likely acceptable -
# f  = forceUpdate <T or F>                               This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# tn = "nameColumn"                                       This sets the column with common names of the species. 
#  Categorical:
# tc = <c("nameOfCategory1","nameOfCategory2")>           This sets a list of category names to be used in the analysis. Species with values not in this list will not be included in the analysis. Include ALL categories to be included (even those merged later).  
#  Continuous:
# k  = < "Diff" or "mean" or "last" >                     This is used for continuous traits, to determine if the metic should be the difference between the nodes (diff), the mean(mean of the two nodes), or last(the downstream value). Note that Mean and Last are not phylogenetically independent, and do not have downstream processing. 

# - Optional - 
# h  = <T>                                                This ARG causes the script to return its ARG list, instead of running the script. 
# ts = "screenCollumn"                                    This is a collumn which must have a value of 1 for the species to be included. 
#  Categorical:
# tr = list(c("replace1", "with1"),c("replace2, with2"))  This combined multiple entries to a new value. The new value can be an existing category, or a new label. 
# tm = list(c("phenotype1", "intophen1"), c("2", "i2"))   This causes combination phenotypes to be merged into the second phenotype, but does not replace standalone phenotypes
# ta = "ancestralTrait"                                   This can be used to set all non-terminal branches to this category. Use be one of the categories in the list. 
# tb = c( "colors", "in ","the", "palette")               This sets the color palette used in the tree graphs. 
# tl = <T or F>                                           This sets if alternate species sets should be created. 
# ta = <number of alternates>                             This sets the number of alternate species sets to generate, if alternate species sets are generated. 
# tz = <T or F>                                           DO NOT USE UNLESS SPECIFIED -- This determines if the liam infrence nodes should be added 
# q  = <ER or SYM or ARD>                                 This sets the model type used to estimate ancestral branches 

# - Pruning - 
# pm = <minimum branch length>                            This sets the minimum branch length for terminal branches in the master tree. Branches shorter than this will be removed. 
# pk = "pruningPrefrenceColumn"                           This sets a column, where if the value is 1, the tip will be preferentially kept. If the value is TRUE, the tip will never be pruned.
# ps = "c('unprunedtip1', 'unprunedtip2')"                This allows you to add a list of specific tips to not be dropped during pruning. Must use the tip name, not common name. 
# pd = "c('prunedtip1', 'prunedtip2')"                    This allows you to manually specify additional branches to be pruned
# --- End arguments list ---
#-------------------------------
}


#---------------------------------------------------------------------
# --- Combined Analysis Control Args  --- 
# --------------------------------------------------------------------
{ # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
#-------------------------------
# --- Command arguments list ---
# - Required - 
# r = filePrefix                                                                REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# m = mainTreeFilename.txt or .rds                                              REQUIRED: This sets the location of the maintrees file
#     #RunRERAndCorrelation.R; RunMergedAnalysis.R; ProcessAlternates.R; CategoricalPermulationsGeneration.R
# g = c("gene", "sets", "file", "paths")                                        This sets the gene sets (file names) that the script should be run on, or "gene" for use on genes. 
#     #EnrichmentAnalysisNew.R; CategoricalCombinePairwiseResults.R; RunAlternateEnrichments.R; CombineAlternates.R; AddAlternatesToCombination.R

# - Required, Defaults likely acceptable -
# f = forceUpdate <T or F>                                                      This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# a = < ["b" or "binary"] or ["c" or "continuous"] or ["g" or "categorical"]>   This prefix is used to set the Analysis Style type of phenotype being supplied; it will detect automatically if not specified.
#     #RunRERAndCorrelation.R; ProcessAlternates.R
# k = < "diff" or "mean" or "last" >                                            This is used for continuous traits, to determine if the metic should be the difference between the nodes (diff), the mean(mean of the two nodes), or last(the downstream value). Note that Mean and Last are not phylogenetically independent, and do not have downstream processing. 
#     #RunRERAndCorrelation.R; ProcessAlternates.R
# p = <T or F or C or A>                                                        This sets if the code should use permulated or unpermulated values. If set to "C" or "A", will use permulations for categorical values, which are stored in the main file. 
#     #EnrichmentAnalysisNew.R; CategoricalCombinePairwiseResults.R; RunAlternateEnrichments.R
# s = "subdirectoryName"                                                        This is used to specify a subdirectory for the analysis to be run in. Primarily used for the components of categorical results. 
#     #EnrichmentAnalysisNew.R; CategoricalCombinePairwiseResults.R; RunAlternateEnrichments.R 
# w = writeData  <T or F>                                                       This sets if the data should be saved, potentially overwriting old combinations
#     #CategoricalCombinePairwiseResults.R; AddAlternatesToCombination.R
# i = usedAlternates                                                            This sets which alternate(s) to run on for parallelization. if left blank, will run on all alternates. Set to "all" to allow a local script to process all alternates.
#     #ProcessAlternates.R; RunAlternateEnrichments.R 
# o = alternateRatio <decimal>                                                  This sets the percentage of alternates which must be significant for a result to be considered robust to alternate species selection. 
#     #AddAlternatesToCombination.R
#   Permulations:
# n = numberOfPermulations                                                      This is the number of permulations to run in the script
#     #CategoricalPermulationsGeneration.R; ContinousPermulations.R
# q = <ER or SYM or ARD>                                                        This sets the model type used to estimate ancestral branches; MUST BE THE SAME AS ONE USED IN MakeCategoricalPhenotypeTree.R
#     #CategoricalPermulationsGeneration.R
# j = runInstanceValuePerms                                                     This is used to generate unique filenames for each instance of the script. Used in parrallelization. Set as j to avoid overlap with alternates.
#     #CategoricalPermulationsGeneration.R
# u = <auto or stationary or flat>                                              This sets the probability of each state at the root of the tree.   
#     #CategoricalPermulationsGeneration.R
# st = <Start number>                                                           This is the permulation number to start at. Used for parrallelization. 
#     #CategoricalPermualtionsCombineAndCalculatePValue.R
# po = outgroupSpecies                                                          This is the outgroup species with the tip name as it appears on mainTrees$masterTree
#     #ContinousPermulations.R


# - Optional -
# h = <T>                                                                       This ARG causes the script to return its ARG list, instead of running the script. 
# l = <min.sp value>                                                            This sets the min.sp value to be used in the correlation. 
#     #RunRERAndCorrelation.R; RunMergedAnalysis.R; ProcessAlternates.R
# e = c("phenotypes", "to", "merge)                                             This sets which phenotype should be combined into a single phenotype in a merged analysis. 
#     #RunMergedAnalysis.R; ProcessAlternates.R
# t = c("nameOfMergedPhenotype")                                                This sets the name of the now merged phenotype
#     #RunMergedAnalysis.R; ProcessAlternates.R
# c = significance cutoff                                                       This sets the cutoff for inclusion in the significant columns
#     #CategoricalCombinePairwiseResults.R; AddAlternatesToCombination.R
#   Permulations:
# x = relaxationValue <int, 0-1>                                                This is a value between 0 and 1, the percentage off of exact a permulation is allowed to be. Increases runspeed, but reduces permulation accuracy. Recommend 0.1 if more than 3 categories
#     #CategoricalPermulationsGeneration.R
# b = b = Metacombination <T or F>                                              This is used to set if the script is targeting metacombined permulations. Called "metacombination". Used for parrallelization. 
#     #CategoricalPermulationsIntermediateCorrelation.R; CategoricalPermualtionsCombineAndCalculatePValue.R
# d = CalculatePValue <T or F>                                                  This value determines whether or not to calculate the pValues.
#     #CategoricalPermualtionsCombineAndCalculatePValue.R; ContinousPermulations.R
# v = v = onlyCalulateValue <T or F>                                            This value is used to run only the p value calculation, and not the combination. 
#     #CategoricalPermualtionsCombineAndCalculatePValue.R

# - Overrides - 
# y = speciesFilterText                                                         OPTIONAL OVERRIDE: This can be used to manually specify a species filer; leave blank for automatic
#     #RunRERAndCorrelation.R;
# z = phenotypeTreeFilename.txt or .rds                                         OPTIONAL OVERRIDE: This can be used to manually override the phenotype tree being used. For continuous analyses, this is the location of the trait vector.
#     #RunRERAndCorrelation.R;
# y = "correlationFileOverride.rds"                                             This can be used to manually set the correlation file. Used primarily to target pairwise comparisons of categorical phenotypes.
#     #EnrichmentAnalysisNew.R; RunAlternateEnrichments.R
# z = "permulationPvalueFileLocation.rds"                                       This is a manual override to specify the script use a specific Permulation p-value file.
#     #EnrichmentAnalysisNew.R; RunAlternateEnrichments.R 
# y = targetFileInstance                                                        This is used to manually set an instance value of the permulations data to be loaded (ie, instance "Dev"). Effects the filename to be loaded in. Defaults to this script's run instance value. Used in parallelization and debugging. 
#     #CategoricalPermulationsIntermediateCorrelation.R
# y = PermulationPrefix                                                         This is used to add a prefix (eg. FastPruned) to permulations if present in the filename.
#     #CategoricalPermualtionsCombineAndCalculatePValue.R
# --- End arguments list ---
#-------------------------------
}


#This function is used to import and arg, and is shared across all imports. 
readArg = function(index, outputVariable = NULL, 
                   successCode = NULL, failureCode = NULL, emptyMessage = NULL,
                   onlyIfUsed = F, usedIndexes = inUse, reportStatus = F){
  if(onlyIfUsed){
    if(usedIndexes == "extant"){
      RunThisArg = exists(outputVariable, envir = .GlobalEnv, inherits = FALSE)
    }else{
      RunThisArg = index %in% usedIndexes
    }
  }else{RunThisArg = T}
  
  if(RunThisArg){
    readArgSuccess = !all(is.na(cmdArgImport(index, T)))
    
    if(readArgSuccess){
      evaluation = cmdArgImport(index, F)
      evaluation = type.convert(evaluation, as.is = TRUE)
      
      assign(outputVariable, evaluation, envir = parent.frame())
      message(paste(outputVariable), " = ", evaluation)
      if(!is.null(successCode)){
        eval(parse(text = successCode), envir = parent.frame())
      }
    }else{

      if(!is.null(failureCode)){
        eval(parse(text = failureCode), envir = parent.frame())
      }else{ #print an empty message if no failure code is run
        if(!is.null(outputVariable) && exists(outputVariable, envir = parent.frame())){
          defaultValue = get(outputVariable, envir = parent.frame(), inherits = TRUE)
          if(is.null(defaultValue)){defaultValue = "NULL"}
          emptyMessage = paste("No", outputVariable, "provided; using", defaultValue)
        }else{
          emptyMessage = paste0("No ", outputVariable, " provided; using default.")
        }
        message(emptyMessage)
      }
      
    }
    message("") #Add a blank line to provide spacing between args 
  }else{
    readArgSuccess = "NotInUse"
  }
  if(reportStatus){return(readArgSuccess)}
}




#-- Figure out how to select which args to use. --
if(exists("usedArgs")){ #if directly specified, use that list. 
  inUse =  usedArgs
}else{ #otherwise, use a check to see if the output variable already exists. 
  inUse = "extant"
}



#This is a framework used to import an arg, shared across arguments. 
# --- 
#Argument: FRAMEWORK
#Description: 
#   []
#   []
#Scripts: 
#   []
Index = 'index'
OutputVariable = "outVar"
successCode = NULL
failureCode = NULL

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)
#---------------------------------------------------------------------
# --- Argument Import Code   --- 
# --------------------------------------------------------------------


# ---- TREE CREATION ----- 

# --- 
#Argument: mainTreesLocation
#Description: 
#   [m  = mainTreeFilename.txt or .rds]
#   [This sets the location of the mainTrees file]
#Scripts: 
#   [MakeCategoricalTree.R; MakeContinousTree.R; RunRERAndCorrelation.R; RunMergedAnalysis.R; ProcessAlternates.R; CategoricalPermulationsCorrelation.R; CategoricalPermualtionsIntermediateCorrelation.R; ContinousPermulations.R]
Index = 'm'
OutputVariable = "mainTreesLocation"
successCode = NULL
failureCode = NULL

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T, reportStatus = F)

# --- 
#Argument: spreadsheet File
#Description: 
#   [td = spreadSheetFilename.csv]
#   [This sets the spreadsheet to read the phenotype data from.]
#Scripts: 
#   [MakeCategoricalTree.R; MakeContinousTree.R]
Index = 'td'
OutputVariable = "spreadSheetLocation"
successCode = NULL
failureCode = NULL
readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Phenotype Column
#Description: 
#   [tp = "phenotypeColumn"]
#   [This sets the spreadsheet to read the phenotype data from.]
#Scripts: 
#   [MakeCategoricalTree.R; MakeContinousTree.R]
Index = 'tp'
OutputVariable = "phenotypeColumn"
successCode = NULL
failureCode = 'stop("ISSUE: SPECIFY FILE PHENOTYPE DATA COLUMN. Aborting.")'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)



# --- 
#Argument: Tip Column
#Description: 
#   [tt = "tipColumn"]
#   [This sets the column with the tip names as they appear in the mainTrees file.]
#Scripts: 
#   [MakeCategoricalTree.R; MakeContinousTree.R]
Index = 'tt'
OutputVariable = "tipCol"
successCode = NULL
failureCode = NULL

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)


# --- 
#Argument: Name Column
#Description: 
#   [tn = "nameColumn"]
#   [This sets the column with common names of the species. ]
#Scripts: 
#   [MakeCategoricalTree.R; MakeContinousTree.R]
Index = 'tn'
OutputVariable = "nameColumn"
successCode = NULL
failureCode = 'message("No name column specified, using in-built list of common name columns.")'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Category Set 
#Description: 
#   [tc = <c("nameOfCategory1","nameOfCategory2")> ]
#   [This sets a list of category names to be used in the analysis. Species with values not in this list will not be included in the analysis. Include ALL categories to be included (even those merged later).  
#Scripts: 
#   [MakeCategoricalTree.R]
Index = 'tc'
OutputVariable = "categoryList"
successCode = NULL
failureCode = 'message("Using all available categories in phenotype column.")'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Screen Column
#Description: 
#   [ts = "screenColumn"]
#   [This is a column which must have a value of 1 for the species to be included.]
#Scripts: 
#   [MakeCategoricalTree.R; MakeContinousTree.R]
Index = 'ts'
OutputVariable = "screenColumn"
successCode = "useScreen = T"
failureCode = 'message("No screen column specified, using all species with valid categories.")'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Substitution list 
#Description: 
#   [tr = list(c("replace1", "with1"),c("replace2, with2"))]
#   [This combined multiple entries to a new value. The new value can be an existing category, or a new label.]
#Scripts: 
#   [MakeCategoricalTree.R]
Index = 'tr'
OutputVariable = "substitutions"
successCode = NULL
failureCode = NULL
isLogical = F

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Merge list 
#Description: 
#   [tm = list(c("phenotype1", "intophen1"), c("2", "i2"))]
#   [This causes combination phenotypes to be merged into the second phenotype, but does not replace standalone phenotypes]
#Scripts: 
#   [MakeCategoricalTree.R]
Index = 'tm'
OutputVariable = "mergeOnlys"
successCode = NULL
failureCode = 'message("No Merge-Mixed-Clasifications Provided; This is normal.")'
isLogical = F

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Ancestral Trait
#Description: 
#   [ta = "ancestralTrait" ]
#   [This can be used to set all non-terminal branches to this category. Use be one of the categories in the list.]
#Scripts: 
#   [MakeCategoricalTree.R]
Index = 'ta'
OutputVariable = "ancestralTrait"
successCode = NULL
failureCode = NULL
isLogical = F

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Color pallete
#Description: 
#   [tb = c( "colors", "in ","the", "palette")  ]
#   [This sets the color palette used in the tree graphs. ]
#Scripts: 
#   [MakeCategoricalTree.R]
Index = 'tb'
OutputVariable = "colorSet"
successCode = NULL
failureCode = NULL
isLogical = F

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Generate alternates
#Description: 
#   [tl = <T or F>]
#   [ This sets if alternate species sets should be created. ]
#Scripts: 
#   [MakeCategoricalTree.R]
Index = 'tl'
OutputVariable = "generateAlternates"
successCode = "numberOfAlternates = 100"
failureCode = NULL
isLogical = T

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Number of Alternates
#Description: 
#   [ta = <number of alternates>   ]
#   [This sets the number of alternate species sets to generate, if alternate species sets are generated. ]
#Scripts: 
#   [MakeCategoricalTree.R]
Index = 'ta'
OutputVariable = "numberOfAlternates"
successCode = NULL
failureCode = NULL
readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: use Liam Inference
#Description: 
#   [tz = <T or F>   ]
#   [ DO NOT USE UNLESS SPECIFIED -- This determines if the liam infrence nodes should be added ]
#Scripts: 
#   [MakeCategoricalTree.R]
Index = 'tz'
OutputVariable = "useLiam"
successCode = 'message("Using Liam inference -- DO NOT DO THIS UNLESS USING A SPECIFIC LIAM PHENOTYPE")'
failureCode = NULL
readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Model Type
#Description: 
#   [q = <ER or SYM or ARD>  ]
#   [This sets the model type used to estimate ancestral branches. ]
#Scripts: 
#   [MakeCategoricalTree.R; CategoricalPermulationsGeneration]
Index = 'q'
OutputVariable = "modelType"
successCode = NULL
failureCode = NULL
readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: ManualPruningSpecies
#Description: 
#   [pd = "c('prunedtip1', 'prunedtip2')" ]
#   [This allows you to manually specify additional branches to be pruned]
#Scripts: 
#   [MakeCategoricalTree.R; MakeContinousTree.R]
Index = 'pm'
OutputVariable = "manualPruningSpecies"
successCode = "usingPruning = T; pruningPrefrenceColumn = NULL; manualPruningProtections = NULL"
failureCode = 'message("No Manual Pruning.")'
readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Pruning cutoff
#Description: 
#   [pm = <minimum branch length>   ]
#   [This sets the minimum branch length for terminal branches in the master tree. Branches shorter than this will be removed. ]
#Scripts: 
#   [MakeCategoricalTree.R; MakeContinousTree.R]
Index = 'pc'
OutputVariable = "pruningCutoff"
successCode = "usingPruning = T; usingAutoPruning = T; pruningPrefrenceColumn = NULL; manualPruningProtections = NULL"
failureCode = 'message("No Automatic Pruing.")'
readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: PruningPrefrenceColumn
#Description: 
#   [pk = "pruningPrefrenceColumn"; MakeContinousTree.R]
#   [This sets a column, where if the value is 1, the tip will be preferentially kept. If the value is TRUE, the tip will never be pruned.]
#Scripts: 
#   [MakeCategoricalTree.R]
Index = 'pk'
OutputVariable = "pruningPrefrenceColumn"
successCode = NULL
failureCode = NULL
readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: ManualPruningProtections
#Description: 
#   [ps = "c('unprunedtip1', 'unprunedtip2')" ]
#   [This allows you to add a list of specific tips to not be dropped during pruning. Must use the tip name, not common name. ]
#Scripts: 
#   [MakeCategoricalTree.R; MakeContinousTree.R]
Index = 'ps'
OutputVariable = "manualPruningProtections"
successCode = NULL
failureCode = NULL
readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)


# - MakeContinousPhenotypeTree.R

# --- 
#Argument: Continuous Metric 
#Description: 
#   [k  = < "Diff" or "mean" or "last" >  ]
#   [This is used for continuous traits, to determine if the metic should be the difference between the nodes (diff), the mean(mean of the two nodes), or last(the downstream value). Note that Mean and Last are not phylogenetically independent. ]
#Scripts: 
#   [MakeContinousTree.R; RunRERAndCorrelation.R; ProcessAlternates.R]
Index = 'k'
OutputVariable = "continousMetric"
successCode = ' validMetrics= c("diff", "mean", "last") 
  if(!any(continousMetric == validMetrics)){
    stop("Contious metric is not a valid option. Use diff, mean, or last")
  }'
failureCode = NULL

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)



# ---- RER ANALYSIS ----- 


# --- 
#Argument: Phenotype Style
#Description: 
#   [a = < ["b" or "binary"] or ["c" or "continuous"] or ["g" or "categorical"]>]
#   [This prefix is used to set the Analysis Style type of phenotype being supplied; it will detect automatically if not specified.]
#Scripts: 
#   [RunRERAndCorrelation.R; ProcessAlternates.R]
Index = 'a'
OutputVariable = "phenotypeStyle"
successCode = 'if(phenotypeStyle == "b" | phenotypeStyle == "B" | phenotypeStyle == "binary"){phenotypeStyle = "Binary"}
  if(phenotypeStyle == "c" | phenotypeStyle == "C" | phenotypeStyle == "continuous"){phenotypeStyle = "Continuous"}
  if(phenotypeStyle == "g" | phenotypeStyle == "G" | phenotypeStyle == "categorical"){phenotypeStyle = "Categorical"}
  message(paste("phenotypeStyle =", phenotypeStyle))'
failureCode = 'message("No phenotype style specified, attempting to detect automatically.")
#Check for categorical 
phenotypeTreeFilename = paste(outputFolderName, filePrefix, "Categorical", "Tree.rds", sep="")
if(file.exists(phenotypeTreeFilename)){message("Found information for Categorical phenotype, using Categorical phenotype"); phenotypeStyle = "Categorical"
}else{
  phenotypeTreeFilename = paste(outputFolderName, filePrefix, "Continuous", "PhenotypeVector.rds", sep="")
  if(file.exists(phenotypeTreeFilename)){message("Found information for continuous phenotype, using continuous phenotype"); phenotypeStyle = "Continuous"
  }else{
    phenotypeTreeFilename = paste(outputFolderName, filePrefix, "Binary", "Tree.rds", sep="")
    if(file.exists(phenotypeTreeFilename)){message("Found information for binary phenotype, using binary phenotype"); phenotypeStyle = "binary"
    }else{
      message("No phenotype specified and no correct tree found.")
      stop()
    }
  } 
}
message("PhenotypeStyle not specified, using continuous")'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: min.sp Value
#Description: 
#   [l = <min.sp value> ]
#   [This sets the min.sp value to be used in the correlation.]
#Scripts: 
#   [RunRERAndCorrelation.R; RunMergedAnalysis.R; ProcessAlternates.R]
Index = 'l'
OutputVariable = "minSpValue"
successCode = NULL
failureCode = NULL

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Species Filter
#Description: 
#   [y = speciesFilterText ]
#   [OPTIONAL OVERRIDE: This can be used to manually specify a species filer; leave blank for automatic ]
#Scripts: 
#   [RunRERAndCorrelation.R]
Index = 'y'
OutputVariable = "speciesFilter"
successCode = NULL
failureCode = 'speciesFilterFileName = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #Make the name of the location a pre-made filter would have to test for it
  if (file.exists(paste(speciesFilterFileName))){speciesFilter = readRDS(speciesFilterFileName); message("Pre-made filter found, using pre-made filter.")
  }else{message("No speciesFilter arg, using all species")}'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Phenotype tree location
#Description: 
#   [z = phenotypeTreeFilename.txt or .rds]
#   [OPTIONAL OVERRIDE: This can be used to manually override the phenotype tree being used. For continuous analyses, this is the location of the trait vector.]
#Scripts: 
#   [RunRERAndCorrelation.R]
Index = 'z'
OutputVariable = "phenotypeTreeLocation"
successCode = NULL
failureCode = 'if(phenotypeStyle == "Continuous"){ phenotypeTreeFilename = paste(outputFolderName, filePrefix, phenotypeStyle, "PhenotypeVector.rds", sep="")
  }else{ phenotypeTreeFilename = paste(outputFolderName, filePrefix, phenotypeStyle, "Tree.rds", sep="")}
  if(file.exists(paste(phenotypeTreeFilename))){                              #if so, use it                
    phenotypeTreeLocation = phenotypeTreeFilename                     
    message("Pre-made Phenotype tree found, using pre-made tree.")
  }else{
    stop("THIS IS AN ISSUE MESSAGE; NO PHENOTYPE TREE/VECTOR FOUND. ENSURE PHENOTYPE TREE/VECTOR FOR SELECTED STYLE AVAILABLE AND PREFIX IS CORRECT.")
  }'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: toMerge
#Description: 
#   [e = c("phenotypes", "to", "merge)]
#   [This sets which phenotype should be combined into a single phenotype in a merged analysis. ]
#Scripts: 
#   [RunMergedAnalysis.R; ProcessAlternates.R]
Index = 'e'
OutputVariable = "toMerge"
successCode = "runMerge = T"
failureCode = 'message("No phenotypes set to merge. Script will effectively do nothing; if you are not running a merged phenotype, this is fine.")
runMerge = F'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: mergeValue
#Description: 
#   [t = c("nameOfMergedPhenotype")]
#   [This sets the name of the now merged phenotype]
#Scripts: 
#   [RunMergedAnalysis.R; ProcessAlternates.R]
Index = 't'
OutputVariable = "mergeValue"
successCode = NULL
failureCode = 'if(!runMerge){stop("ERROR: Attempting to run merge, but the merged phenotype is unnamed. Exiting.")}'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Gmt file location
#Description: 
#   [ g = c("gene", "sets", "file", "paths")]
#   [This sets the gene sets (file names) that the script should be run on, or "gene" for use on genes. ]
#Scripts: 
#   [RunEnrichmentAnalysis.R; CategoricalCombinePairwiseResults.R; RunAlternateEnrichments.R; CombineAlternates.R; AddAlternatesToCombination.R]
Index = 'g'
OutputVariable = "geneSets"
successCode = 'if(is.null(geneSets) | any(geneSets %in% validGeneArguments)){
    usingGene = T
    message("Null geneset specified or gene specified, performing gene correlation.")
  }else{
    usingGene = F
  }
  
  if(any(!geneSets %in% validGeneArguments)){
    usingGo = T
  }else{
    usingGo = F
  }
'
failureCode = 'message("No geneset specified, performing gene correlation."); usingGo = F; usingGene = T'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Permulation use
#Description: 
#   [p = permulations use <T or F or C or A>]
#   [This sets if the code should use permulated or unpermulated values. If set to "C" or "A", will use permulations for categorical values, which are stored in the main file. ]
#Scripts: 
#   [RunEnrichmentAnalysis.R; CategoricalCombinePairwiseResults.R; RunAlternateEnrichments.R]
Index = 'p'
OutputVariable = "usePermulations"
successCode = '  if(usePermulations %in% c("C","A")){
    usePermulations = TRUE
    useCategoricalPermulations = TRUE
    if(usePermulations == "A"){
      addPermulations = T
    }
  }'
failureCode = NULL

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Import subdirectory
#Description: 
#   [s = c("subdirectory", "Names")]
#   [This is used to specify a subdirectory for the analysis to be run in. Primarily used for the components of categorical results. ]
#Scripts: 
#   [RunEnrichmentAnalysis.R; CategoricalCombinePairwiseResults.R; RunAlternateEnrichments.R]
Index = 's'
OutputVariable = "subdirectoryValueList"
successCode = 'useSubdirectory = TRUE; if(length(subdirectoryValueList ==1)){outputFolderName = paste(outputFolderName, subdirectoryValueList[1], "/", sep=""); subdirectoryValue = subdirectoryValueList[1]}'
failureCode = 'message("No subdirectory specified.")
  folderNames = basename(list.dirs(outputFolderName, full.names = TRUE, recursive = FALSE))
  if(length(folderNames) >0 ){
    useSubdirectory = TRUE
    subdirectoryValueList = folderNames
    subdirectoryValueList = subdirectoryValueList[!(subdirectoryValueList == "Alternates")]
    message("Folders not specified, using all folders but Alternates.")
  }else{
    message("Using main directory")
  }'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Import permulation p-value Override 
#Description: 
#   [z = "permulationPvalueFileLocation.rds"]
#   [This is a manual override to specify the script use a specific Permulation p-value file.]
#Scripts: 
#   [RunEnrichmentAnalysis.R; RunAlternateEnrichments.R]
Index = 'z'
OutputVariable = "permulationPValOverride"
successCode = 'usePermulationPValOverride = TRUE; message("Using Manually specified permulation p-value file.")'
failureCode = 'if(!useCategoricalPermulations){
      message("No Permulation pValue override, using standard CombinedBinaryPermulationsPValue.rds.")
    }'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Import correlation file override
#Description: 
#   [y = "correlationFileOverride.rds"]
#   [This can be used to manually set the correlation file. Used primarily to target pairwise comparisons of categorical phenotypes.]
#Scripts: 
#   [RunEnrichmentAnalysis.R; RunAlternateEnrichments.R]
Index = 'y'
OutputVariable = "correlationOverride"
successCode = 'useCorrelationOverride = TRUE; message("Using Manually specified correlation file.")'
failureCode = NULL

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: writeData
#Description: 
#   [w = writeData  <T or F>]
#   [This sets if the data should be saved, potentially overwriting old combinations. ]
#Scripts: 
#   [CategoricalCombinePairwiseResults.R; AddAlternatesToCombination.R]
Index = 'w'
OutputVariable = "writeData"
successCode = NULL
failureCode = NULL

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Significance Cutoff
#Description: 
#   [c = significance cutoff]
#   [This sets the cutoff for inclusion in the significant columns ]
#Scripts: 
#   [CategoricalCombinePairwiseResults.R; AddAlternatesToCombination.R]
Index = 'c'
OutputVariable = "significanceCutoff"
successCode = NULL
failureCode = NULL

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: UsedAlternates
#Description: 
#   [i = run instance  ]
#   [This sets which alternate(s) to run on for parallelization. if left blank, will run on all alternates. Set to "all" to allow a local script to process all alternates.]
#Scripts: 
#   [ProcessAlternates.R; RunAlternateEnrichments.R]
Index = 'i'
OutputVariable = "usedAlternates"
successCode = 'validLocalRunAllAlternates = c("local", "sure", "all")
    if(tolower(usedAlternates) %in% tolower(validLocalRunAllAlternates)){
    message("Run all alternates on local machine CONFIRMED.")
    message("Caution: May have long run times.")
    usedAlternates = NULL
    alternatesSpecified = F
  }else{
    alternatesSpecified = T
  }'
failureCode = 'message("Instance not specified, running on all alternates.")
  if(!clusterRun){
    stop("Attempting to process all alternates on a local machine! This may have run times in excess of 24 hours. Aborting. \n If you wish to avoid this abort, use i=all")
  }'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Alternate cuttof ratio 
#Description: 
#   [o = alternateRatio <decimal>]
#   [This sets the percentage of alternates which must be significant for a result to be considered robust to alternate species selection.]
#Scripts: 
#   [AddAlternatesToCombination.R]
Index = 'o'
OutputVariable = "alternateCutoffRatio"
successCode = NULL
failureCode = NULL

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)



# --- Permulations ---




# --- 
#Argument: Number of permulations
#Description: 
#   [n = numberOfPermulations]
#   [This is the number of permulations to run in the script]
#Scripts: 
#   [CategoricalPermulationsGeneration.R; ContinousPermulations.R]
Index = 'n'
OutputVariable = "permulationAmount"
successCode = NULL
failureCode = NULL

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Run instance value
#Description: 
#   [j = runInstanceValuePerms]
#   [This is used to generate unique filenames for each instance of the script. Used in parrallelization. Set as j to avoid overlap with alternates.]
#Scripts: 
#   [CategoricalPermulationsGeneration.R; CategoricalPermulationsIntermediateCorrelation.R; CategoricalPermulationsCombineAndCalculatePValue.R; ContinousPermulations.R]
Index = 'j'
OutputVariable = "runInstanceValue"
successCode = NULL
failureCode = 'message("This script does not have a run instance value")'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Root Probability
#Description: 
#   [u = <auto or stationary or flat>]
#   [This sets the probability of each state at the root of the tree.]
#Scripts: 
#   [CategoricalPermulationsGeneration.R]
Index = 'u'
OutputVariable = "rootProbability"
successCode = NULL
failureCode = NULL

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Relaxation
#Description: 
#   [x = relaxationValue <int, 0-1>]
#   [This is a value between 0 and 1, the percentage off of exact a permulation is allowed to be. Increases runspeed, but reduces permulation accuracy. Recommend 0.1 if more than 3 categories]
#Scripts: 
#   [CategoricalPermulationsGeneration.R]
Index = 'x'
OutputVariable = "relaxationValue"
successCode = 'useRelaxation = TRUE'
failureCode = 'message("By default, permulations are relaxed to require only 90% matching.")'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Metacombination
#Description: 
#   [b = Metacombination <T or F> ]
#   [This is used to set if the script is targeting metacombined permulations. Called "metacombination". Used for parrallelization. ]
#Scripts: 
#   [CategoricalPermulationsIntermediateCorrelation.R; CategoricalPermulationsCombineAndCalculatePValue.R]
Index = 'b'
OutputVariable = "metacombineValue"
successCode = 'metacombineValue = as.logical(metacombineValue)
  if(is.na(metacombineValue)){
    metacombineValue = FALSE
    message("Metacombination value not interpretable as logical. Did you remember to capitalize? Using FALSE.")
  }'
failureCode = 'message("Metacombination value not specified, using FALSE. If you are not parrallelizing, do not worry about this.")
'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: target File Instance
#Description: 
#   [y = targetFileInstance ]
#   [This is used to manually set an instance value of the permulations data to be loaded (ie, instance "Dev"). Effects the filename to be loaded in. Defaults to this script's run instance value. Used in parallelization and debugging. ]
#Scripts: 
#   [CategoricalPermulationsIntermediateCorrelation.R]
Index = 'y'
OutputVariable = "targetFileInstance"
successCode = NULL
failureCode = '  targetFileInstance = runInstanceValue
  message("No target file Instance specified, defaulting to this scripts run instance.")'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Number of permulations files to combine
#Description: 
#   [nc = <number Of Permulations Files>]
#   [This is the number of permulation files the script will try to combine]
#Scripts: 
#   [CategoricalPermulationsCombineAndCalculatePValue.R]
Index = 'nc'
OutputVariable = "permulationNumberValue"
successCode = NULL
failureCode = '  message("Number of permulations not specified, running on all.")
  intermediateFilenames = list.files(paste0(outputFolderName), pattern = "CategoricalPermulationsIntermediates")
  intermediateFilenames = intermediateFilenames[-grep("Combined", intermediateFilenames)]
  permulationNumberValue = length(intermediateFilenames)'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Start Number
#Description: 
#   [st = <Start number>  ]
#   [This is the permulation number to start at. Used for parrallelization.]
#Scripts: 
#   [CategoricalPermulationsCombineAndCalculatePValue.R]
Index = 'st'
OutputVariable = "startValue"
successCode = NULL
failureCode = NULL

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Calculate P Value
#Description: 
#   [d = CalculatePValue <T or F>]
#   [This value determines whether or not to calculate the pValues.]
#Scripts: 
#   [CategoricalPermulationsCombineAndCalculatePValue.R; ContinousPermulations.R]
Index = 'd'
OutputVariable = "calulateValue"
successCode = '  calulateValue = as.logical(calulateValue)
  if(is.na(calulateValue)){
    calulateValue = FALSE
    message("p-value calculation value not interpretable as logical. Did you remember to capitalize? Using FALSE.")
  }'
failureCode = 'message("p-value calculation value not specified, using TRUE. If you are not parrallelizing, do not worry about this.")'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Only Calculate P Value
#Description: 
#   [v = onlyCalulateValue <T or F> ]
#   [This value is used to run only the p value calculation, and not the combination. ]
#Scripts: 
#   [CategoricalPermulationsCombineAndCalculatePValue.R]
Index = 'v'
OutputVariable = "onlyCalulateValue"
successCode = 'onlyCalulateValue = as.logical(calulateValue)
  if(is.na(calulateValue)){
    onlyCalulateValue = FALSE
    message("Only p-value calculation value not interpretable as logical. Did you remember to capitalize? Using FALSE.")
  }'
failureCode = 'message("Only p-value calculation value not specified, using FALSE. If you are not parrallelizing, do not worry about this.")'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Permulation Prefix
#Description: 
#   [y = "PermulationPrefix" ]
#   [This is used to add a prefix (eg. FastPruned) to permulations if present in the filename.]
#Scripts: 
#   [CategoricalPermulationsCombineAndCalculatePValue.R]
Index = 'y'
OutputVariable = "permulationPrefix"
successCode = NULL
failureCode = NULL

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: Outgroup Species
#Description: 
#   [po = outgroupSpecies ]
#   [This is the outgroup species with the tip name as it appears on mainTrees$masterTree]
#Scripts: 
#   [ContinousPermulations.R]
Index = 'po'
OutputVariable = "outgroupVal"
successCode = NULL
failureCode = 'stop("An outgroup is required, specify outgroup species")'

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)

# --- 
#Argument: FRAMEWORK
#Description: 
#   []
#   []
#Scripts: 
#   []
Index = 'index'
OutputVariable = "outVar"
successCode = NULL
failureCode = NULL

readArg(Index, OutputVariable, successCode, failureCode, onlyIfUsed = T)











