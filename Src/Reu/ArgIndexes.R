#---------------------------------------------------------------------
# --- Combined Tree Generation Control Args  --- 
# --------------------------------------------------------------------
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
# tz = <T or F>                                           DO NOT USE UNLESS SPECIFIED -- This determines if the liam infrence nodes should be added 
# q  = <ER or SYM or ARD>                                 This sets the model type used to estimate ancestral branches 

# - Pruning - 
# pm = <minimum branch length>                            This sets the minimum branch length for terminal branches in the master tree. Branches shorter than this will be removed. 
# pk = "pruningPrefrenceColumn"                           This sets a column, where if the value is 1, the tip will be preferentially kept. If the value is TRUE, the tip will never be pruned.
# ps = "c('unprunedtip1', 'unprunedtip2')"                This allows you to add a list of specific tips to not be dropped during pruning. Must use the tip name, not common name. 
# pd = "c('prunedtip1', 'prunedtip2')"                    This allows you to manually specify additional branches to be pruned
# --- End arguments list ---
#-------------------------------



#---------------------------------------------------------------------
# --- Combined Analysis Control Args  --- 
# --------------------------------------------------------------------
#-------------------------------
# --- Command arguments list ---
# - Required - 
# r = filePrefix                                                                REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# m = mainTreeFilename.txt or .rds                                              REQUIRED: This sets the location of the maintrees file
#     #RunRERAndCorrelation.R; RunMergedAnalysis.R; ProcessAlternates.R; CategoricalPermulationsGeneration.R
# a = < ["b" or "binary"] or ["c" or "continuous"] or ["g" or "categorical"]>   This prefix is used to set the Analysis Style type of phenotype being supplied; it will detect automatically if not specified.
#     #RunRERAndCorrelation.R; ProcessAlternates.R
# g = c("gene", "sets", "file", "paths")                                        This sets the gene sets (file names) that the script should be run on, or "gene" for use on genes. 
#     #EnrichmentAnalysisNew.R; CategoricalCombinePairwiseResults.R; RunAlternateEnrichments.R; CombineAlternates.R; AddAlternatesToCombination.R

# - Required, Defaults likely acceptable -
# f = forceUpdate <T or F>                                                      This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# k = < "diff" or "mean" or "last" >                                            This is used for continuous traits, to determine if the metic should be the difference between the nodes (diff), the mean(mean of the two nodes), or last(the downstream value). Note that Mean and Last are not phylogenetically independent, and do not have downstream processing. 
#     #RunRERAndCorrelation.R; ProcessAlternates.R
# p = <T or F or C or A>                                                        This sets if the code should use permulated or unpermulated values. If set to "C" or "A", will use permulations for categorical values, which are stored in the main file. 
#     #EnrichmentAnalysisNew.R; CategoricalCombinePairwiseResults.R; RunAlternateEnrichments.R
# s = "subdirectoryName"                                                        This is used to specify a subdirectory for the analysis to be run in. Primarily used for the components of categorical results. 
#     #EnrichmentAnalysisNew.R; CategoricalCombinePairwiseResults.R; RunAlternateEnrichments.R 
# w = writeData  <T or F>                                                       This sets if the data should be saved, potentially overwriting old combinations
#     #CategoricalCombinePairwiseResults.R; AddAlternatesToCombination.R
# i = run instance                                                              This sets which alternate(s) to run on for parallelization. if left blank, will run on all alternates. Set to "all" to allow a local script to process all alternates.
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
# b = <T or F>                                                                  This is used to set if the script is targeting metacombined permulations. Called "metacombination". Used for parrallelization. 
#     #CategoricalPermulationsIntermediateCorrelation.R; CategoricalPermualtionsCombineAndCalculatePValue.R
# d = <T or F>                                                                  This value determines whether or not to calculate the pValues.
#     #CategoricalPermualtionsCombineAndCalculatePValue.R; ContinousPermulations.R
# v = <T or F>                                                                  This value is used to run only the p value calculation, and not the combination. 
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
# y = targetFileInstance              This is used to manually set an instance value of the permulations data to be loaded (ie, instance "Dev"). Effects the filename to be loaded in. Defaults to this script's run instance value. Used in parallelization and debugging. 
#     #CategoricalPermulationsIntermediateCorrelation.R
# y = <PermulationPrefix>                     This is used to add a prefix (eg. FastPruned) to permulations if present in the filename.
#     #CategoricalPermualtionsCombineAndCalculatePValue.R
# --- End arguments list ---
#-------------------------------