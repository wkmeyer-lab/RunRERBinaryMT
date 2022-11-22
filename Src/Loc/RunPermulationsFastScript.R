#Library setup 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/convertLogiToNumeric.R")
source("Src/Reu/fast_bin_perm.r")

# --- USAGE ---
#Used to run permulations on a given phenotype. Will use pre-exisitng RER values for it's prefix, if they exist. 
#At the moment, the sisters-list must be generated manually. Changing this to automatic is a possible future feature.

#Note that this script includes something which converts any logical vectors (values of entirely NA, typically) from the permulations to numeric vectors. 

# ARUGMENTS: 
#If an argument contains a '(' it is evaluated as code.
# 'r="filePrefix"' This is the prefix attached to all files a required argument. 
# 'm=mainTreeFilename.txt or .rds' This is the location of the maintree file. Accepts .txt or .rds. 
# 'f=speciesFilterText'  This is the text of a species filter. Expects character string. Will use pre-made file, if one exists. 
# 't=rootSpeciesName'     This is the name of the root species, if not using REFERENCE(human)
# 'n=numberOfPermulations' This is the number of permulations to run in the script 
# 'i=runInstanceValue'    This is used to generate unique filenames for each instance of the script. Typically fed in by for loop used to run script in parallel. 
# 'a=<T OR F>'            This stands for "automatic" and if FALSE forces the script to use the manual lists 
# 'e=<integer>'           This specifies the number of internal nodes, overrides an automatic one. 
# 'p=<T or F>'            This specifies if the trees used for permulation should be pruned. This significantly speeds up the permulations but maybe affect results (unclear) 






#testing args 
args = c('r=allInsectivory','n=1','m=Data/RemadeTreesAllZoonomiaSpecies.rds')
args = c('r=carnvHerbs','n=1','m=Data/RemadeTreesAllZoonomiaSpecies.rds')

#Get start time of the script 
timeStart = Sys.time()

#default values: 
#maintree and phylo location:
mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"
#mainTrees = readRDS(mainTreesLocation)

#file prefix: 
filePrefix = NULL

#species filter
speciesFilter = NULL

#root species 
rootSpeciesValue = "REFERENCE"

#Number of permulations
permulationNumberValue = 100

#Run instance value 
runInstanceValue = NULL

#Use automatic lists 
useAutomatic = T

#Manual internal number
useManualInternalNumber = F

#Tree Pruning 
willPruneTree = T

# --- Import prefix ----
args = commandArgs(trailingOnly = TRUE)
paste(args)
message(args)

#File Prefix
if(!is.na(cmdArgImport('r'))){
  filePrefix = cmdArgImport('r')
}else{
  stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
}



#------ Make Output directory -----

#Make output directory if it does not exist
if(!dir.exists("Output")){
  dir.create("Output")
}
#Make a specific subdirectory if it does not exist 
outputFolderNameNoSlash = paste("Output/",filePrefix, sep = "")
#create that directory if it does not exist
if(!dir.exists(outputFolderNameNoSlash)){
  dir.create(outputFolderNameNoSlash)
}
outputFolderName = paste("Output/",filePrefix,"/", sep = "")

#------ Command args import ------

#MainTree location
if(!is.na(cmdArgImport('m'))){
  mainTreesLocation = cmdArgImport('m')
}else{
  paste("No maintrees arg, using default")                          #Report using default
  message("No maintrees arg, using default")
}
mainTrees = readRDS(mainTreesLocation)


#speciesFilter
speciesFilterFileName = paste(outputFolderName, filePrefix, "SpeciesFilter.rds",sep="") #Make the name of the location a pre-made filter would have to test for it

if(!is.na(cmdArgImport('f'))){
  speciesFilter = cmdArgImport('f')
}else if (file.exists(paste(speciesFilterFileName))){                  #See if a pre-made filter for this prefix exists 
  speciesFilter = readRDS(speciesFilterFileName)                       #if so, use it 
  paste("Pre-made filter found, using pre-made filter.")
}else{                                                    
  paste("No speciesFilter arg, using NULL")                           #if not, use no filter
}

#Root species
if(!is.na(cmdArgImport('t'))){
  rootSpeciesValue = cmdArgImport('t')
}else{
  paste("No root species specified, using 'REFERENCE'")
}

#Number of permulations
if(!is.na(cmdArgImport('n'))){
  permulationNumberValue = cmdArgImport('n')
  permulationNumberValue = as.numeric(permulationNumberValue)
}else{
  paste("Number of permulations not specified, using 100")
}

#instance of the script 
if(!is.na(cmdArgImport('i'))){
  runInstanceValue = cmdArgImport('i')
}else{
  paste("This script does not have a run instance value")
}
print(permulationNumberValue)
str(permulationNumberValue)

#use automatic sisterlist 
if(!is.na(cmdArgImport('a'))){
  useAutomatic = cmdArgImport('a')
  useAutomatic = as.logical(useAutomatic)
}else{
  paste("Manual list not being forced, using automatic if available")
}

#number of internal nodes 
if(!is.na(cmdArgImport('e'))){
  useManualInternalNumber = T
  manualInternalNumber = cmdArgImport('e')
}else{
  paste("Manual internal number not specified, using automatic if available")
}

#Should the tree be pruned
if(!is.na(cmdArgImport('p'))){
  willPruneTree = cmdArgImport('p')
  willPruneTree = as.logical(willPruneTree)
}else{
  paste("Tree pruning not specified, pruning tree.")
}
# --------------------------------- BEGIN SCRIPT ---------------------


# ---- RERS -----
#Also allows import of RERs from the non-permulations version of the script 
RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")
if(!file.exists(paste(RERFileName))){
  RERObject = getAllResiduals(mainTrees, useSpecies = speciesFilter, plot = F)
  saveRDS(RERObject, file = RERFileName)
}else{
  RERObject = readRDS(RERFileName)
}

# --- phentotypeVector ---
phenotypeVectorFilename = paste(outputFolderName, filePrefix, "phenotypeVector.rds", sep="")
if(file.exists(phenotypeVectorFilename = paste(outputFolderName, filePrefix, "phenotypeVector.rds", sep=""))){
  phenotypeVector = readRDS(phenotypeVectorFilename)
}else{
  stop("THIS IS AN ISSUE MESSAGE, GENERATE A PHENTOTYPEVECTOR (sistersListGeneration.R)")
}

# --- Number of internal nodes --- 
internalNodeFilename = paste(outputFolderName, filePrefix, "internalNodeNumber.rds", sep="")
if(useManualInternalNumber){
  internalNumber = manualInternalNumber
}else if(file.exists(internalNodeFilename)){
  internalNumber = readRDS(internalNodeFilename)
}else{
  stop("THIS IS AN ISSUE MESSAGE: NO MANUAL INTERNAL NUMBER SPECIFIED, AND NO AUTOMATIC FOUND. EXITING.")
}


# ----get Permulations  step ------

#get the root species value from definition above
rootSpecies = rootSpeciesValue

#get number of permulations from definition above
permulationNumber = permulationNumberValue

#master tree
masterTree = mainTrees$masterTree

#print the current runtime
timeBefore = Sys.time()
runTimeBefore = timeBefore - timeStart
message("Runtime: ", runTimeBefore)


# ----- Using the fast permulations ---- 


#Make a rooted version of the master tree
rootNode = which(masterTree$tip.label %in% rootSpeciesValue)
rootedMasterTree = multi2di(masterTree)
plot(masterTree)
plot(rootedMasterTree)
if(willPruneTree){
  prunedMasterTree = pruneTree(rootedMasterTree, names(phenotypeVector))
  rootedMasterTree = prunedMasterTree
  plot(rootedMasterTree)
}
#The function used for each permulation:
computeCorrelationOnePermulation = function(rootedMasterTree, phenotypeVector, mainTrees, RERObject, min.sp =35, internalNumber){
  permulatedForeground = fastSimBinPhenoVecReport(tree=rootedMasterTree, phenvec=phenotypeVector, internal=internalNumber)                                     #generate a null foreground via permulation
  permulatedTree = foreground2Tree(permulatedForeground, mainTrees, plotTree=F, clade="all", transition="bidirectional", useSpecies=speciesFilter) #generate a tree using that foregound 
  permulatedPaths = tree2Paths(permulatedTree, mainTrees, binarize=T, useSpecies=speciesFilter)                                                    #generate a path from that tree
  permulatedCorrelations = correlateWithBinaryPhenotype(RERObject, permulatedPaths, min.sp=min.sp)                                                 #Use that path to get a coreelation of the null phenotype to genes (this is the outbut of a Get PermsBinary run)
}

#run repeated permulations
correlationList = list()
for(i in 1:permulationNumber){                                                  #Repeat for the number of permulations
  singlePermCorrelation = computeCorrelationOnePermulation(rootedMasterTree, phenotypeVector, mainTrees, RERObject) #run one permulation each time
  correlationList = append(correlationList, list(singlePermCorrelation))        #add it to a growing list of the dataframes outputted from CorrelateWithBinaryPhenotype
  message("Completed permulation: ", i)                                         #report completed the permulation
}

#Convert the fast output into the getPermsBinary output format
convertPermulationFormat = function(permulationCorList, RERObject = RERObject, permulationNumber = permulationNumber){
  permulationCorList
  permPvals = data.frame(matrix(ncol = permulationNumber, nrow = nrow(RERObject)))
  rownames(permPvals) = rownames(RERObject)
  permRhovals = data.frame(matrix(ncol = permulationNumber, nrow = nrow(RERObject)))
  rownames(permRhovals) = rownames(RERObject)
  permStatvals = data.frame(matrix(ncol = permulationNumber, nrow = nrow(RERObject)))
  rownames(permStatvals) = rownames(RERObject)
  for (i in 1:length(permulationCorList)) {
    permPvals[, i] = permulationCorList[[i]]$P
    permRhovals[, i] = permulationCorList[[i]]$Rho
    permStatvals[, i] = sign(permulationCorList[[i]]$Rho) * -log10(permulationCorList[[i]]$P)
  }
  output = vector("list", 3)
  output[[1]] = permPvals
  output[[2]] = permRhovals
  output[[3]] = permStatvals
  names(output) = c("corP", "corRho", "corStat")
  output
}
convertedPermulations = convertPermulationFormat(correlationList)


# ----- end of using fast permulations ----


#Get time spent on permulations
timeAfter = Sys.time()
runTimeAfter = timeAfter - timeStart
runTimeOfPerms = timeAfter - timeBefore
message("Total runtime: ", runTimeAfter)
message("Permulation runtime: ", runTimeOfPerms)


#-save the permulations- 
#Make different filenames based on if the tree is pruned or not
if(willPruneTree){
  permualationsDataFileName = paste(outputFolderName, filePrefix, "PrunedFastPermulationsData", runInstanceValue, ".rds", sep= "")
}else{
  permualationsDataFileName = paste(outputFolderName, filePrefix, "UnprunedFastPermulationsData", runInstanceValue, ".rds", sep= "")
}
saveRDS(convertedPermulations, file = permualationsDataFileName)

#get time spent on saving the file 
timePostSave = Sys.time()
runTimeAfter = timePostSave - timeStart
runTimeOfSaving = timePostSave - timeAfter
message("Total runtime: ", runTimeAfter)
message("Saving runtime: ", runTimeOfSaving)












# ---- DISABLED: Calculate the permulation P values -----
##Disabled to allow for pValue calculation on combined permutation datasets

#Calculate the permulations P values
#permulationPValues = permpvalcor(cladesCorrelation, permulationsCCVersion)

#save the permaulations p values
#permulationPValueFileName = paste(outputFolderName, filePrefix, "PermulationsPValue.rds", sep= "")
#saveRDS(permulationPValues, file = permulationPValueFileName)


# ----- Permulations with enrichments -----
#Currently disabled
#load annotations
#annotations=RERconverge::read.gmt("Data/gmtfile.gmt")
#annotationslist=list(annotations)
#names(annotationslist)="MSigDBpathways"

#permsCCWithEnrichment = getPermsBinary(permulationNumber, foregroundString, sistersList, rootSpecies, RERObject, mainTrees, mastertree =  masterTree, permmode = "cc", calculateenrich = T, annotlist =  annotationslist)




# ----- Limit the number of genes permulations used 
#if(limit = T){
#  #get a list of genes with low correlation p values
#  corellationOrdered = sort(cladesCorrelation)
#  correlationCutoff = correlationOrdered[correlationOrdered < 0.1]
#  #add the master tree to that list 
#  correlationCutoff = append(correlationCutoff, 1, 0)
#  names(correlationCutoff)[1] = masterTree
#  correlationCutoff

length(which(rootedMasterTree$tip.label %in% names(phenotypeVector)))
length(phenotypeVector)

inputTree = masterTree
inputTree = rootedMasterTree
#dev.off(); dev.new(); dev.new(); testplot2 = plotTreeHighlightBranches(inputTree,hlspecies=which(inputTree$edge.length== 3),hlcols="blue", main="Marine mammals trait tree"); edgelabels(cex = 0.7, frame="none", font=2, adj=c(0,-0.2), col="blue"); nodelabels(cex = 0.7, frame="none", font=2, adj=c(-0.2,0.3), col="dark green"); tiplabels(cex = 0.8, frame="none", font=2, adj=c(0.2,0), col="dark red")



data(bird.orders)
plot(bird.orders)
plot(root(bird.orders, 1))
plot(root(bird.orders, 1:5))

tr <- root(bird.orders, 1)
is.rooted(bird.orders) # yes
is.rooted(tr)          # no

is.rooted(multi2di(masterTree))
is.rooted(root(bird.orders, 1, r = TRUE))

plot(multi2di(tr))
plot(root(bird.orders, 1, r = TRUE))
plot(multi2di(masterTree))


"vs_aotNan1" %in% names(phenotypeVector)
missingNames = masterTree$tip.label[(!masterTree$tip.label %in% names(phenotypeVector))]
missingNames
expandedPhenotypeVector = phenotypeVector
for(i in 1:length(missingNames)){
  addSpecies = 0 
  names(addSpecies) = missingNames[i]
  expandedPhenotypeVector = append(expandedPhenotypeVector, addSpecies )
}
expandedPhenotypeVector
masterTree$tip.label %in% names(expandedPhenotypeVector)
permulatedForeground = fastSimBinPhenoVec(tree=rootedMasterTree, phenvec=expandedPhenotypeVector, internal=internalNumber) 

timeBefore = Sys.time()
permulatedForeground = fastSimBinPhenoVec(tree=rootedMasterTree, phenvec=phenotypeVector, internal=internalNumber)                                     #generate a null foreground via permulation
timeAfter = Sys.time()
runTimeOfPerms = timeAfter - timeBefore
message("Permulation runtime: ", runTimeOfPerms)
