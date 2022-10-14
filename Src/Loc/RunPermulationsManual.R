#Library setup 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")

#testing args 
#args = c('r=allInsectivory','n=648360')

#default values: 
#maintree and phylo location:
mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"
mainTrees = readRDS(mainTreesLocation)

#file prefix: 
filePrefix = "allInsectivory"

#species filter
speciesFilter = NULL

#root species 
rootSpeciesValue = "REFERENCE"

#Number of permulations
permulationNumberValue = 100

#Run instance value 
runInstanceValue = NULL

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
}else{
  paste("Number of permulations not specified, using 100")
}

#instance of the script 
if(!is.na(cmdArgImport('i'))){
  runInstanceValue = cmdArgImport('i')
}else{
  paste("This script does not have a run instance value")
}


# --------------------------------- MANUAL PORTION ---------------------
#Setup foreground (manual)

#Pull the list of foreground species from the premade tree
#PremadeTree = readRDS("Results/allInsectivoryBinaryForegroundTree.rds")
#pdf(width = 20, height = 200)
#plotTree(PremadeTree)
#dev.off()

#manual foreground
foregroundString = c("vs_eptFus1", "vs_HLpipPip2",
                     "vs_HLlasBor1",
                     "vs_myoBra1", "vs_HLmyoLuc1",
                     "vs_mypDav1", "vs_HLmyoMyo6",
                     "vs_HLmurAurFea1", 
                     "vs_HLminSch1", "vs_HLminNat1", 
                     "vs_HLtadBra1", 
                     "vs_prePar1", "vs_HLmorBla1",
                     "vs_HLmicHir1",
                     "vs_HLhipGal1", "vs_HLhipArm1",
                     "vs_HLrhiSin1",
                     "vs_HLcraTho1",
                     "vs_HLmanPen2", "vs_HLmanJav2",
                     "vs_HLmunMug1"
                     )
#manual making of sister list
sistersList = list("clade1" = c("vs_eptFus1", "vs_HLpipPip2"), 
                   "clade2" = c("clade1","vs_HLlasBor1"), 
                   "clade3" = c("vs_myoBra1", "vs_HLmyoLuc1"),
                   "clade4" = c("vs_mypDav1", "vs_HLmyoMyo6"),
                   "clade5" = c("clade3","clade4"),
                   "clade6" = c("clade5", "vs_HLmurAurFea1"),
                   "clade7" = c("clade6", "clade2"),
                   "clade8" = c("vs_HLminSch1", "vs_HLminNat1"),
                   "clade9" = c("clade7", "clade8"),
                   "clade10"= c("clade9", "vs_HLtadBra1"),
                   "clade11"= c("vs_prePar1", "vs_HLmorBla1"),
                   "clade12"= c("clade11", "vs_HLmicHir1"),
                   "clade13"= c("clade12", "clade10"),
                   "clade14"= c("vs_HLhipGal1", "vs_HLhipArm1"), 
                   "clade15"= c("clade14", "vs_HLrhiSin1"), 
                   "clade16"= c("clade15", "vs_HLcraTho1"),
                   "clade17"= c("clade16", "clade13"),
                   "clade18"= c("vs_HLmanPen2", "vs_HLmanJav2"),
                   "clade19"= c("clade18", "vs_HLmunMug1")
                    )
# ---------------------------- End Manual Portion ----------------------------

# --- Print a copy of the foreground Clades tree
fgCladeTreeFilename = paste(outputFolderName, filePrefix, "CladesForegroundTreeFile.pdf", sep= "")
pdf(file = fgCladeTreeFilename )
foregroundCladeTree = foreground2TreeClades(foregroundString, sistersList, mainTrees, plotTree = T, )
dev.off()


# ------ Clades Paths -----
cladesPathsFileName = paste(outputFolderName, filePrefix, "CladesPathsFile.rds", sep= "")
if(!file.exists(paste(cladesPathsFileName))){
  pathCladesObject = tree2PathsClades(foregroundCladeTree, mainTrees)
  saveRDS(pathCladesObject, file = cladesPathsFileName)
}else{
  pathCladesObject = readRDS(cladesPathsFileName)
}

# ---- RERS -----
#Also allows import of RERs from the non-permulations version of the script 
RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")
if(!file.exists(paste(RERFileName))){
  RERObject = getAllResiduals(mainTrees, useSpecies = speciesFilter, plot = F)
  saveRDS(RERObject, file = RERFileName)
}else{
  RERObject = readRDS(RERFileName)
}

# --- Clades Correlation ---
#This correlation uses the Clades version of the path, and thus cannot be imported from the normal RER script. 

cladesCorellationFileName = paste(outputFolderName, filePrefix, "CladesCorrelationFile", sep= "")
if(!file.exists(paste(cladesCorellationFileName, ".rds", sep=""))){
  cladesCorrelation = correlateWithBinaryPhenotype(RERObject, pathCladesObject, min.sp =35)
  write.csv(cladesCorrelation, file= paste(cladesCorellationFileName, ".csv", sep =""), row.names = T, quote = F)
  saveRDS(cladesCorrelation, file= paste(cladesCorellationFileName, ".csv", sep=""))
}else{
  cladesCorrelation = readRDS(paste(cladesCorellationFileName, ".csv", sep=""))
}


# ----get PermsBinary step ------

#get the root species value from definition above
rootSpecies = rootSpeciesValue

#get number of permulations from definition above
permulationNumber = permulationNumberValue

#master tree
masterTree = mainTrees$masterTree

#set permulation output filename


permulationsCCVersion = getPermsBinary(permulationNumber, foregroundString, sistersList, rootSpecies, RERObject, mainTrees, mastertree =  masterTree, permmode = "cc")

#save the permulations 
permualationsDataFileName = paste(outputFolderName, filePrefix, "PermulationsData", runInstanceValue, ".rds", sep= "")
saveRDS()


permulationPValues = permpvalcor(cladesCorrelation, permulationsCCVersion)

#save the permaulations p values
permulationPValueFileName = paste(outputFolderName, filePrefix, "PermulationsPValue.rds", sep= "")
saveRDS(permulationPValues, file = permulationPValueFileName)


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








#}else{
#  trimmedTree = mainTrees
#}