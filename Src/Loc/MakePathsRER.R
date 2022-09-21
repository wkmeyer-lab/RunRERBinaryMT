.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)
library(tools)


## INSERT MAIN TREE FILE LOCATION HERE -- UPDATE THIS
mainTreesLocation = "../../RemadeTreesAllZoonomiaSpecies.rds"
#

## CHOOSE ANNOTATION COLLUMN NAME ## 
annotCollumn = "Ins_v_herbs"

##Choose file prefix

filePrefix = "test"

# -- Import a tree --
if(file_ext(mainTreesLocation) == "rds"){
  mainTrees = readRDS(mainTreesLocation)
}else{
  mainTrees = readTrees(mainTreesLocation) 
}

#make a Results directory if one doesn't exist
if(!dir.exists("Results")){
  dir.create("Results")
}

# -- Import manual annotations CSv -- 

manualAnnots = read.csv("Data/manualAnnotationsSheet1.csv")

# --- subset the manual annots to only those with data in the collumn
relevantSpecies = manualAnnots[ manualAnnots[[annotCollumn]] %in% c(0,1),]

relevantSpeciesNames = relevantSpecies$FaName

# -- output a species filter for this fileprefix 

speciesFilterFilename = paste("Results/", filePrefix, "SpeciesFilter.rds")
saveRDS(relevantSpeciesNames, file = speciesFilterFilename)

# -- Setup foreground Species -- 

foregroundSpeciesAnnot = relevantSpecies[ relevantSpecies[[annotCollumn]] %in% 1,]

foregroundNames = foregroundSpeciesAnnot$FaName

#Make binary tree output 

binaryForegroundTreeOutput = foreground2Tree(foregroundNames, mainTrees, useSpecies = relevantSpeciesNames)

#Save that output 
binaryTreeFilename = paste("Results/", filePrefix, "BinaryForegroundTree.rds")
saveRDS(binaryForegroundTreeOutput, file = binaryTreeFilename)


