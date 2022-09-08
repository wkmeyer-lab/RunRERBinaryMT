#Tesing code
#mainTrees = readTrees(paste(find.package('RERconverge'),"/extdata/","subsetMammalGeneTrees.txt",sep=""), max.read = 200)
#marineb=read.tree(paste(rerpath,"/extdata/MarineTreeBinCommonNames_noCGM.txt",sep=""))
#binaryPhenotypeTree = marineb



#Library setup 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)

# ---- USAGE,README ----
# Command line arguments: 
#m=mainTreeFilename.txt
#p=phenotypeTreeFilename
#r=filePrefix
#f=speciesFilterText




# ---- MAIN ----

#Read in main tree file 
mainTrees = readTrees("mainTreeFilename.txt") #This is assuming a filetree outside the project on sol 

#Read in the phenotype tree -- branch length is phenotype value, binary 
binaryPhenotypeTree = read.tree("binaryTreeFilename.txt")

#prefix for output files. Typically the phenotype of interest. 
filePrefix = "test" 

#Put a filter for species here. Default is NULL, which means useSpecies is not applied, and all species are used.
#copy of default code: speciesFilter = NULL
speciesFilter = NULL

# ---- Command Line Imports ----

args = commandArgs(TRUE)

#Main Tree Location
mTreesCommandline = grep("m=",args, value = TRUE)
if(!is.null(mTreesCommandline)){
  mainTrees = mTreesCommandline
}

#phenotype tree location
pTreesCommandline = grep("p=",args, value = TRUE)
if(!is.null(pTreesCommandline)){
  binaryPhentotypeTree = pTreesCommandline
}else{
  paste("THIS IS AN ERROR MESSAGE; SPECIFY PHENOTYPE TREE")
}

#File Prefix
fPrefixCommandLine = grep("r=", args, value = TRUE)
if(!is.null(fPrefixCommandLine)){
  filePrefix = fPrefixCommandLine
}else{
  paste("THIS IS AN ERROR MESSAGE; SPECIFY FILE PREFIX")
}

#speciesFilter
sFilterCommandLine = grep("f=", args, value = TRUE)
if(!is.null(sFilterCommandLine)){
  speciesFilter = sFilterCommandLine
}

# ---- RERs ----

RERFileName = paste("output/", filePrefix, "RERFile.rds", sep= "")



if(!file.exists(paste(RERFileName))){
  RERObject = getAllResiduals(mainTrees, useSpecies = speciesFilter, plot = F)
  saveRDS(RERObject, file = RERFileName)
}else{
  RERObject = readRDS(RERFileName)
}

# ---- PATHS ----



pathsFileName = paste("output/", filePrefix, "PathsFile.rds", sep= "")

if(!file.exists(paste(pathsFileName))){
  pathsObject = tree2Paths(binaryPhenotypeTree, mainTrees, binarize=T, useSpecies = speciesFilter)
  saveRDS(pathsObject, file = pathsFileName)
}else{
  RERObject = readRDS(pathsFileName)
}


# ---- Correlate ----

outputFileName = paste("output/", filePrefix, "CorrelationFile", sep= "")

correl = correlateWithBinaryPhenotype(RERObject, pathsObject, min.sp =35)

#head(correl[order(correl$p.adj),])

write.csv(correl, file= paste(outputFileName, ".csv", sep=""), row.names = T, quote = F)
saveRDS(correl, paste(outputFileName, ".rds", sep=""))
   
   
# ----- Arguments integration testing -----

?commandArgs()


