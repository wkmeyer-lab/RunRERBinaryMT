#Tesing code
#mainTrees = readTrees(paste(find.package('RERconverge'),"/extdata/","subsetMammalGeneTrees.txt",sep=""), max.read = 200)
#marineb=read.tree(paste(rerpath,"/extdata/MarineTreeBinCommonNames_noCGM.txt",sep=""))
#binaryPhenotypeTree = marineb



#Library setup 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package

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
   
   

