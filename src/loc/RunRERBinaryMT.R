#Tesing code
#rerpath = find.package('RERconverge')
#mainTrees = readTrees(paste(find.package('RERconverge'),"/extdata/","subsetMammalGeneTrees.txt",sep=""), max.read = 200)
#marineb=read.tree(paste(rerpath,"/extdata/MarineTreeBinCommonNames_noCGM.txt",sep=""))
#binaryPhenotypeTree = marineb
#data("logAdultWeightcm")
#args = c("m=C:/Users/Michael/AppData/Local/R/win-library/4.2/RERconverge/extdata/subsetMammalGeneTrees.txt","p=C:/Users/Michael/AppData/Local/R/win-library/4.2/RERconverge/extdata/MarineTreeBinCommonNames_noCGM.txt", "r=Command", 'f=names(logAdultWeightcm)')
#args = c("m=C:/Users/Michael/AppData/Local/R/win-library/4.2/RERconverge/extdata/subsetMammalGeneTrees.txt","p=C:/Users/Michael/AppData/Local/R/win-library/4.2/RERconverge/extdata/MarineTreeBinCommonNames_noCGM.txt", 'f=names(logAdultWeightcm)')
#testTreePath = paste(find.package('RERconverge'),"/extdata/","subsetMammalGeneTrees.txt",sep="")

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


# ---- Default values if no arguments

#default maintree and phylo location:
mainTreesLocation = "filename.txt"
binaryPhenotypeTreeLocation = ""

#prefix for output files. Typically the phenotype of interest. 
filePrefix = "test" 

#Put a filter for species here. Default is NULL, which means useSpecies is not applied, and all species are used.
#copy of default code: speciesFilter = NULL
speciesFilter = NULL

# ---- Command Line Imports ----

args = commandArgs(trailingOnly = TRUE)

#Main Tree Location
mTreesCommandline = grep("m=",args, value = TRUE) #get a string based on the identifier
if(length(mTreesCommandline) != 0){                      #If the string is not empty:
  mainTreesLocation = substring(mTreesCommandline, 3)    #set to a string without the identifier
}else{
  paste("No maintrees arg, using default")
}

#phenotype tree location
pTreesCommandline = grep("p=",args, value = TRUE)
if(length(pTreesCommandline) != 0){
  binaryPhenotypeTreeLocation = substring(pTreesCommandline,3)
}else{
  #paste("THIS IS AN ERROR MESSAGE; SPECIFY PHENOTYPE TREE")
  stop("THIS IS AN ERROR MESSAGE; SPECIFY PHENOTYPE TREE")
}

#File Prefix
fPrefixCommandLine = grep("r=", args, value = TRUE)
if(length(fPrefixCommandLine) != 0){
  filePrefix = substring(fPrefixCommandLine, 3)
}else{
  stop("THIS IS AN ERROR MESSAGE; SPECIFY FILE PREFIX")
}

#speciesFilter
sFilterCommandLine = grep("f=", args, value = TRUE)      #get a string based on the identifier
if(length(sFilterCommandLine) != 0){                        #If the string is not empty:
  sFilterCommandString = (substring(sFilterCommandLine, 3))  #get a string without the identifier
  speciesFilter = eval((str2lang(sFilterCommandString)))     #convert that string to code, then evaluate that code
}else{
  paste("No speciesFilter arg, using default")
}


# ---- MAIN ----

#Read in main tree file 
mainTrees = readTrees(mainTreesLocation) 

#Read in the phenotype tree -- branch length is phenotype value, binary 
binaryPhenotypeTree = read.tree(binaryPhenotypeTreeLocation)





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


