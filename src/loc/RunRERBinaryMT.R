#Tesing code
#rerpath = find.package('RERconverge')
#mainTrees = readTrees(paste(find.package('RERconverge'),"/extdata/","subsetMammalGeneTrees.txt",sep=""), max.read = 200)
#marineb=read.tree(paste(rerpath,"/extdata/MarineTreeBinCommonNames_noCGM.txt",sep=""))
#binaryPhenotypeTree = marineb
#data("logAdultWeightcm")
#args = c("m=C:/Users/Michael/AppData/Local/R/win-library/4.2/RERconverge/extdata/subsetMammalGeneTrees.txt","p=C:/Users/Michael/AppData/Local/R/win-library/4.2/RERconverge/extdata/MarineTreeBinCommonNames_noCGM.txt", "r=Command", 'f=names(logAdultWeightcm)')
#args = c("m=C:/Users/Michael/AppData/Local/R/win-library/4.2/RERconverge/extdata/subsetMammalGeneTrees.txt","p=C:/Users/Michael/AppData/Local/R/win-library/4.2/RERconverge/extdata/MarineTreeBinCommonNames_noCGM.txt", 'f=names(logAdultWeightcm)')
#testTreePath = paste(find.package('RERconverge'),"/extdata/","subsetMammalGeneTrees.txt",sep="")
#args = c('m=paste(find.package("RERconverge"),"/extdata/","subsetMammalGeneTrees.txt",sep="")', 'p=paste(find.package("RERconverge"),"/extdata/MarineTreeBinCommonNames_noCGM.txt",sep="")','r="Command"', 'f=names(logAdultWeightcm)')

# sol args:   'm=paste(find.package("RERconverge"),"/extdata/","subsetMammalGeneTrees.txt",sep="")' 'p=paste(find.package("RERconverge"),"/extdata/MarineTreeBinCommonNames_noCGM.txt",sep="")' 'r="Command"' 'f=names(logAdultWeightcm)'

#Library setup 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)

# ---- USAGE,README ----
# Command line arguments: 
#All strings must be in quotes. Arguments are evaluated as code.
#m=mainTreeFilename.txt
#p=phenotypeTreeFilename
#r="filePrefix"
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

sink("debugmessages.txt")

# ---- Command Line Imports ----

args = commandArgs(trailingOnly = TRUE)
paste(args)
message(args)

#Main Tree Location
mTreesCommandline = grep("^m=",args, value = TRUE) #get a string based on the identifier
if(length(mTreesCommandline) != 0){                      #If the string is not empty:
  mainTreesLocationString = substring(mTreesCommandline, 3)    #make a string without the identifier
  mainTreesLocation = eval(str2lang(mainTreesLocationString))  #convert that string to code, then evaluate that code
}else{
  paste("No maintrees arg, using default")
  message("No maintrees arg, using default")
}

#phenotype tree location
pTreesCommandline = grep("^p=",args, value = TRUE)
if(length(pTreesCommandline) != 0){
  binaryPhenotypeTreeLocationString = substring(pTreesCommandline,3)
  binaryPhenotypeTreeLocation = eval(str2lang(binaryPhenotypeTreeLocationString))
}else{
  #paste("THIS IS AN ERROR MESSAGE; SPECIFY PHENOTYPE TREE")
  stop("THIS IS AN ERROR MESSAGE; SPECIFY PHENOTYPE TREE")
}

#File Prefix
fPrefixCommandLine = grep("^r=", args, value = TRUE)
if(length(fPrefixCommandLine) != 0){
  filePrefixString = substring(fPrefixCommandLine, 3)
  filePrefix = eval(str2lang(filePrefixString))
}else{
  stop("THIS IS AN ERROR MESSAGE; SPECIFY FILE PREFIX")
}

#speciesFilter
sFilterCommandLine = grep("^f=", args, value = TRUE)      #get a string based on the identifier
if(length(sFilterCommandLine) != 0){                        #If the string is not empty:
  sFilterCommandString = (substring(sFilterCommandLine, 3))  #get a string without the identifier
  speciesFilter = eval(str2lang(sFilterCommandString))     #convert that string to code, then evaluate that code
}else{
  paste("No speciesFilter arg, using default")
}


# ---- MAIN ----

#Read in main tree file 
mainTrees = readTrees(mainTreesLocation, max.read = 200) 

#Read in the phenotype tree -- branch length is phenotype value, binary 
binaryPhenotypeTree = read.tree(binaryPhenotypeTreeLocation)

#make an output directory if one doesn't exist
if(!dir.exists("output")){
  dir.create("output")
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
  pathsObject = readRDS(pathsFileName)
}


# ---- Correlate ----

outputFileName = paste("output/", filePrefix, "CorrelationFile", sep= "")

correl = correlateWithBinaryPhenotype(RERObject, pathsObject, min.sp =35)

#head(correl[order(correl$p.adj),])

write.csv(correl, file= paste(outputFileName, ".csv", sep=""), row.names = T, quote = F)
saveRDS(correl, paste(outputFileName, ".rds", sep=""))
   
paste("Script successful.")
message("Script successful.")
# ----- Arguments integration testing -----




