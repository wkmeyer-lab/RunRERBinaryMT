.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)
library(tools)

# ---- USAGE,README ----
# Command line arguments: 
#All strings must be in quotes. 
#If an argument contains a '(' it is evaluated as code.
# 'm=mainTreeFilename.txt or .rds'
# 'a="annotCollumn"' 
# 'r="filePrefix"'


# ---- Default Arguments ----
#Default main tree file location

mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"

#local computer debug version:
# mainTreesLocation = "../../RemadeTreesAllZoonomiaSpecies.rds"

#Other defaults if not specified
annotCollumn = "ERRORDEFAULT"
filePrefix = "ERRORDEFAULT"


# ---- Command Line Imports ----

args = commandArgs(trailingOnly = TRUE)
paste(args)
message(args)

write.csv(args, file = "Output/MakePathsArgs.csv")

#Main Tree Location
mTreesCommandline = grep("^m=",args, value = TRUE) #get a string based on the identifier

#Process the input
if(length(mTreesCommandline) != 0){                                 #If the string is not empty:
  mainTreesLocationString = substring(mTreesCommandline, 3)         #make a string without the identifier
  if(grepl('(', mTreesCommandline, fixed = TRUE)){                  # if the input is code -- has a '(' in it
    mainTreesLocation = eval(str2lang(mainTreesLocationString))     #convert that string to code, then evaluate that code
  }else{                                                            #if it is not code
    mainTreesLocation = mainTreesLocationString                     #use the string directly 
  } 
  message(mainTreesLocation)                                        #Report the result
}else{                                                              #if it is empty 
  paste("No maintrees arg, using default")                          #Report using default
  message("No maintrees arg, using default")
}

#File Prefix
fPrefixCommandLine = grep("^r=", args, value = TRUE)
if(length(fPrefixCommandLine) != 0){
  filePrefixString = substring(fPrefixCommandLine, 3)
  if(grepl('(', fPrefixCommandLine, fixed = TRUE)){
    filePrefix = eval(str2lang(filePrefixString))
  }else{
    filePrefix = filePrefixString
  }
  message(filePrefix)
}else{
  stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
}

#Annots Collumn
annotsCommandLine = grep("^a=", args, value = TRUE)
if(length(annotsCommandLine) != 0){
  annotsCommandString = substring(annotsCommandLine, 3)
  if(grepl('(', annotsCommandLine, fixed = TRUE)){
    annotCollumn = eval(str2lang(annotsCommandString))
  }else{
    annotCollumn = annotsCommandString
  }
  message(annotCollumn)
}else{
  stop("THIS IS AN ISSUE MESSAGE; SPECIFY ANNOTATION COLLUMN")
}


# -------- Make Paths Main Code ---------------


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

speciesFilterFilename = paste("Results/", filePrefix, "SpeciesFilter.rds",sep="")
saveRDS(relevantSpeciesNames, file = speciesFilterFilename)

# -- Setup foreground Species -- 

foregroundSpeciesAnnot = relevantSpecies[ relevantSpecies[[annotCollumn]] %in% 1,]

foregroundNames = foregroundSpeciesAnnot$FaName

#Make binary tree output 

binaryForegroundTreeOutput = foreground2Tree(foregroundNames, mainTrees, useSpecies = relevantSpeciesNames)

#Save that output 
binaryTreeFilename = paste("Results/", filePrefix, "BinaryForegroundTree.rds", sep="")
saveRDS(binaryForegroundTreeOutput, file = binaryTreeFilename)


