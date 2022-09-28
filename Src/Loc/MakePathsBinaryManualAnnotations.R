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
# 't=<"uni"/"bi>"
# 'c=<"ancestral"/"all"/"terminal">'
# 'w=<T/F>"


# ---- Default Arguments ----
#Default main tree file location

mainTreesLocation = "/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/RemadeTreesAllZoonomiaSpecies.rds"

#local computer debug version:
# mainTreesLocation = "data/RemadeTreesAllZoonomiaSpecies.rds"
#args = c("m=data/RemadeTreesAllZoonomiaSpecies.rds", "r=Insectivory", "a=Ins_v_herbs", "t=bi", "c=terminal", "w=F")

#Other defaults if not specified
annotCollumn = "ERRORDEFAULT"
filePrefix = "ERRORDEFAULT"
transisitionValue = "Default"
cladeValue = "Default"
weightValue = FALSE


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

#Transition value
transitionCommandLine = grep("^t=", args, value = TRUE)
if(length(transitionCommandLine) != 0){
  transitionCommandString = substring(transitionCommandLine, 3)
  if(grepl('(', transitionCommandLine, fixed = TRUE)){
    transitionValue = eval(str2lang(transitionCommandString))
  }else{
    transitionValue = transitionCommandString
  }
  message(transitionValue)
}else{
  message("Using default unidirectional transistion")
}

#clade value
cladeCommandLine = grep("^c=", args, value = TRUE)
if(length(cladeCommandLine) != 0){
  cladeCommandString = substring(cladeCommandLine, 3)
  if(grepl('(', cladeCommandLine, fixed = TRUE)){
    cladeValue = eval(str2lang(cladeCommandString))
  }else{
    cladeValue = cladeCommandString
  }
  message(cladeValue)
}else{
  message("Using default all clade")
}

#Weight Value
weightCommandLine = grep("^w=", args, value = TRUE)
if(length(weightCommandLine) != 0){
  weightCommandString = substring(weightCommandLine, 3)
  if(grepl('(', weightCommandLine, fixed = TRUE)){
    weightValue = eval(str2lang(weightCommandString))
  }else{
    weightValue = weightCommandString
  }
  message(weightValue)
}else{
  message("Weight = false")
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


# -- set arguments for foreground2Trees --
f2tInputList = list(foregroundNames, mainTrees, useSpecies = relevantSpeciesNames)
#Transition
if(transitionValue == "bi"){
  f2tInputList[["transition"]]= "bidirectional"
  message("Bidirectional transition")
}else{
  f2tInputList[["transition"]]= "unidirectional"
  message("Unidirectional transition")
}

#clade
if(cladeValue == "ancestral"){
  f2tInputList[["clade"]] = "ancestral"
  message("ancestral clade")
}else if(cladeValue == "terminal"){
  f2tInputList[["clade"]] = "terminal"
  message("terminal clade")
}else{
  f2tInputList[["clade"]] = "all"
  message("all clade")
}

#Weight
if(weightValue ==TRUE || weightValue == 't'){
  f2tInputList[["weighted"]] = TRUE
  message("Weighted = TRUE")
}else{
  f2tInputList[["weighted"]] = FALSE
  message("Weighted = False")
}


#Make binary tree output 
#paste(f2tInputList)
f2tInputList
binaryForegroundTreeOutput = do.call(foreground2Tree, f2tInputList)



#Save that output 
binaryTreeFilename = paste("Results/", filePrefix, "BinaryForegroundTree.rds", sep="")
saveRDS(binaryForegroundTreeOutput, file = binaryTreeFilename)
