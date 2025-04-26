clusterRun = F
clusterRun = T
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations
library(seqinr)
library(RERconverge)
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/paths2Tree.R")
source("Src/Reu/customSeqinrFunctions.R")
source("Src/Reu/makePhenMasterTree.R")

# -- Command arguments list
# r = filePrefix                                                               This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                                                                 This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# m = mainTreeFilename.txt or .rds                                             This sets the location of the maintrees file
# g = geneName                                                                 This sets the gene for hyphy to be run on 
# p = phenotypeTreeFilename.txt or .rds                                        This can be used to manually override the phenotype tree being used. For continuous analyses, this is the location of the trait vector.
# a = fastafileLocation                                                        This is the location of the fasta alignment
# f = foregroundBranchIndicator                                                This converts a specific path value to "FOREGROUND" in the tree output


#Argument sets
#geneName = "EHHADH"; fileprefix = "CategoricalInsVertivoreTree"; useManualTree = F; fastaLocation = "Results/ENST00000231887.EHHADH.filt.fa"; mainTreesLocation = 'data/zoonomiaAllMammalsTrees.rds'; foregroundCategory = "1"; phenotypeTreeLocation = "Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds"
args = c("g=EHHADH", "r=CategoricalInsVertivoreTree", "a=Results/ENST00000231887.EHHADH.filt.fa")
args = c("g=SDS", "r=CategoricalInsVertivoreTree", "a=Results/ENST00000231887.EHHADH.filt.fa")
args = c("g=EHHADH", "r=CategoricalInsVertivoreTree", "a=Results/ENST00000231887.EHHADHTEST.filt.fa")


# --- Standard start-up code ---
if(clusterRun)args = commandArgs(trailingOnly = TRUE)
{  # Bracket used for collapsing purposes
  #File Prefix
  if(!is.na(cmdArgImport('r'))){
    filePrefix = cmdArgImport('r')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
  }
  
  #  Output Directory 
  if(!dir.exists("Output")){                                      #Make output directory if it does not exist
    dir.create("Output")
  }
  outputFolderNameNoSlash = paste("Output/",filePrefix, sep = "") #Set the prefix sub directory
  if(!dir.exists(outputFolderNameNoSlash)){                       #create that directory if it does not exist
    dir.create(outputFolderNameNoSlash)
  }
  outputFolderName = paste("Output/",filePrefix,"/", sep = "")
  
  #  Force update argument
  forceUpdate = FALSE
  if(!is.na(cmdArgImport('v'))){                                 #Import if update being forced with argument 
    forceUpdate = cmdArgImport('v')
    forceUpdate = as.logical(forceUpdate)
  }else{
    message("Force update not specified, not forcing update")
  }
}

#Make hyphySpecific output folder
outputHyphyFolderNameNoSlash = paste0(outputFolderName, "Hyphy")
if(!dir.exists(outputHyphyFolderNameNoSlash)){                       #create that directory if it does not exist
  dir.create(outputHyphyFolderNameNoSlash)
}
outputHyphyFolderName = paste0(outputHyphyFolderNameNoSlash, "/")


# --- Argument Imports ---
# Defaults
geneName = NULL
fileprefix = NULL
fastaLocation = NULL

if(!clusterRun){mainTreesLocation = '../RunRER/Data/zoonomiaAllMammalsTrees.rds'  #This targets the runRERBInary version of the data by default, regardless of if being run in a separate project or not. 
}
if(clusterRun){mainTreesLocation = '../RunRERBinaryMT/Data/zoonomiaAllMammalsTrees.rds'  #This targets the runRERBInary version of the data by default, regardless of if being run in a separate project or not. 
}

useManualTree = F
phenotypeTreeLocation = NULL
foregroundCategory = "dshakgldskagkjshadgkhalgh" #this is a random string which will likely never be found in the wild; used so that an actual NULL doesn't cause problems

{ # Bracket used for collapsing purposes
  
  #MainTrees Location
  if(!is.na(cmdArgImport('m'))){
    mainTreesLocation = cmdArgImport('m')
  }else{
    message("No maintrees arg, using default")
  }
  
  #phenotype tree location
  if(!is.na(cmdArgImport('p'))){
    phenotypeTreeLocation = cmdArgImport('p')
    useManualTree = T
  }else{                                                                        #See if a pre-made tree for this prefix and style exists 
    message("No manual tree specified, using tree in output folder.")
  }
  
  #Target gene
  if(!is.na(cmdArgImport('g'))){
    geneName = cmdArgImport('g')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY TARGET GENE")
  }
  
  #Fasta File
  if(!is.na(cmdArgImport('a'))){
    fastaLocation = cmdArgImport('a')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY FASTA FILE")
  }
  
  #Foreground Replacement
  if(!is.na(cmdArgImport('f'))){
    foregroundCategory = cmdArgImport('f')
  }else{
    message("No foreground replacement selected, not replacing any path values with FOREGROUND. Note: This is standard behavior.")
  }

}



phenMasterTree = makePhenMasterTree(geneName, filePrefix, manualPhenotypeTreeLocation = phenotypeTreeLocation)


# - Read Fasta file - 
message(fastaLocation)
fasta = read.fasta(fastaLocation)

fastaTipHeaders = names(fasta)
fastaTipHeaders = sub("\\t.*", "", fastaTipHeaders)
fastaTipHeaders[fastaTipHeaders == "REFERENCE"] = "vs_hg38"

# trim files to match eachother 
noDataTips = phenMasterTree$tip.label[!originalTipValues %in% fastaTipHeaders]
phenMasterTree = drop.tip(phenMasterTree, noDataTips)

fastaToDrop = which(!fastaTipHeaders %in% originalTipValues)
fasta = fasta[-fastaToDrop]
fastaTipHeaders = fastaTipHeaders[-fastaToDrop]
names(fasta) = fastaTipHeaders

#Write output file
write.fastaToVar(fasta, names = names(fasta), file.out = "fastaLinesVar")
fastaLinesVar = gsub("!", "", fastaLinesVar) #this removes ! found in some alignment files which breaks hyphy
treeOut = write.tree(phenMasterTree)

combinedContent = c(fastaLinesVar, treeOut)

fastaOutputFilename = paste0(outputHyphyFolderName, filePrefix, geneName, "HyphyInputFile.fna")
writeLines(combinedContent, fastaOutputFilename)
writeLines(combinedContent, "Results/TempHyphyInputFile.fna")


