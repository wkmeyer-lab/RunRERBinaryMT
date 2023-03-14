.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")

#---- USAGE -----
#Used to generate various plots for data analysis. 
#Makes: 
#P value table; Q value table; Histograms of p-values for positive Rho entires, histogram of p-values for negative Rho entires. 

# ARUGMENTS: 
#If an argument contains a '(' it is evaluated as code.
# 'r="filePrefix"'                          This is the prefix attached to all files; a required argument. 
#  m = gmtFileLocation.gmt                  This is the location of the main gmt file 
#testing args: 
args = c('r=CVHRemake', 'g=F')


{
  #---- Prefix Setup -----
  #------------------------
  
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
  
  #------------
  
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
  
  #----------
}

{
  #----Startup----
  #----------------
  
  #--Default Values-- 
  usePermulations = TRUE
  gmtFileLocation = "Data/enrichmentGmtFile.gmt"
  
  
  #--Import arguments--
  #import gmt file location
  if(!is.na(cmdArgImport('m'))){
    gmtFileLocation = cmdArgImport('m')
  }else{
    paste("No arugment, using Data/enrichmentGmtFile.gmt")                          #Report using default
    message("No maintrees arg, using Data/enrichmentGmtFile.gmt")
  }
  
  
  #Import permulation use
  if(!is.na(cmdArgImport('p'))){
    usePermulations = cmdArgImport('p')
    usePermulations = as.logical(usePermulations)
    if(is.na(usePermulations)){
      usePermulations = TRUE
      message("Use Permualtions value not interpretable as logical. Did you remember to capitalize? Using TRUE.")
    }
  }else{
    message("Use Permulations value not specified, using TRUE.")
  }

  
  #--
}

#---- MAIN CODE ----
#-------------------

#Load correlation file
correlationFileLocation = paste(outputFolderName, filePrefix, "CorrelationFile.rds", sep= "")
correlationData = readRDS(correlationFileLocation)                            #Import the correlation data (non-permulated)
rerStats = getStat(correlationData)

#Load the gmt annotations 
gmtAnnotations = read.gmt(gmtFileLocation)
annotationsList = list(gmtAnnotations)
names(annotationsList) = "MSigDBPathways"

enrichmentResult = fastwilcoxGMTall(rerStats, annotationsList, outputGeneVals = F)



# --- Visualize the enrichment ----

enrichmentResult[1]
enrichmentResult = enrichmentResult[order(enrichmentResult$MSigDBPathways$p.adj)]

enrichmentResult2 = enrichmentResult$MSigDBPathways

enrichmentResult2 = enrichmentResult2[order(enrichmentResult2$p.adj),]

enrichmentResult2[,c(1,3,4)]

library(stringr)
library(insight)
data = enrichmentResult2
collumn = enrichmentResult2$p.adj
makeTextPlot = function(data, collumn){
  ValueHead = head(data[order(collumn),], n=40)
  ValueHead$num.genes = as.character(ValueHead$num.genes)
  ValueHead = format_table(ValueHead, pretty_names = F, digits = "scientific5")
  #ValueHead$N = str_sub(ValueHead$N, end = -9)
  ValueHead
}

enrichHead = makeTextPlot(enrichmentResult2, enrichmentResult2$p.adj)
textplot(pValueHead, mar = c(0,0,2,0), cmar = 1.5)
title(main = paste("Top genes by",  targetcolumn))



object.size(enrichmentResult)

library(pathview)


library(help=pathview)

pathview(enrichmentResult)















if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
BiocManager::install("pathview")

