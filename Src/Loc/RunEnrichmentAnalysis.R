.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")

#---- USAGE -----
#Used to generate GO Values from non-permulated RERs.  

# ARUGMENTS: 
#If an argument contains a '(' it is evaluated as code.
# 'r="filePrefix"'                          This is the prefix attached to all files; a required argument. 
#  m = gmtFileLocation.gmt                  This is the location of the main gmt file
# 'a = <annotationListName>'                This is an override for the enrichment annotation list name
#testing args: 
args = c('r=CVHRemake', 'g=F')

#---- Prefix Setup -----
#------------------------
{
  
  #- Import prefix -
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
  
  #- Make Output directory -
  
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
  gmtFileLocation = "Data/enrichmentGmtFile.gmt"
  enrichmentAnnotationListName = "MSigDBPathways"
  
  #--Import arguments--
  #import gmt file location
  if(!is.na(cmdArgImport('m'))){
    gmtFileLocation = cmdArgImport('m')
  }else{
    paste("No gmt location arugment, using Data/enrichmentGmtFile.gmt")                          #Report using default
    message("No gmt location arugment, using Data/enrichmentGmtFile.gmt")
  }
  
  #Annotation List Name 
  if(!is.na(cmdArgImport('a'))){
    enrichmentAnnotationListName = cmdArgImport('a')
  }else{
    paste("No enrichment Annotation List Name argument, using 'MSigDBPathways'")                          #Report using default
    message("No enrichment Annotation List Name argument, using 'MSigDBPathways'")
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
names(annotationsList) = enrichmentAnnotationListName

enrichmentResult = fastwilcoxGMTall(rerStats, annotationsList, outputGeneVals = F)

#save the enrichment output
enrichmentFileName = paste(outputFolderName, filePrefix, "EnrichmentFile.rds", sep= "")
saveRDS(enrichmentResult, enrichmentFileName)



# --- Visualize the enrichment ----
{
  #This is manual only -- run-as-script does not accept a visualize output because no way to output result. 
  #For a script version, use PvQvGoVisualize.R 
  visualize = F
  
  if(visualize){
    library(stringr)
    library(insight)
    enrichmentResult2 = enrichmentResult$enrichmentAnnotationListName
    makeGOTable = function(data, collumn){
      ValueHead = head(data[order(collumn),], n=40)
      ValueHead$num.genes = as.character(ValueHead$num.genes)
      ValueHead$stat = round(ValueHead$stat, digits = 5)
      ValueHead$stat = as.character(ValueHead$stat)
      ValueHead = format_table(ValueHead, pretty_names = F, digits = "scientific5")
      ValueHead
    }
    
    enrichHead = makeGOTable(enrichmentResult2, enrichmentResult2$p.adj)
    enrichHead
    textplot(enrichHead, mar = c(0,0,2,0), cmar = 1.5)
    title(main = paste("Top pathways by non-permulation"))
  }
}
