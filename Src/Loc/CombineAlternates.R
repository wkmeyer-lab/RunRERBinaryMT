clusterRun = F
clusterRun = T
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations

library(RERconverge)
library(tools)
source("Src/Reu/cmdArgImport.R")
library(data.table)

# -- Command arguments list
# r = filePrefix                                                               This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                                                                 This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# g = geneSet                                                                  This is the gmtset(as shown in the filename) of enrichment runs to be combined. 

#----------------
args = c('r=ComplexDietCentralAnalysis', 's=g', 'g=KeggReactome', 'l=170', 'i=1')


geneSet = NULL
useGO = F

#Gmt set
if(!any(is.na(cmdArgImport('g')))){
  geneSet = cmdArgImport('g')
  useGO = T
}else{
  paste("No gmt  arugment, not combining enrichments")                          #Report using default
  message("No gmt  arugment, not combining enrichments")
}


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

primaryOutputFolderName = outputFolderName
primaryFilePrefix = filePrefix

#---
#load correlations
CorrelationsFilenames = list.files(paste0(outputFolderName, "Alternates/"), pattern = "PairwiseCorrelationFile.rds")
CorrelationsFilenames = CorrelationsFilenames[-grep("Merged", CorrelationsFilenames)]

CorrelationsFilenames = paste0(outputFolderName, "Alternates/", CorrelationsFilenames)

correlationRuns = lapply(CorrelationsFilenames, readRDS)


#----
#Combine correlations
#----


combinedCorrelations <- lapply(names(correlationRuns[[1]]), function(nm) {
  
  # Extract list of data frames for this comparison
  dfs <- lapply(correlationRuns, function(x) x[[nm]])
  
  n <- nrow(dfs[[1]])
  
  data.frame(
    Rho   = I(lapply(seq_len(n), function(i) sapply(dfs, function(df) df$Rho[i]))),
    P     = I(lapply(seq_len(n), function(i) sapply(dfs, function(df) df$P[i]))),
    p.adj = I(lapply(seq_len(n), function(i) sapply(dfs, function(df) df$p.adj[i])))
  )
})
names(combinedCorrelations) <- names(correlationRuns[[1]])




combinedCorrelations <- lapply(combinedCorrelations, function(df) {
  
  df$PNumSignificant<- sapply(df$P, function(x) sum(x < 0.05, na.rm = TRUE))
  
  df$PadjNumSignificant <- sapply(df$p.adj, function(x) sum(x < 0.05, na.rm = TRUE))
  
  df$RhoMaxDiff <- sapply(df$Rho, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  })
  
  df$PMaxDiff <- sapply(df$P, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  })
  
  df$PadjMaxDiff <- sapply(df$p.adj, function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  })
  
  df
})


saveRDS(combinedCorrelations, paste0(primaryOutputFolderName, primaryFilePrefix, "AlternatesCombinedCorrelations.rds"))


#----
#Combine enrichments
#----

if(useGO){
  alternateSets = readRDS(paste0(outputFolderName, filePrefix, "AlternatePruningSpecies.rds"))
  
  
  subdirectories = list.dirs(paste0(outputFolderName, "Alternates"))
  subdirectories = subdirectories[-which(subdirectories == paste0(outputFolderName, "Alternates"))]
  subdirectories = gsub(paste0(outputFolderName, "Alternates/"), "", subdirectories)
  
  pairwiseSets = subdirectories
  
  

  outputFolderName = paste0(outputFolderName, "Alternates/")

  CombineGo = T
  saveCombinedData = T
  
  
  if(CombineGo){
    for(j in 1:length(alternateSets)){
      currentAlternate = j
      filePrefix = paste0("Alternate", j, primaryFilePrefix)
  
        GOResults = list()
        for(i in 1:length(pairwiseSets)){
          currentSet = pairwiseSets[i]
          correlationSubsetName = gsub("-", " - ", currentSet)
          correlationPrefix = paste(substr(strsplit(currentSet, split = "-")[[1]],1,1), collapse = '')
          
          goFilename = paste0(outputFolderName, currentSet, "/", filePrefix, currentSet,"Enrichment-", geneSet, ".rds")
          currentGoData = readRDS(goFilename)[[1]]
          
          
          
          #names(currentGoData) = paste0(correlationPrefix, "-", names(currentGoData))
          
    
          GOResults[[i]] = currentGoData
          names(GOResults)[i] = currentSet
          
        }
        rm(currentGoData)
        

        #rm(GOResults)
        
       
        
        if(saveCombinedData){
          combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResults-", geneSet)
          saveRDS(GOResults, paste0(combinedGODataFilename, ".rds"))
        }
    }
  }
  
  
  #---- 
  #combine alternates
  
  EnrichmentFilenames = list.files(paste0(primaryOutputFolderName, "Alternates/"), pattern = paste0("combinedGOResults-", geneSet, ".rds"))
  
  EnrichmentFilenames = paste0(primaryOutputFolderName, "Alternates/", EnrichmentFilenames)
  
  EnrichmentRuns = lapply(EnrichmentFilenames[1:4], readRDS)
  
  combinedEnrichment <- lapply(seq_along(EnrichmentRuns[[1]]), function(g) {
    
    group_name <- names(EnrichmentRuns)[g]
    
    dfs <- lapply(EnrichmentRuns, function(run) run[[g]])
    
    n <- nrow(dfs[[1]])
    
    data.frame(
      stat = I(lapply(seq_len(n), function(i)
        sapply(dfs, function(df) df$stat[i])
      )),
      
      pval = I(lapply(seq_len(n), function(i)
        sapply(dfs, function(df) df$pval[i])
      )),
      
      p.adj = I(lapply(seq_len(n), function(i)
        sapply(dfs, function(df) df$p.adj[i])
      ))
    )
    
  })
  
  names(combinedEnrichment) <- names(EnrichmentRuns[[1]])
  
  
  combinedEnrichment <- lapply(combinedEnrichment, function(df) {
    
    df$pNumSignificant <- sapply(df$pval, function(x) sum(x < 0.05, na.rm = TRUE))
    df$padjNUmSignificant <- sapply(df$p.adj, function(x) sum(x < 0.05, na.rm = TRUE))
    
    df$statMaxDiff <- sapply(df$stat, function(x)
      if (all(is.na(x))) NA_real_ else diff(range(x, na.rm = TRUE))
    )
    
    df$pMaxDiff <- sapply(df$pval, function(x)
      if (all(is.na(x))) NA_real_ else diff(range(x, na.rm = TRUE))
    )
    
    df$padjMaxDiff <- sapply(df$p.adj, function(x)
      if (all(is.na(x))) NA_real_ else diff(range(x, na.rm = TRUE))
    )
    
    df
  })
  
  saveRDS(combinedEnrichment, paste0(primaryOutputFolderName, primaryFilePrefix, "AlternatesCombinedEnrichments-" , geneSet, ".rds"))
  
}
