# -- Usage:
#This script combines all of the pairwise comparisons, and the overall comparison, of a categorical RERConverge analysis into a single file. 
#An older version of the script added additional information (comparisons between the pairwise analyses), but these additions were removed due to making the output file absurdly log.
#For the older version of the script which included overlap, delta, and driver analyses, see CategoricalCombinePairwiseResults-OLD.R

# -- Libraries 
source("Src/Reu/loadLibraries.R")
loadLibraries() # Note that this library load includes cluster run detection. 
source("Src/Reu/trimGmtFilenameToSetname.R")
#-------------------------------
# --- Command arguments list ---
# - Required - 
# r = filePrefix                              REQUIRED: This is a prefix used to organize and separate files by analysis run. Always required. 
# g = c("gene", "sets")                       This sets the gene set (set name or file name) that the script should be run on, or "gene" for use on genes. 

# - Required, Defaults likely acceptable -
# f = forceUpdate <T or F>                    This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# s = c("subdirectory", "Names")              This is used to specify a subdirectory for the analysis to be run in. Primarily used for the components of categorical results. 
# p = permulations use <T or F or C or A>     This sets if the code should use permulated or unpermulated values. If set to "A", it will add permulations data, but use the unpermulationed data for significance. It also accepts "C" as "T" for compatibility with upstream functions which use "C" to denote the categorical permulations used by this script.
# w = writeData  <T or F>                     This sets if the data should be saved, potentially overwriting old combinations. 
# - Optional - 
# h = <T>                                     This ARG causes the script to return its ARG list, instead of running the script. 
# c = significance cutoff                     This sets the cutoff for inclusion in the significant columns
# - Overrides - 
# --- End arguments list ---
#-------------------------------
# -- Section to include manually-specified args --
if(!exists("args", envir = environment(), inherits = FALSE)){
  #This is a space where you can manually write in arguments for the script. 
  #The above if statement ensures these lines will not override any existing args if the script is run whole-sale. 
  args = c('r=ComplexDietCentralAnalysisSimplify', 'g=c("gene", "KeggReactome")')
}
#-------------------------------

# --- Standard start-up code; identical in all scripts---
source("Src/Reu/StandardizedStartupCode.R") #This runs the standardized code block at this location, ensuring identity across scripts. 
#This code reads the "report args" help ARG, and returns the args if it is active. 
#It then reads the filePrefix and sets the output directory to that file prefix, and reads the forceUpdate arg. 


# --- Argument Imports ---
{ # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
  # Defaults
  geneSets = NULL
  validGeneArguments = c("NULL", "gene", "Gene")
  pairwiseSets
  subdirectoryValueList = NULL
  usePermulations = F
  addPermulations = F
  writeData = T
  significanceCutoff = 0.05
  
  
  #Load from args 
  source("Src/Reu/ArgumentImport.R")
  
  #Have this script only run on folders with the pairwise comparison format
  pairwiseSets = pairwiseSets[grep("-", pairwiseSets)]
  #Convert geneset paths to geneset indicators
  geneSets = trimGmtFilenameToSetname(geneSets)

}

# ---- Main Code ----- 



# -- Make central RERData object 

getComparisionDifference = function(dataframe, colOne, colTwo){
  colOneIndex = names(dataframe)[which(names(dataframe) == colOne)]
  colTwoIndex = names(dataframe)[which(names(dataframe) == colTwo)]
  distanceFromEqual = abs(dataframe[colOneIndex] - dataframe[colTwoIndex]) / sqrt(2)
  distanceFromEqual
}

if(!usePermulations){
  pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "PairwiseCorrelationFile.rds", sep= "") #make a name for the pairwise comparisons based on prefix
  overallCorrelationFilename = paste(outputFolderName, filePrefix, "CorrelationFile.rds", sep= "") #make a name for the pairwise comparisons based on prefix
}else{
  pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "PermulationPairwiseCorrelationFile.rds", sep= "") #make a name for the pairwise comparisons based on prefix
  overallCorrelationFilename = paste(outputFolderName, filePrefix, "PermulationCorrelationFile.rds", sep= "") #make a name for the pairwise comparisons based on prefix
}

correlationResults = readRDS(pairwiseCorrelationFileName)
overallResults = readRDS(overallCorrelationFilename)
overall = list(overallResults)
names(overall)[[1]] = "Overall"
correlationResults = append(correlationResults, overall)
subdirectoryValueList = append(subdirectoryValueList, "Overall")

if(usingGene){
  combinedResults = NA
  combinedDrivers = NA
  combinedBinaries = NA
  prefixSet = NULL
  prefixList = NULL
  
  { # load in the data 
    for(i in 1:length(subdirectoryValueList)){
      currentSet = subdirectoryValueList[i]
      if(grepl(" ", currentSet)){
        correlationSubsetName = gsub("-", " - ", currentSet)
      }else{
        correlationSubsetName = currentSet
      }
      
      correlationPrefix = paste(substr(strsplit(currentSet, split = "-")[[1]],1,1), collapse = '')
      prefixEntry = correlationPrefix; names(prefixEntry) = currentSet; prefixSet = append(prefixSet, prefixEntry)
      for(j in 1:2){
        prefixSingle = strsplit(correlationPrefix, split = "")[[1]][j]; names(prefixSingle) = strsplit(currentSet, split="-")[[1]][j]; prefixList = append(prefixList, prefixSingle)
      }
      currentDataframe = which(names(correlationResults) == correlationSubsetName) 
      
      currentResults = correlationResults[[currentDataframe]]
      if(usePermulations){
        if(exists(currentResults$permP.adj)){
          currentResults$permSignificant = currentResults$permP.adj < significanceCutoff
          currentResults$unpermSignificant = currentResults$p.adj < significanceCutoff
        }else{
          currentResults$permSignificant = currentResults$permP < significanceCutoff
          currentResults$unpermSignificant = currentResults$p.adj < significanceCutoff
        }
        if(addPermulations){
          currentResults$significant = currentResults$unpermSignificant
        }else{
          currentResults$significant = currentResults$permSignificant
        }
      }else{
        currentResults$significant = currentResults$p.adj < significanceCutoff
      }
      
      
      names(currentResults) = paste0(correlationPrefix, "-", names(currentResults))
      combinedResults = cbind(combinedResults, currentResults)
      
      rm(currentResults)
      
    }
  }
  
  combinedResults = combinedResults[,-1]
  
  # -- save combination -- 
  if(writeData){
    combinedDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResults")
    write.csv(combinedResults, paste0(combinedDataFilename, ".csv"))
    saveRDS(combinedResults, paste0(combinedDataFilename, ".rds"))
  }
}


# --- Import GO Data --- 
if(usingGo){
  goSets = geneSets[-which(geneSets %in% validGeneArguments)]
  for(j in goSets){
    geneSet = j
    
    GOResults = list()
    for(i in 1:length(subdirectoryValueList)){
      currentSet = subdirectoryValueList[i]
      correlationSubsetName = gsub("-", " - ", currentSet)
      correlationPrefix = paste(substr(strsplit(currentSet, split = "-")[[1]],1,1), collapse = '')
      
      goFilename = paste0(outputFolderName, currentSet, "/", filePrefix, currentSet,"Enrichment-", geneSet, ".rds")
      tryCatch({
        currentGoData = readRDS(goFilename)[[1]]

        if(usePermulations){
          if(exists(currentGoData$permP.adj)){
            currentGoData$permSignificant = currentGoData$permP.adj < significanceCutoff
            currentGoData$unpermSignificant = currentGoData$p.adj < significanceCutoff        
          }else{
            currentGoData$permSignificant = currentGoData$permP < significanceCutoff
            currentGoData$unpermSignificant = currentGoData$p.adj < significanceCutoff
          }
          if(addPermulations){
            currentGoData$significant = currentGoData$unpermSignificant
          }else{
            currentGoData$significant = currentGoData$permSignificant
          }
        }else{
          currentGoData$significant = currentGoData$p.adj < significanceCutoff
        }
        
        
        
        names(currentGoData) = paste0(correlationPrefix, "-", names(currentGoData))
        
        
        GOResults[[i]] = currentGoData
        names(GOResults)[i] = correlationPrefix
      }, error = function(e) {
        message("Skipping file: ", goFilename, " | Error: ", conditionMessage(e))
      })
    }
    rm(currentGoData)
    missingGoResults = which(sapply(GOResults, is.null)) # Remove any empty lists generated by missing enrichment files 
    GOResults = GOResults[-missingGoResults]
    
    
    GoCombinedResults = NA
    longestGOResult = which(sapply(GOResults, nrow) == max(sapply(GOResults, nrow)))[1]
    for(i in 1:length(GOResults)){
      if(all(rownames(GOResults[[longestGOResult]]) == rownames(GOResults[[i]]))){
        cat("Combining GO Data", i, "\n")
        GoCombinedResults = cbind(GoCombinedResults, GOResults[[i]])
      }else{
        tryCatch({
          cat("Mismatch in order; reordering current data to first data; Combining GO Data", i, "\n")
          
          which(!rownames(GOResults[[longestGOResult]]) %in% rownames(GOResults[[i]]))
          
          matchedOrder = match(rownames(GOResults[[longestGOResult]]), rownames(GOResults[[i]]))
          #matchedOrder = matchedOrder[-which(is.na(matchedOrder))]
          GOResults[[i]] = GOResults[[i]][matchedOrder,]
          GoCombinedResults = cbind(GoCombinedResults, GOResults[[i]])
        },error = function(e){      
          stop("Rownames of GO results are not the same and cannot be reordered, there is an issue with the GO data.")
        }
        
        )
      }
    }
    GoCombinedResults = GoCombinedResults[,-1]
    rm(GOResults)
    
    
    
    if(writeData){
      combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResults-", geneSet)
      write.csv(GoCombinedResults, paste0(combinedGODataFilename, ".csv"))
      saveRDS(GoCombinedResults, paste0(combinedGODataFilename, ".rds"))
    }
  }
}

