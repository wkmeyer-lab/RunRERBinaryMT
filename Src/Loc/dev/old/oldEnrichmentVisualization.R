# Outdated code cut directly from the enrichment script. 


# --- Visualize the enrichment ----

{
  #This is manual only -- run-as-script does not accept a visualize output because no way to output result. 
  #For a script version, use PvQvGoVisualize.R 
  visualize = T
  visualize = F
  clean = T
  clean = F
  
  if(visualize){
    library(stringr)
    library(insight)
    enrichmentResult2 = enrichmentResult[[1]]
    makeGOTable = function(data, collumn){
      ValueHead = head(data[order(collumn, decreasing = T),], n=20)
      ValueHead$num.genes = as.character(ValueHead$num.genes)
      ValueHead$stat = round(ValueHead$stat, digits = 5)
      ValueHead$stat = as.character(ValueHead$stat)
      ValueHead = format_table(ValueHead, pretty_names = F, digits = "scientific5")
      ValueHead
    }
    enrichHead = makeGOTable(enrichmentResult2, abs(enrichmentResult2$stat))
    enrichHead
    textplot(enrichHead[1:4], mar = c(0,0,2,0), cmar = 1.5)
    if(usePermulations){
      title(main = paste("Top pathways by permulation"))
    }else{
      title(main = paste("Top pathways by non-permulation"))
    }
    
    if(clean){
      # ---- enrichment cleaning ----
      pathwaysToRemove = grep("CANCER", rownames(enrichHead))
      rowsToKeep = (!1:nrow(enrichHead) %in% pathwaysToRemove)
      cleanedHead = enrichHead[rowsToKeep,]
      textplot(cleanedHead[1:4], mar = c(0,0,2,0), cmar = 1.5)
      if(usePermulations){
        title(main = paste("Top pathways by permulation"))
      }else{
        title(main = paste("Top pathways by non-permulation"))
      }
      CleanheadFileName = paste(outputFolderName, filePrefix, "CleanedEnrichmentHead.csv", sep= "")
      write.csv(cleanedHead, CleanheadFileName)
    }
  }
}  