source("Src/Reu/ZoonomTreeNameToCommon.R")

autopruner= function(masterTree, dropPercent = NA, dropValue = 0.01, tipsToKeep = NA, returnEdgeTable = F, procedurePlot = F, nameConversionColumn = NA, nameConversionData = "Data/mergedData.csv", preDroppedTips = NA, originalTree = NA){
  if(all(is.na(originalTree))){originalTree = masterTree}
  
  # -- determine the branch length cutoff -- 
  if(!is.na(dropPercent)){
    if(dropPercent > 1){dropPercent = dropPercent/100}
    lengthCutoff = quantile(masterTree$edge.length, dropPercent)  
  }else{
    lengthCutoff = dropValue
  }
  
  # -- drop tips until the branch lengths are within the cutoff -- 
  trimmedTree = masterTree 
  terminalEdges = which(trimmedTree$edge[,2] <= length(trimmedTree$tip.label))
  terminalEdges = which(trimmedTree$edge[,2] <= length(trimmedTree$tip.label))
  terminalEdgeTable = data.frame(edgeNumber = terminalEdges, edgeLength = trimmedTree$edge.length[terminalEdges], edgeStart = trimmedTree$edge[terminalEdges,1], edgeEnd = trimmedTree$edge[terminalEdges,2])
  terminalEdgeTable$edgeTip = trimmedTree$tip.label[terminalEdgeTable$edgeEnd]
  terminalEdgeTable = terminalEdgeTable[order(terminalEdgeTable$edgeLength),]
  if(returnEdgeTable){
    return(terminalEdgeTable)
    stop()
  }
  
  #message("Species in order of shortest tips:")
  #message(paste(terminalEdgeTable$edgeTip, collapse="\n"))
  
  droppedTips = vector()
  skippedTips = vector()
  
  while(min(terminalEdgeTable$edgeLength) <= lengthCutoff){
    {
    selectingBranch = T; i = 1
    while(selectingBranch){
      shortestTerminalBranch = terminalEdgeTable[i,]
      if(shortestTerminalBranch$edgeTip %in% tipsToKeep){
        if(!shortestTerminalBranch$edgeTip %in% skippedTips){
          message(paste("Skipping", shortestTerminalBranch$edgeTip))
          skippedTips = append(skippedTips, shortestTerminalBranch$edgeTip)
        }
        i = i+1
      }else{
        message(paste("Pruning", shortestTerminalBranch$edgeTip, " -- length:", shortestTerminalBranch$edgeLength))
        i = 1
        droppedTips = append(droppedTips, shortestTerminalBranch$edgeTip)
        selectingBranch = F
      }
    }
    trimmedTree = drop.tip(trimmedTree, shortestTerminalBranch$edgeTip)
    
    terminalEdges = which(trimmedTree$edge[,2] <= length(trimmedTree$tip.label))
    terminalEdgeTable = data.frame(edgeNumber = terminalEdges, edgeLength = trimmedTree$edge.length[terminalEdges], edgeStart = trimmedTree$edge[terminalEdges,1], edgeEnd = trimmedTree$edge[terminalEdges,2])
    terminalEdgeTable$edgeTip = trimmedTree$tip.label[terminalEdgeTable$edgeEnd]
    if(!length(which(terminalEdgeTable$edgeTip %in% skippedTips))==0){
      terminalEdgeTable = terminalEdgeTable[-which(terminalEdgeTable$edgeTip %in% skippedTips),] #This prevents the while loop from continuing due to a skipped tip 
    }
    terminalEdgeTable = terminalEdgeTable[order(terminalEdgeTable$edgeLength),]
    
    if(procedurePlot){
      tipDropPlot(nameConversionColumn, nameConversionData, originalTree, trimmedTree, droppedTips, skippedTips, preDroppedTips, message = F)
    }
    } 
    min(trimmedTree$edge.length[terminalEdges])
  }

  
  message(paste("Skipped tips:"))
  message(paste(skippedTips, collaspe = ", ", sep=""))
  message(paste("Dropped tips:"))
  message((paste(droppedTips, collaspe = ", ", sep="")))  
  
  if(!is.na(nameConversionColumn)){tipDropPlot(nameConversionColumn, nameConversionData, originalTree, trimmedTree, droppedTips, skippedTips, preDroppedTips, message = T)}
  
  par(mfrow = c(1,2))
  plotTreeHighlightBranches(originalTree, hlspecies = droppedTips, hlcols = "red")
  plot.phylo(trimmedTree)
  par(mfrow = c(1,1))
  
  droppedTips = append(droppedTips, preDroppedTips)
  droppedTips <<- droppedTips
  trimmedTree
}

tipDropPlot = function(nameConversionColumn, nameConversionData, masterTree, trimmedTree, droppedTips, skippedTips, preDroppedTips, message = F){
  
  if(!all(is.na(preDroppedTips))){
    droppedTips = append(droppedTips, preDroppedTips)
  }
  
  if(!is.na(nameConversionColumn)){
    plotMasterTree = ZoonomTreeNameToCommon(masterTree, manualAnnotLocation = nameConversionData, tipCol = nameConversionColumn, plot = F)
    plotTrimmedTree = ZoonomTreeNameToCommon(trimmedTree, manualAnnotLocation = nameConversionData, tipCol = nameConversionColumn, plot = F)
    plotDroppedTips = ZonomNameConvertVectorCommon(droppedTips, annotationLocation = nameConversionData, tipCol = nameConversionColumn)
    plotSkippedTips = ZonomNameConvertVectorCommon(skippedTips, annotationLocation = nameConversionData, tipCol = nameConversionColumn)
  }else{
    plotMasterTree = masterTree
    plotTrimmedTree = trimmedTree
    plotDroppedTips = droppedTips
    plotSkippedTips = skippedTips
  }
  
  par(mfrow = c(1,2))
  plotTreeHighlightBranches(plotMasterTree, hlspecies = plotDroppedTips, hlcols = "red")
  plot.phylo(plotTrimmedTree)
  par(mfrow = c(1,1))
  
  if(message){
    message(paste("Skipped tips:"))
    message(paste(plotSkippedTips, collaspe = ", ", sep=""))
    message(paste("Dropped tips:"))
    message((paste(plotDroppedTips, collaspe = ", ", sep="")))
  }
}

