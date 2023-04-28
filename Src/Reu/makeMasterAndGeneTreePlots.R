makeMasterVsGeneTreePlots = function(mainTrees, RERObject, geneInQuestion, foregroundVector, colors = "blue"){
  source("Src/Reu/ZonomNameConvertVector.R")
  source("Src/Reu/ZoonomTreeNameToCommon.R")
  source("Src/Reu/GetForegroundEdges.R")
  foregroundVector = ZonomNameConvertVectorCommon(foregroundVector)
  RERObject[geneInQuestion,]
  
  namesToKeep = names(RERObject[geneInQuestion,][!is.na(RERObject[geneInQuestion,])])
  namesToKeep = namesToKeep[!is.na(namesToKeep)]
  
  masterTree = mainTrees$masterTree
  prunedMaster = drop.tip(masterTree, which(!masterTree$tip.label %in% namesToKeep))
  geneTree = mainTrees$trees[[geneInQuestion]]
  prunedGeneTree = drop.tip(geneTree, which(!geneTree$tip.label %in% namesToKeep))
  
  commonMaster = ZoonomTreeNameToCommon(prunedMaster)
  commonGene = ZoonomTreeNameToCommon(prunedGeneTree)
  
  par(mfrow = c(1,2))
  masterFGEdges = getForegroundEdges(commonMaster, foregroundVector)
  plotTreeHighlightBranches(commonMaster, hlspecies = masterFGEdges, main = "Overall Genome", hlcols = colors)
  
  geneFGEdges = getForegroundEdges(commonGene, foregroundVector)
  plotTreeHighlightBranches(commonGene, main = geneInQuestion, hlspecies = geneFGEdges, hlcols = colors)
  
}
