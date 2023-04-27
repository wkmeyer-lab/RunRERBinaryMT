makeMasterVsGeneTreePlots = function(mainTrees, RERObject, geneInQuestion, foregroundVector){
  source("Src/Reu/ZonomNameConvertVector.R")
  source("Src/Reu/ZoonomTreeNameToCommon.R")
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
  plotTreeHighlightBranches(commonMaster, hlspecies = foregroundVector, main = "Overall Genome", hlcols = "blue")
  title()
  plotTreeHighlightBranches(commonGene, main = geneInQuestion, hlspecies = foregroundVector, hlcols = "blue")
  title(main = geneInQuestion)
}
