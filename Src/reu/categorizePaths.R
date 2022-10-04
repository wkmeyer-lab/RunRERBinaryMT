categorizePaths = function(binaryPhenotypeTree, zonomMasterObject, prefix, overwrite = F){
  filename = paste("Results/", prefix, "ManualFGTree.rds", sep='')
  message(filename)
  if(file.exists(filename) && overwrite == F){
    manualFGTree = readRDS(filename)
  }else{
    source("Src/Reu/ZonomNameConvertVector.R")
    binaryTree = binaryPhenotypeTree
    renamedTree = binaryTree
    renamedTree$tip.label = ZonomNameConvertVectorCommon(binaryTree$tip.label)
    manualFGTree = click_select_foreground_branches(renamedTree)
    manualFGTree$tip.label = binaryTree$tip.label
    saveRDS(manualFGTree, filename)
  }
  manualPaths = tree2Paths(manualFGTree, zonomMasterObject)
  return(manualPaths)
}