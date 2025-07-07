#USAGE: 
#Input a tree you want tip labels changed to common name 
#Chose if you want a plot 
#Specify if foreground tree (only used for plotting)
#Will output tree of same shape with tips renamed

#Will expect the manual annotations spreadsheet to be placed in "Data/manualAnnotationsSheet.csv" by default, if not there, specify location 

#For manual use: 
#treeToConvertLocation = "Data/CVHRemakeBinaryForegroundTree.rds"
#inputTree = readRDS(treeToConvertLocation)

source("Src/Reu/ZonomNameConvertVectorCommon.R")

ZoonomTreeNameToCommon = function(tree, plot = T, isForegroundTree = T, manualAnnotLocation = "Data/mergedData.csv", hlcol = "blue", bgcol = "black", fontSize = 0.8, scientific = F, scientificCol = "ScientificName", commonCol = "CommonName", tipCol = "tipName"){
  
  inputTree = tree
  tipNames = inputTree$tip.label
  
  tipNames = ZonomNameConvertVectorCommon(tipNames, annotationLocation = manualAnnotLocation, toScientific = scientific, scientificColumn = scientificCol, commonColumn = commonCol, tipColumn = tipCol)
  
  inputTree$tip.label = tipNames
  
  if(plot){
    if(isForegroundTree){
      readableTree = inputTree
      readableTree$edge.length[readableTree$edge.length == 0] = 1
      plotTreeHighlightBranches2(readableTree, hlspecies = which(inputTree$edge.length == 1), hlcols = hlcol, bgcol = bgcol, fontSize = fontSize)
    }else{
      plotTree(inputTree)
    }
  }
  return(inputTree)
}




