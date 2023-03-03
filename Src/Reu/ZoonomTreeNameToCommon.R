#USAGE: 
#Input a tree you want tip labels changed to common name 
#Chose if you want a plot 
#Specify if foreground tree (only used for plotting)
#Will output tree of same shape with tips renamed

#Will expect the manual annotations spreadsheet to be placed in "Data/manualAnnotationsSheet.csv" by default, if not there, specify location 

#For manual use: 
#treeToConvertLocation = "Data/CVHRemakeBinaryForegroundTree.rds"
#inputTree = readRDS(treeToConvertLocation)

ZoonomTreeNameToCommon = function(tree, plot = T, isForegroundTree = T, manualAnnotLocation = "Data/manualAnnotationsSheet.csv"){
  library(RERconverge)
  
  manualAnnot = read.csv(manualAnnotLocation) 
  inputTree = tree
  tipNames = inputTree$tip.label
  for(i in 1:length(tipNames)){                                                    #for each name: 
    currentName = tipNames[i]                                                      #use the 'i'th name in the list
    currentRow = manualAnnot[manualAnnot$FaName %in% currentName, ]             #find a row with the zonom name that matches the current name 
    currentSize = dim(currentRow)                                               #Part of "does row exist check": get the dimensions of the currentRow dataframe
    obsNumber = currentSize[1]                                                  #set "size" equal to the number of observations in 'currentRow'; which is the number of matches to the current name. If none exist it will be 0, if more than one it will be greater than 1. 
    if(obsNumber == 1){                                                         #if only one match exists:
      currentName = currentRow$Common.Name.or.Group                             #get the name from that row
    }else{
      currentName = tipNames[i]                                                    #Otherwise, keep the name the same
    }
    tipNames[i] = currentName                                                      #update the main name list with the name chose
  }
  inputTree$tip.label = tipNames
  
  if(plot){
    if(isForegroundTree){
      readableTree = inputTree
      readableTree$edge.length[readableTree$edge.length == 0] = 1
      plotTreeHighlightBranches(readableTree, hlspecies = which(inputTree$edge.length == 1), hlcols = "red")
    }else{
      plotTree(inputTree)
    }
  }
  return(inputTree)
}
