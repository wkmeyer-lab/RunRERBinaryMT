#USAGE: 
# Convert names of one type to another type based on two columns in a specified spreadsheet. 





ZonomNameConvertVectorCommon = function(namesVector, 
                                        toScientific = F, 
                                        annotationLocation = "Data/mergedData.csv", 
                                        scientificColumn = NULL, 
                                        commonColumn = NULL, 
                                        tipColumn = NULL){
 
  #Handle null inputs of column names by using the default sets 
  {
  defaultScientificColumns =  c("ScientificName","Scientific Name","Scientific.Name",
                               "ScientificNames", "Scientific Names","Scientific.Names",
                               "ScientificNameFull", "Scientific Name Full", "Scientific.Name.Full",  
                               "scientificBinomial", "scientific Binomial", "scientific.Binomial", 
                               "scientificBinomials", "scientific Binomials", "scientific.Binomials")
  defaultCommonColumns =      c("CommonName", "Common Name", "Common.Name",
                                "CommonNames", "Common Names",  "Common.Names")
  defaultTipColumns =         c("tipName", "tip Name", "tip.Name",
                                "tipNames", "tip Names", "tip.Names",
                                "treeName", "tree Name", "tree.Name",
                                "treeNames", "tree Names", "tree.Names")
  
  if(is.null(scientificColumn)){scientificColumn = defaultScientificColumns}
  if(is.null(commonColumn)){commonColumn = defaultCommonColumns}
  if(is.null(tipColumn)){tipColumn = defaultTipColumns}
  }
  
  names = namesVector                                                    #make a vector of the names
  manualAnnot = read.csv(annotationLocation)                     #import manual annots file
  tipColumnNumber = which(names(manualAnnot) == tipColumn)
  if(toScientific){
    replaceColumnNumber = which(tolower(names(manualAnnot)) %in% tolower(scientificColumn))[1]
    nametype = "scientific"
  }else{
    replaceColumnNumber = which(tolower(names(manualAnnot)) %in% tolower(commonColumn))[1]
    nametype = "common"
  }
  
  if(!is.na(replaceColumnNumber)){
    for(i in 1:length(names)){                                                    #for each name: 
      currentName = names[i]                                                      #use the 'i'th name in the list
      currentRow = manualAnnot[manualAnnot[,tipColumnNumber] %in% currentName, ]             #find a row with the zonom name that matches the current name 
      currentSize = dim(currentRow)                                               #Part of "does row exist check": get the dimensions of the currentRow dataframe
      obsNumber = currentSize[1]                                                  #set "size" equal to the number of observations in 'currentRow'; which is the number of matches to the current name. If none exist it will be 0, if more than one it will be greater than 1. 
      if(obsNumber == 1){                                                         #if only one match exists:
        currentName = currentRow[[replaceColumnNumber]]                            #get the name from that row
      }else{
        currentName = names[i]                                                    #Otherwise, keep the name the same
      }
      names[i] = currentName                                                      #update the main name list with the name chose
    }
  }else{
    message(paste("ISSUE: No", nametype, "column found. This means either that no column within the defaults could be detected, or the supplied column was incorrect. \n The common tips have been kept identical to the main tips."))
  }
  

  #colnames(nMatrix) = names                                                     #update the matrix with the new names. 
  return(names)
}
