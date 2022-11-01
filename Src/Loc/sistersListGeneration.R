#Library setup 
#.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge) #load RERconverge package
library(RERconverge)
library("tools")
library(paleotree)
source("Src/Reu/cmdArgImport.R")
source("Src/Reu/convertLogiToNumeric.R")


#Setup toy debug file
rerpath = find.package('RERconverge')
toytreefile = "subsetMammalGeneTrees.txt" 
toyTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)

marineExtantForground = c("Walrus", "Seal", "Killer_whale", "Dolphin", "Manatee")

#






#
marineAll = foreground2Tree(marineExtantForground, toyTrees, clade = "all", useSpecies = names(logAdultWeightcm))
complexFG = c("Golden_hamster","Chinese_hamster","Vole","Walrus", "Seal", "Killer_whale", "Dolphin", "Manatee", "Rat", "Mouse")

complexTree  = foreground2Tree(complexFG, toyTrees, clade = "all", useSpecies = names(logAdultWeightcm))
marineplot2 = plotTreeHighlightBranches(marineAll,
                                        hlspecies=which(marineAll$edge.length== 3),
                                        hlcols="blue", main="Marine mammals trait tree")

complexplot2 = plotTreeHighlightBranches(complexTree,
                                        hlspecies=which(complexTree$edge.length== 3),
                                        hlcols="blue", main="Marine mammals trait tree")

edgelabels(cex = 0.7, frame="none", font=2, adj=c(0,-0.7), col="blue")
nodelabels(cex = 0.7, frame="none", font=2, adj=c(-0.2,0.3), col="dark green")
tiplabels(cex = 0.8, frame="none", font=2, adj=c(0.2,0), col="dark red")

batTree = readRDS("Output/allInsectivory/allInsectivoryBinaryForegroundTree.rds")
batplot2 = plotTreeHighlightBranches(batTree,
                                         hlspecies=which(batTree$edge.length== 3),
                                         hlcols="blue", main="Marine mammals trait tree")


inputTree= batTree
fgEdges = which(inputTree$edge.length==1)
fgEdges
fgEdgeObjects = inputTree$edge[fgEdges,]
fgEdgeObjects



#Loop startup code
cladNumber = 1
compeltedStartNodes = NULL
firstWrapperEdges = NULL
metaWrapperEdges = NULL
run =1 
for(i in 1:10){cladName = paste("clade",i,sep=""); rm(list=cladName)}
cladList = NULL

#Loop code
for(i in fgEdges){
  message(i)
  message(run)
  run=run+1
  startNode = inputTree$edge[i,1]
  startNode
  if(startNode %in% compeltedStartNodes){message("repeat"); next}
  compeltedStartNodes = append(compeltedStartNodes, startNode)
  
  endNodes = fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2]
  endNodes
  
  currentCladName = paste("clade", cladNumber, sep = "")
  currentCladName
  
  if(!(length(endNodes)==1 | length(endNodes)==2)){
    message(paste("Something has gone wrong; foreground branch ", i, "has incorrect number of child nodes."))
  }
  
  
  
  #classify each edge as either being a solo clade, starter clade or a wrapper clade
  cladType = NULL
  if(length(endNodes) == 1 & all(endNodes<length(inputTree$tip.label))){
    cladType = "solo"
  }else if(all(endNodes <= length(inputTree$tip.label))){
    cladType = "starter"
  }else if(any(endNodes <= length(inputTree$tip.label)) & any(endNodes > length(inputTree$tip.label))){
    cladType = "firstWrapper"
  }else if(length(endNodes) == 1 & all(endNodes>length(inputTree$tip.label))){
    cladType= "internal"
  }else if(all(endNodes > length(inputTree$tip.label))){
    cladType = "metaWrapper"
  }
  message(cladType)
  
  #Generate a clade vector based on clade type
  if(cladType == "internal"){
    #These clades are only produced via the terminal/all/ancestral setup, and thus don't need to be factored in.
    next
  }
  
  if(cladType == "solo"){
    speciesOne = inputTree$tip.label[endNodes[1]]
    assign(currentCladName, c(speciesOne))
  }
  if(cladType == "starter"){
    speciesOne = inputTree$tip.label[endNodes[1]]
    speciesTwo = inputTree$tip.label[endNodes[2]]
    assign(currentCladName, c(speciesOne, speciesTwo))
    cladEntry = startNode
    names(cladEntry) = currentCladName
    message(cladEntry)
    cladList = append(cladList, cladEntry)
  }
  if(cladType == "firstWrapper"){
    firstWrapperEdges = append(firstWrapperEdges,i)     #These have to be done after all of the starters, so the index is saved.
    next
  }
  if(cladType == "metaWrapper"){
    metaWrapperEdges = append(metaWrapperEdges,i)     #These have to be done after all of the firstwrappers, so the index is saved.
    next
  }
  cladNumber = cladNumber+1
}

#loop for first wrappers
for(i in firstWrapperEdges){
  message(run)
  run=run+1
  startNode = inputTree$edge[i,1]
  startNode
  
  endNodes = fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2]
  endNodes
  
  currentCladName = paste("clade", cladNumber, sep = "")
  currentCladName
  
  if(!(length(endNodes)==1 | length(endNodes)==2)){
    message(paste("Something has gone wrong; foreground branch ", i, "has incorrect number of child nodes."))
  }
  
  #making new clade
  speciesNode = endNodes[which(endNodes <= length(inputTree$tip.label))]
  speciesNode
  speciesName = inputTree$tip.label[speciesNode]
  speciesName
  cladNode = endNodes[which(endNodes > length(inputTree$tip.label))]
  message(cladNode)
  cladName = names(which(cladList == cladNode))
  assign(currentCladName, c(cladName, speciesName))
  
  cladEntry = startNode
  names(cladEntry) = currentCladName
  message(cladEntry)
  cladList = append(cladList, cladEntry)
  
  cladNumber = cladNumber+1
}

#loop for metaWrappers
remainingMetaWrapperEdges = metaWrapperEdges
stuckCycles = 0 

if(length(remainingMetaWrapperEdges >0 & stuckCycles < 10)){
  message(paste("Length of remaining Edge Wrappers = ", length(remainingMetaWrapperEdges)))
  startRemainingEdges = remainingMetaWrapperEdges
  for(i in remainingMetaWrapperEdges){
    message(run)
    run=run+1
    startNode = inputTree$edge[i,1]
    startNode
    endNodes = fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2]
    endNodes
    currentCladName = paste("clade", cladNumber, sep = "")
    currentCladName
    if(!(length(endNodes)==1 | length(endNodes)==2)){
      message(paste("Something has gone wrong; foreground branch ", i, "has incorrect number of child nodes."))
    }
    firstNode = endNodes[1]
    secondNode = endNodes[2]
    
    #if the two child nodes have not bee proccessed yet, skip this node for now
    if(!(firstNode %in% cladList & secondNode %in% cladList)){next}
    
    #otherwise, remove this node from the to-do list, and continue
    remainingMetaWrapperEdges = metaWrapperEdges[!metaWrapperEdges %in% i]
    
    firstCladName = names(which(cladList == firstNode))
    secondCladName = names(which(cladList == secondNode))
    
    assign(currentCladName, c(firstCladName, secondCladName))
    
    cladEntry = startNode
    names(cladEntry) = currentCladName
    message(cladEntry)
    cladList = append(cladList, cladEntry)
  }
  endRemainingEdges = remainingMetaWrapperEdges
  if(startRemainingEdges == endRemainingEdges){
    stuckCycles = stuckCycles+1
  }else{
    stuckCycles = 0
  }
  
}
cladObjectSet = ls(pattern = "clade")
sistersListExport = list(cladObjectsList)
#
 
# ?ls
