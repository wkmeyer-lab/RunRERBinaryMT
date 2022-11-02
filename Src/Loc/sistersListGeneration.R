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

edgelabels(cex = 0.7, frame="none", font=2, adj=c(0,-0.2), col="blue")
nodelabels(cex = 0.7, frame="none", font=2, adj=c(-0.2,0.3), col="dark green")
tiplabels(cex = 0.8, frame="none", font=2, adj=c(0.2,0), col="dark red")


batplot2 = plotTreeHighlightBranches(batTree,
                                         hlspecies=which(batTree$edge.length== 3),
                                         hlcols="blue", main="Marine mammals trait tree")


batTree = readRDS("Output/allInsectivory/allInsectivoryBinaryForegroundTree.rds")

inputTree= batTree
fgEdges = which(inputTree$edge.length==1)
fgEdges
fgEdgeObjects = inputTree$edge[fgEdges,]
fgEdgeObjects



#Loop startup code
cladNumber = 1
compeltedStartNodes = NULL
firstWrapperEdges = NULL
WrapperEdges = NULL
run =1 
for(i in 1:40){cladName = paste("clade",i,sep=""); rm(list=cladName)}
cladList = NULL

#Loop code
for(i in fgEdges){
  message("-----")
  message("First step run # ", run)
  message("Edge name ", i)
  run=run+1
  startNode = inputTree$edge[i,1]
  
  
  if(startNode %in% compeltedStartNodes){message("repeat"); next}
  compeltedStartNodes = append(compeltedStartNodes, startNode)
  
  message("Start node: ", startNode)
  
  endNodes = fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2]
  message("End nodes:", endNodes[1], " ", endNodes[2])
  
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
    
    cladEntry = startNode
    names(cladEntry) = currentCladName
    message("Clade entry:", cladEntry)
    message("Clade number: ", cladNumber)
    cladList = append(cladList, cladEntry)
  }
  if(cladType == "starter"){
    speciesOne = inputTree$tip.label[endNodes[1]]
    speciesTwo = inputTree$tip.label[endNodes[2]]
    assign(currentCladName, c(speciesOne, speciesTwo))
    
    cladEntry = startNode
    names(cladEntry) = currentCladName
    message("Clade entry:", cladEntry)
    message("Clade number: ", cladNumber)
    cladList = append(cladList, cladEntry)
  }
  if(cladType == "firstWrapper"){
    WrapperEdges = append(WrapperEdges, i)     #These have to be done after all of the starters, so the index is saved.
    next
  }
  if(cladType == "metaWrapper"){
    WrapperEdges = append(WrapperEdges, i)     #These have to be done after all of the firstwrappers, so the index is saved.
    next
  }
  cladNumber = cladNumber+1
}



#loop for metaWrappers
remainingWrapperEdges = WrapperEdges
stuckCycles = 0 
run = 1
while(length(remainingWrapperEdges >0) & stuckCycles < 20){
  for(i in remainingWrapperEdges){
    message("-----")
    message("Second step Run # ", run)
    message("Edge name: ", i)
    
    startRemainingEdges = length(remainingWrapperEdges)
    message(paste("Length of remaining Edge Wrappers = ", startRemainingEdges))
    message(paste("number of stuck cycles = ", stuckCycles))
    
    
    
    run=run+1
    startNode = inputTree$edge[i,1]
    message("startnode: ", startNode)
    endNodes = fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2]
    message("Endnodes: ", endNodes[1], " ", endNodes[2])
    currentCladName = paste("clade", cladNumber, sep = "")
    currentCladName
    
    if(!(length(endNodes)==1 | length(endNodes)==2)){
      message(paste("Something has gone wrong; foreground branch ", i, "has incorrect number of child nodes."))
    }
    
    firstNode = endNodes[1]
    firstNode
    firstNodeValid = (firstNode <= length(inputTree$tip.label) | firstNode %in% cladList )
    message("first node valid: ", firstNodeValid)
      
    secondNode = endNodes[2]
    secondNode
    secondNodeValid = (secondNode <= length(inputTree$tip.label) | secondNode %in% cladList )
    message("secondNodeValid: ", secondNodeValid)
    
    #if the two child nodes have not bee proccessed yet, skip this node for now
    if((firstNodeValid & secondNodeValid)){
    
    #otherwise, remove this node from the to-do list, and continue
    remainingWrapperEdges = remainingWrapperEdges[!remainingWrapperEdges %in% i]
    remainingWrapperEdges
    
    #get the first node's name
    if(firstNode <= length(inputTree$tip.label)){
      firstNodeName = inputTree$tip.label[firstNode]
    }else{
      firstNodeName = names(which(cladList == firstNode))
    }
    message("firstNodeName: ", firstNodeName)
    
    #get the second node's name
    if(secondNode <= length(inputTree$tip.label)){
      secondNodeName = inputTree$tip.label[secondNode]
    }else{
      secondNodeName = names(which(cladList == secondNode))
    }
    message("secondNodeName: ", secondNodeName)
    
    #Make the clade object
    assign(currentCladName, c(firstNodeName, secondNodeName))
    
    #Add the clade to the list
    cladEntry = startNode
    names(cladEntry) = currentCladName
    message("clade entry: ", cladEntry)
    message("clade number: ", cladNumber)
    cladList = append(cladList, cladEntry)
    
    cladNumber = cladNumber+1
    
    
    }else{
      message("Node not valid, skipping branch")
    }
    
    #check if it was stuck during this loop
    endRemainingEdges = length(remainingWrapperEdges)
    endRemainingEdges
    if(startRemainingEdges == endRemainingEdges){
      stuckCycles = stuckCycles+1
      stuckCycles
    }else{
      stuckCycles = 0
    }
  }
}
stuckCycles = 0 

#check for internal branches which are causing errors
run = 1
stuckCycles = 0 
while(length(remainingWrapperEdges >0) & stuckCycles < 20){
  for(i in remainingWrapperEdges){
    message("-----")
    message("Resolve Internals Run # ", run)
    message("Edge name: ", i)
    
    startRemainingEdges = length(remainingWrapperEdges)
    message(paste("Remaining Edges = ", remainingWrapperEdges))
    message(paste("number of stuck cycles = ", stuckCycles))
    
    
    
    run=run+1
    startNode = inputTree$edge[i,1]
    message("startnode: ", startNode)
    endNodes = fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2]
    message("Endnodes: ", endNodes[1], " ", endNodes[2])
    currentCladName = paste("clade", cladNumber, sep = "")
    currentCladName
    
    if(!(length(endNodes)==1 | length(endNodes)==2)){
      message(paste("Something has gone wrong; foreground branch ", i, "has incorrect number of child nodes."))
    }
    
    message("")
    firstNode = endNodes[1]
    message("First Node: ", firstNode)
    firstNodeChildren = inputTree$edge[which(inputTree$edge[,1] == firstNode),2]
    message("First Node Children: ", firstNodeChildren[1], " ", firstNodeChildren[2])
    
    if(any(!(firstNodeChildren %in% fgEdgeObjects))){
      if((!(firstNodeChildren[1] %in% fgEdgeObjects))){
        message(firstNodeChildren[1], " is not foreground")
        firstNode = firstNodeChildren[2]
        fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2][1] = firstNodeChildren[2]
      }
      if((!(firstNodeChildren[2] %in% fgEdgeObjects))){
        message(firstNodeChildren[2], " is not foreground")
        firstNode = firstNodeChildren[1]
        fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2][1] = firstNodeChildren[1]
      }
    }
    
    message("")
    secondNode = endNodes[2]
    message("Second Node: ", secondNode)
    secondNodeChildren = inputTree$edge[which(inputTree$edge[,1] == secondNode),2]
    message("Second Node Children: ", secondNodeChildren[1], " ", secondNodeChildren[2])
    
    if(any(!(secondNodeChildren %in% fgEdgeObjects))){
      if((!(secondNodeChildren[1] %in% fgEdgeObjects))){
        message(secondNodeChildren[1], " is not foreground")
        secondNode = secondNodeChildren[2]
        fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2][2] = secondNodeChildren[2]
        }
      if((!(secondNodeChildren[2] %in% fgEdgeObjects))){
        message(secondNodeChildren[2], " is not foreground")
        secondNode = secondNodeChildren[1]
        fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2][2] = secondNodeChildren[1]
        }
    }
    
    firstNodeValid = (firstNode <= length(inputTree$tip.label) | firstNode %in% cladList )
    message("first node valid: ", firstNodeValid)
  
    secondNodeValid = (secondNode <= length(inputTree$tip.label) | secondNode %in% cladList )
    message("secondNodeValid: ", secondNodeValid)
    
    #if the two child nodes have not bee proccessed yet, skip this node for now
    if((firstNodeValid & secondNodeValid)){
      
      #otherwise, remove this node from the to-do list, and continue
      remainingWrapperEdges = remainingWrapperEdges[!remainingWrapperEdges %in% i]
      remainingWrapperEdges
      
      #get the first node's name
      if(firstNode <= length(inputTree$tip.label)){
        firstNodeName = inputTree$tip.label[firstNode]
      }else{
        firstNodeName = names(which(cladList == firstNode))
      }
      message("firstNodeName: ", firstNodeName)
      
      #get the second node's name
      if(secondNode <= length(inputTree$tip.label)){
        secondNodeName = inputTree$tip.label[secondNode]
      }else{
        secondNodeName = names(which(cladList == secondNode))
      }
      message("secondNodeName: ", secondNodeName)
      
      #Make the clade object
      assign(currentCladName, c(firstNodeName, secondNodeName))
      
      #Add the clade to the list
      cladEntry = startNode
      names(cladEntry) = currentCladName
      message("clade entry: ", cladEntry)
      message("clade number: ", cladNumber)
      cladList = append(cladList, cladEntry)
      
      cladNumber = cladNumber+1
      
      
    }else{
      message("Node not valid, skipping branch")
    }
    
    #check if it was stuck during this loop
    endRemainingEdges = length(remainingWrapperEdges)
    endRemainingEdges
    if(startRemainingEdges == endRemainingEdges){
      stuckCycles = stuckCycles+1
      stuckCycles
    }else{
      stuckCycles = 0
    }
  }
}

cladObjectSet = ls(pattern = "clade")
sistersListExport = list(cladObjectsList)
#










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
 
# ?ls
