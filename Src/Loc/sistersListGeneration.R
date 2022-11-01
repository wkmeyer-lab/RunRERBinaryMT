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


inputTree= complexTree

fgEdges = which(complexTree$edge.length==1)
fgEdges
fgEdgeObjects = complexTree$edge[fgEdges,]
fgEdgeObjects

i=5


#Loop startup code
cladeNumber = 1
compeltedStartNodes = NULL
firstWrapperEdges = NULL
metaWrapperEdges = NULL
run =1 
for(i in 1:10){cladeName = paste("clade",i,sep=""); rm(list=cladeName)}
cladeList = NULL

#Loop code
for(i in fgEdges){
  #message(i)
  message(run)
  run=run+1
  startNode = inputTree$edge[i,1]
  startNode
  if(startNode %in% compeltedStartNodes){message("repeat"); next}
  compeltedStartNodes = append(compeltedStartNodes, startNode)
  
  endNodes = fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2]
  endNodes
  
  currentCladeName = paste("clade", cladeNumber, sep = "")
  currentCladeName
  
  if(!(length(endNodes)==1 | length(endNodes)==2)){
    message(paste("Something has gone wrong; foreground branch ", i, "has incorrect number of child nodes."))
  }
  
  
  
  #classify each edge as either being a solo clade, starter clade or a wrapper clade
  cladeType = NULL
  if(length(endNodes) == 1 & all(endNodes<length(inputTree$tip.label))){
    cladeType = "solo"
  }else if(all(endNodes < length(inputTree$tip.label))){
    cladeType = "starter"
  }else if(any(endNodes < length(inputTree$tip.label)) & any(endNodes > length(inputTree$tip.label))){
    cladeType = "firstWrapper"
  }else if(length(endNodes) == 1 & all(endNodes>length(inputTree$tip.label))){
    cladeType= "internal"
  }else if(all(endNodes > length(inputTree$tip.label))){
    cladeType = "metaWrapper"
  }
  message(cladeType)
  
  #Generate a clade vector based on clade type
  if(cladeType == "internal"){
    #These clades are only produced via the terminal/all/ancestral setup, and thus don't need to be factored in.
    next
  }
  
  if(cladeType == "solo"){
    speciesOne = inputTree$tip.label[endNodes[1]]
    assign(currentCladeName, c(speciesOne))
  }
  if(cladeType == "starter"){
    speciesOne = inputTree$tip.label[endNodes[1]]
    speciesTwo = inputTree$tip.label[endNodes[2]]
    assign(currentCladeName, c(speciesOne, speciesTwo))
    cladeEntry = startNode
    names(cladeEntry) = currentCladeName
    message(cladeEntry)
    cladeList = append(cladeList, cladeEntry)
  }
  if(cladeType == "firstWrapper"){
    firstWrapperEdges = append(firstWrapperEdges,i)     #These have to be done after all of the starters, so the index is saved.
    next
  }
  if(cladeType == "metaWrapper"){
    metaWrapperEdges = append(metaWrapperEdges,i)     #These have to be done after all of the firstwrappers, so the index is saved.
    next
  }
  cladeNumber = cladeNumber+1
}

#loop for first wrappers
for(i in firstWrapperEdges){
  message(run)
  run=run+1
  startNode = inputTree$edge[i,1]
  startNode
  
  endNodes = fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2]
  endNodes
  
  currentCladeName = paste("clade", cladeNumber, sep = "")
  currentCladeName
  
  if(!(length(endNodes)==1 | length(endNodes)==2)){
    message(paste("Something has gone wrong; foreground branch ", i, "has incorrect number of child nodes."))
  }
  
  #making new clade
  speciesNode = endNodes[which(endNodes<length(inputTree$tip.label))]
  speciesNode
  speciesName = inputTree$tip.label[speciesNode]
  speciesName
  cladeNode = endNodes[which(endNodes>length(inputTree$tip.label))]
  message(cladeNode)
  cladeName = names(which(cladeList == cladeNode))
  assign(currentCladeName, c(cladeName, speciesName))
  
  cladeEntry = startNode
  names(cladeEntry) = currentCladeName
  message(cladeEntry)
  cladeList = append(cladeList, cladeEntry)
  
  cladeNumber = cladeNumber+1
}

#loop for metaWrappers
remainingMetaWrapperEdges = metaWrapperEdges

while(length(remainingMetaWrapperEdges >0)){
  for(i in remainingMetaWrapperEdges){
    message(run)
    run=run+1
    startNode = inputTree$edge[i,1]
    startNode
    endNodes = fgEdgeObjects[which(fgEdgeObjects[,1] == startNode),2]
    endNodes
    currentCladeName = paste("clade", cladeNumber, sep = "")
    currentCladeName
    if(!(length(endNodes)==1 | length(endNodes)==2)){
      message(paste("Something has gone wrong; foreground branch ", i, "has incorrect number of child nodes."))
    }
    firstNode = endNodes[1]
    secondNode = endNodes[2]
    
    #if the two child nodes have not bee proccessed yet, skip this node for now
    if(!(firstNode %in% cladeList & secondNode %in% cladeList)){next}
    
    #otherwise, remove this node from the to-do list, and continue
    remainingMetaWrapperEdges = metaWrapperEdges[!metaWrapperEdges %in% i]
    
    firstCladeName = names(which(cladeList == firstNode))
    secondCladeName = names(which(cladeList == secondNode))
    
    assign(currentCladeName, c(firstCladeName, secondCladeName))
    
    cladeEntry = startNode
    names(cladeEntry) = currentCladeName
    message(cladeEntry)
    cladeList = append(cladeList, cladeEntry)
  }
}

#
 
# 

 # 
  
  
longEntry = c(1,2,3,4,5,60,70)
trimmedEntry = longEntry[!longEntry %in% 60]

  
if(any(endNodes < 31)){
  message("yes")
}else{
  message("no")
}
#


#










sisters_marine = list("clade1"=c("Killer_whale", "Dolphin"))

marineFgTree = foreground2TreeClades(marineFg,sisters_marine,toyTrees,plotTree=F)

marineFgTree$edge.length[45]=5

marineFg = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee")
marineplot1 = plotTreeHighlightBranches(marineFgTree,
                                        hlspecies=which(marineFgTree$edge.length==5),
                                        hlcols="blue", main="Marine mammals trait tree")
#
edgelabels(cex = 0.7, frame="none", font=2, adj=c(0,-0.7), col="blue")
dev.new(height=1000, width=30)
#
viewableMarineAll = minBranchLength(marineAll, 0.2)
plot(viewableMarineAll, edge.width = 2, label.offset = 0.1, type = "cladogram",)













which(!marineFgTree$tip.label %in% marineAll$tip.label)


which(marineAll$edge.length == 1 )
marineAll$edge.length[19]
marineAll$edge[19,]
marineAll$tip.label[9]

marineAll$edge.length[45]
marineAll$edge[45,]
marineAll$tip.label[9]

names(marineAll$e)

marineFgTree$edge[45]


which(marineAll$edge[,1] >marineAll$edge[,2])
which.max(marineAll$edge)
which(marineAll$tip.label == "Bushbaby")
which(marineAll$edge[,2] == 50)
marineAll$edge[104, ]






tree = read.tree(text = "(((A,B),(C,D)),E);")
plot(tree, type = "cladogram", edge.width = 2)

tree$edge
