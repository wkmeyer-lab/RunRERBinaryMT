#This script takes an input multiphylo and gene name, and outputs 1) a tree based on the master-tree branch lengths with the gene tree topology and 2) the gene tree. 
  #If an RER object is provided, it will instead trim to only tips which have RER values. 
  #If provided a foreground vector, it can color the foreground and background different colors. 
    #This script is no longer dependent on RERConverge
#This script will by default convert names from zoonomia names to common names. This requires "Data/manualAnnotationsSheet.csv". This can be toggled off using convertNames = F=
#This script can also make a plot of overall genome length vs gene length by toggling correlationPlot = T.

makeMasterAndGeneTreePlots = function(mainTrees, geneInQuestion, RERObject = NULL, foregroundVector = NULL, fgcols = "blue", correlationPlot = F, bgcolor = "black", rmlabels = NULL, convertNames = T, twoInOne = T){
  masterTree = mainTrees$masterTree
  geneTree = mainTrees$trees[[geneInQuestion]]
  
  if(!is.null(RERObject)){
  namesToKeep = names(RERObject[geneInQuestion,][!is.na(RERObject[geneInQuestion,])])
  namesToKeep = namesToKeep[!is.na(namesToKeep)]
  }else{
    namesToKeep = geneTree$tip.label
  }
  
  prunedMaster = drop.tip(masterTree, which(!masterTree$tip.label %in% namesToKeep))
  prunedGeneTree = drop.tip(geneTree, which(!geneTree$tip.label %in% namesToKeep))
  
  if(convertNames){
    if(file.exists("Src/Reu/ZonomNameConvertVector.R")){source("Src/Reu/ZonomNameConvertVector.R")}
    if(file.exists("Src/Reu/ZoonomTreeNameToCommon.R")){source("Src/Reu/ZoonomTreeNameToCommon.R")}
    commonMaster = ZoonomTreeNameToCommon(prunedMaster, plot = F)
    commonGene = ZoonomTreeNameToCommon(prunedGeneTree, plot = F)
    plotMaster = commonMaster
    plotGene = commonGene
  }else{
    plotMaster = prunedMaster
    plotGene = prunedGeneTree
  }
  
  if(!is.null(foregroundVector)){
    if(file.exists("Src/Reu/GetForegroundEdges.R")){source("Src/Reu/GetForegroundEdges.R")}
    if(convertNames){foregroundVector = ZonomNameConvertVectorCommon(foregroundVector)}
    masterFGEdges = getForegroundEdges(plotMaster, foregroundVector)
    geneFGEdges = getForegroundEdges(plotGene, foregroundVector)
  }else{
    masterFGEdges = ""
    geneFGEdges = ""
  }
  
  if(twoInOne){
    par(mfrow = c(1,2), mai = c(0.5, 0.1, 0.2, 0.1))
  }
  plotTreeHighlightBranches2(plotMaster, hlspecies = masterFGEdges, main = "Overall Genome Average", hlcols = fgcols, bgcol = bgcolor)
  plotTreeHighlightBranches2(plotGene, main = paste("Gene:", geneInQuestion), hlspecies = geneFGEdges, hlcols = fgcols, bgcol = bgcolor)
  
  if(correlationPlot & !is.null(foregroundVector)){
    if(all.equal(masterFGEdges, geneFGEdges)){ #Only run this code if the trees are the same shape
      
      if(!is.null(rmlabels)){
        commonMaster$tip.label[which(commonMaster$tip.label %in% rmlabels)] = "   "
        commonGene$tip.label[which(commonGene$tip.label %in% rmlabels)] = "   "
      }
      
      edgeGround = rep(NA, length(commonMaster$edge.length))
      edgeGround[masterFGEdges] = 1
      edgeGround[!1:length(edgeGround) %in% masterFGEdges] = 0
      
      endNodes = commonMaster$edge[,2]
      edgeNames = commonMaster$tip.label[endNodes]
      
      edgeCorrealtions = data.frame(commonMaster$edge.length, commonGene$edge.length, as.factor(edgeGround), edgeNames)
      names(edgeCorrealtions) = c("Master", "Gene", "Ground", "Names")
      
      correlPlot <- ggplot(edgeCorrealtions, 
                  aes(x = Master, 
                      y = Gene, 
                      col = Ground, 
                      label = Names)
                  ) + 
                  scale_size_manual(values = c(1, 1, 1, 1)) + 
                  geom_point(aes(size = Ground)) + 
                  scale_color_manual(values = c(bgcolor, fgcols)) + 
                  coord_fixed() + 
                  expand_limits(x = max(edgeCorrealtions$Gene), y = max(edgeCorrealtions$Master))+
                  geom_text(hjust = "bottom", size = 4, check_overlap = T, hjust = 1) + 
                  ylab("Gene-Specific Branch Length") + 
                  xlab("Genome Average Branch Length") + 
                  #ggtitle(plottitle) + 
                  geom_abline(intercept = 0, slope =1, linetype = "dotted") + 
                  theme(
                      #axis.ticks.y = element_blank(), 
                      #axis.text.y = element_blank(), 
                      legend.position = "none", 
                      panel.background = element_blank(), 
                      axis.text = element_text(size = 12, colour = "black"), 
                      axis.title = element_text(size = 20, face = "bold"), 
                      plot.title = element_text(size = 24, face = "bold")
                      ) +
                  theme(axis.line = element_line(colour = "black", size = 0.2))  
                  #theme(axis.line.y = element_blank())
      
      print(correlPlot)
    }
  }
  
}


plotTreeHighlightBranches2 =function (tree, outgroup = NULL, hlspecies, hlcols = NULL, main = "", 
                                      useGG = FALSE, bgcol = "black") 
{
  if (is.null(hlcols)) {
    hlcols <- c(2:(length(hlspecies) + 1))
  }
  if (length(hlcols) < length(hlspecies)) {
    hlcols <- rep_len(hlcols, length(hlspecies))
  }
  if (!is.null(outgroup)) {
    outgroup <- outgroup[outgroup %in% tree$tip.label]
    if (length(outgroup) > 0) {
      if (is.numeric(hlspecies)) {
        dummytree <- tree
        dummytree$edge.length <- c(rep(1, nrow(dummytree$edge)))
        for (i in c(1:length(hlspecies))) {
          dummytree$edge.length[hlspecies] <- i + 1
        }
        dummyrooted <- root(dummytree, outgroup)
      }
      rooted <- root(tree, outgroup)
    }
    else {
      print("No members of requested outgroup found in tree; keeping unrooted.")
      rooted <- tree
      outgroup <- NULL
    }
  }
  else {
    rooted <- tree
  }
  if (useGG) {
    if (is.numeric(hlcols)) {
      hlcols <- rep_len("#0000ff", length(hlspecies))
    }
    hlspecies_named <- vector(mode = "character")
    if (is.numeric(hlspecies)) {
      for (i in 1:length(hlspecies)) {
        hlspecies_named[i] <- tree$tip.label[hlspecies[i]]
      }
    }
    else hlspecies_named <- hlspecies
    rooted2 <- rooted
    mm <- min(rooted2$edge.length[rooted2$edge.length > 0])
    rooted2$edge.length[rooted2$edge.length == 0] <- max(0.02, 
                                                         mm/20)
    tipCols <- vector(mode = "character")
    x <- 1
    for (i in 1:length(tree$tip.label)) {
      if (tree$tip.label[i] %in% hlspecies_named) {
        tipCols[i] <- hlcols[x]
        x <- x + 1
      }
      else tipCols[i] <- bgcol
    }
    nlabel <- rooted2$Nnode + length(rooted2$tip.label)
    edgeCols <- vector(mode = "character", length = nlabel)
    x <- 1
    for (i in 1:nlabel) {
      if (i < length(rooted2$tip.label)) {
        if (rooted2$tip.label[i] %in% hlspecies_named) {
          edgeCols[i] <- hlcols[x]
          x <- x + 1
        }
        else edgeCols[i] <- bgcol
      }
      else edgeCols[i] <- bgcol
    }
    plotobj = ggtree(rooted2, color = edgeCols)
    plotobj = plotobj + geom_tiplab(color = tipCols, geom = "text", 
                                    cex = 3) + labs(title = main)
    return(plotobj)
  }
  else {
    colMaster <- c(rep(bgcol, nrow(rooted$edge)))
    if (is.numeric(hlspecies)) {
      if (!is.null(outgroup)) {
        hlcols <- c(bgcol, hlcols)
        colMaster <- hlcols[dummyrooted$edge.length]
      }
      else {
        for (i in 1:length(hlspecies)) {
          colMaster[hlspecies[i]] <- hlcols[i]
        }
      }
    }
    else {
      wspmr <- rooted$tip.label[rooted$edge[, 2]]
      for (i in 1:length(hlspecies)) {
        colMaster[which(wspmr == hlspecies[i])] <- hlcols[i]
      }
    }
    termedge <- order(rooted$edge[, 2])[1:length(rooted$tip.label)]
    colMasterTip <- colMaster[termedge]
    rooted2 = rooted
    mm = min(rooted2$edge.length[rooted2$edge.length > 0])
    rooted2$edge.length[rooted2$edge.length == 0] = max(0.02, 
                                                        mm/20)
    plotobj = plot.phylo(rooted2, main = main, edge.color = colMaster, 
                         tip.color = colMasterTip, edge.width = 2, cex = 0.8)
    return(plotobj)
  }
}


ZonomNameConvertVectorCommon = function(namesVector){
  names = namesVector                                                    #make a vector of the names
  manualAnnot = read.csv("Data/manualAnnotationsSheet.csv")                     #improt manual annots file
  for(i in 1:length(names)){                                                    #for each name: 
    currentName = names[i]                                                      #use the 'i'th name in the list
    currentRow = manualAnnot[manualAnnot$FaName %in% currentName, ]             #find a row with the zonom name that matches the current name 
    currentSize = dim(currentRow)                                               #Part of "does row exist check": get the dimensions of the currentRow dataframe
    obsNumber = currentSize[1]                                                  #set "size" equal to the number of observations in 'currentRow'; which is the number of matches to the current name. If none exist it will be 0, if more than one it will be greater than 1. 
    if(obsNumber == 1){                                                         #if only one match exists:
      currentName = currentRow$Common.Name.or.Group                             #get the name from that row
    }else{
      currentName = names[i]                                                    #Otherwise, keep the name the same
    }
    names[i] = currentName                                                      #update the main name list with the name chose
  }
  #colnames(nMatrix) = names                                                     #update the matrix with the new names. 
  return(names)
}

#USAGE: 
#Input a tree you want tip labels changed to common name 
#Chose if you want a plot 
#Specify if foreground tree (only used for plotting)
#Will output tree of same shape with tips renamed

#Will expect the manual annotations spreadsheet to be placed in "Data/manualAnnotationsSheet.csv" by default, if not there, specify location 

#For manual use: 
#treeToConvertLocation = "Data/CVHRemakeBinaryForegroundTree.rds"
#inputTree = readRDS(treeToConvertLocation)

plotTreeHighlightBranches2 =function (tree, outgroup = NULL, hlspecies, hlcols = NULL, main = "", 
                                      useGG = FALSE, bgcol = "black") 
{
  if (is.null(hlcols)) {
    hlcols <- c(2:(length(hlspecies) + 1))
  }
  if (length(hlcols) < length(hlspecies)) {
    hlcols <- rep_len(hlcols, length(hlspecies))
  }
  if (!is.null(outgroup)) {
    outgroup <- outgroup[outgroup %in% tree$tip.label]
    if (length(outgroup) > 0) {
      if (is.numeric(hlspecies)) {
        dummytree <- tree
        dummytree$edge.length <- c(rep(1, nrow(dummytree$edge)))
        for (i in c(1:length(hlspecies))) {
          dummytree$edge.length[hlspecies] <- i + 1
        }
        dummyrooted <- root(dummytree, outgroup)
      }
      rooted <- root(tree, outgroup)
    }
    else {
      print("No members of requested outgroup found in tree; keeping unrooted.")
      rooted <- tree
      outgroup <- NULL
    }
  }
  else {
    rooted <- tree
  }
  if (useGG) {
    if (is.numeric(hlcols)) {
      hlcols <- rep_len("#0000ff", length(hlspecies))
    }
    hlspecies_named <- vector(mode = "character")
    if (is.numeric(hlspecies)) {
      for (i in 1:length(hlspecies)) {
        hlspecies_named[i] <- tree$tip.label[hlspecies[i]]
      }
    }
    else hlspecies_named <- hlspecies
    rooted2 <- rooted
    mm <- min(rooted2$edge.length[rooted2$edge.length > 0])
    rooted2$edge.length[rooted2$edge.length == 0] <- max(0.02, 
                                                         mm/20)
    tipCols <- vector(mode = "character")
    x <- 1
    for (i in 1:length(tree$tip.label)) {
      if (tree$tip.label[i] %in% hlspecies_named) {
        tipCols[i] <- hlcols[x]
        x <- x + 1
      }
      else tipCols[i] <- bgcol
    }
    nlabel <- rooted2$Nnode + length(rooted2$tip.label)
    edgeCols <- vector(mode = "character", length = nlabel)
    x <- 1
    for (i in 1:nlabel) {
      if (i < length(rooted2$tip.label)) {
        if (rooted2$tip.label[i] %in% hlspecies_named) {
          edgeCols[i] <- hlcols[x]
          x <- x + 1
        }
        else edgeCols[i] <- bgcol
      }
      else edgeCols[i] <- bgcol
    }
    plotobj = ggtree(rooted2, color = edgeCols)
    plotobj = plotobj + geom_tiplab(color = tipCols, geom = "text", 
                                    cex = 3) + labs(title = main)
    return(plotobj)
  }
  else {
    colMaster <- c(rep(bgcol, nrow(rooted$edge)))
    if (is.numeric(hlspecies)) {
      if (!is.null(outgroup)) {
        hlcols <- c(bgcol, hlcols)
        colMaster <- hlcols[dummyrooted$edge.length]
      }
      else {
        for (i in 1:length(hlspecies)) {
          colMaster[hlspecies[i]] <- hlcols[i]
        }
      }
    }
    else {
      wspmr <- rooted$tip.label[rooted$edge[, 2]]
      for (i in 1:length(hlspecies)) {
        colMaster[which(wspmr == hlspecies[i])] <- hlcols[i]
      }
    }
    termedge <- order(rooted$edge[, 2])[1:length(rooted$tip.label)]
    colMasterTip <- colMaster[termedge]
    rooted2 = rooted
    mm = min(rooted2$edge.length[rooted2$edge.length > 0])
    rooted2$edge.length[rooted2$edge.length == 0] = max(0.02, 
                                                        mm/20)
    plotobj = plot.phylo(rooted2, main = main, edge.color = colMaster, 
                         tip.color = colMasterTip, edge.width = 2, cex = 0.8)
    return(plotobj)
  }
}


ZoonomTreeNameToCommon = function(tree, plot = T, isForegroundTree = T, manualAnnotLocation = "Data/manualAnnotationsSheet.csv", hlcol = "blue", bgcol = "black"){
  
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
      plotTreeHighlightBranches2(readableTree, hlspecies = which(inputTree$edge.length == 1), hlcols = hlcol, bgcol = bgcol)
    }else{
      plotTree(inputTree)
    }
  }
  return(inputTree)
}


ZonomNameConvertVectorCommon = function(namesVector){
  names = namesVector                                                    #make a vector of the names
  manualAnnot = read.csv("Data/manualAnnotationsSheet.csv")                     #improt manual annots file
  for(i in 1:length(names)){                                                    #for each name: 
    currentName = names[i]                                                      #use the 'i'th name in the list
    currentRow = manualAnnot[manualAnnot$FaName %in% currentName, ]             #find a row with the zonom name that matches the current name 
    currentSize = dim(currentRow)                                               #Part of "does row exist check": get the dimensions of the currentRow dataframe
    obsNumber = currentSize[1]                                                  #set "size" equal to the number of observations in 'currentRow'; which is the number of matches to the current name. If none exist it will be 0, if more than one it will be greater than 1. 
    if(obsNumber == 1){                                                         #if only one match exists:
      currentName = currentRow$Common.Name.or.Group                             #get the name from that row
    }else{
      currentName = names[i]                                                    #Otherwise, keep the name the same
    }
    names[i] = currentName                                                      #update the main name list with the name chose
  }
  #colnames(nMatrix) = names                                                     #update the matrix with the new names. 
  return(names)
}

# this script outputs the foreground edges for an "all" clades foreground given a tree and a foreground species vector. 

getForegroundEdges = function(inputTree, foregroundVector, plot = F){
  startnodes = min(inputTree$edge[,1]):max(inputTree$edge[,1])
  edges = inputTree$edge
  foregroundTips = which(inputTree$tip.label %in% foregroundVector)
  backgroundTips = which(!1:length(inputTree$tip.label) %in% foregroundTips)
  unsureNodes = startnodes
  for(i in 1:10){
    for(i in unsureNodes){
      endNodes = edges[which(edges[,1] == i),2]
      if(all(endNodes %in% foregroundTips)){
        foregroundTips = append(foregroundTips, i)
        unsureNodes = unsureNodes[unsureNodes != i]
      }else if(any(endNodes %in% backgroundTips)){
        backgroundTips = append(backgroundTips, i)
        unsureNodes = unsureNodes[unsureNodes != i]
      }else{
        unsureNodes = append(unsureNodes, i)
      }
      #message(endNodes)
    }
  }
  foregroundEdges = which(edges[,2] %in% foregroundTips)
  if(plot){
    plotTreeHighlightBranches(inputTree, hlspecies = foregroundEdges, hlcols = "blue")
  }
  return(foregroundEdges)
}