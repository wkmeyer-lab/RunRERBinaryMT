a = b #this is to prevent accidental full runs

palette(c("yellowgreen", "darkgray", "yellow", "darkgreen", "darkblue", "lightblue", "gold", "black", "pink", "red"))
palette(c("yellowgreen", "yellow", "darkgreen", "darkblue", "lightblue", "gold", "black", "pink", "red"))


palette(c("yellow", "darkgreen", "darkblue", "lightblue", "black", "pink", "red"))
palette(c( "darkgreen", "darkblue", "lightblue", "black", "red"))
palette(c( "darkgreen", "darkblue", "lightblue", "black", "pink", "red"))
palette(c( "darkgreen", "darkblue", "black", "red"))
palette(c(  "red", "darkgreen", "black"))


palette(c( "darkgreen", "black", "darkblue", "red"))


library(RERconverge)
# ---------------------------------------
length(phenotypeVector)



#-------------------------------------------
manualAnnotsTrimmed = manualAnnots
which(manualAnnots$ZoonomiaTip %in% names(phenotypeVector))
manualAnnotsTrimmed = manualAnnotsTrimmed[which(manualAnnots$ZoonomiaTip %in% names(phenotypeVector)), ]

table(manualAnnotsTrimmed$MSWC_Family)

# ----------------------------
lowCategoryGeneDropper(mainTrees, phenotypeVector)

# ----------------------------
categoricalCorrelation = correlateWithCategoricalPhenotype(RERObject, pathsObject, min.sp = 400, min.pos = 2) #Calculate with categorical, min 2 species per category 
overalCategorical = categoricalCorrelation[[1]]                               #select the results relating to overall difference between all categories
correlation = overalCategorical                                               # and classify it as the main correlation file

#process the pairwise outputs
pairwiseCategorical = categoricalCorrelation[[2]]                             #select the group of pairwise comparisons

phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #select the phenotype vector based on prefix
phenotypeVector = readRDS(phenotypeVectorFilename)                            #load in the phenotype vector 
categories = map_to_state_space(phenotypeVector)                              #and use it to connect branch lengths to phenotype name
categoryNames = categories$name2index                                         #store the length-phenotype connection

pairwiseTableNames = names(pairwiseCategorical)                               #Prepare to repalce the number-number titles with phenotype-phenotype titles
for(i in 1:length(categoryNames)){                                            #for each phenotype
  pairwiseTableNames= gsub(i, names(categoryNames)[i], pairwiseTableNames)                        #replace the number with the phenotype name  
}
names(pairwiseCategorical) = pairwiseTableNames                               #update the dataframe titles

pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "PairwiseCorrelationFile", sep= "") #make a name for the pairwise comparisons based on prefix
write.csv(pairwiseCategorical, file= paste(pairwiseCorrelationFileName, ".csv", sep=""), row.names = T, quote = F) #save the correlations as a csv
saveRDS(pairwiseCategorical, paste(pairwiseCorrelationFileName, ".rds", sep="")) #and as an rds 

combinedCategoricalCorrelationFilename = pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "CombinedCategoricalCorrelationFile", sep= "") # make this file for later functions that want it in combo
saveRDS(categoricalCorrelation, paste(combinedCategoricalCorrelationFilename, ".rds", sep="")) #and as an rds 

#save the outputs to subdirectories 
outputSubdirectoryNoslash = paste(outputFolderName, "Overall", sep = "")
if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
  dir.create(outputSubdirectoryNoslash)
}
outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")

correlationsOverallFilename = paste(outputSubdirectory, filePrefix, "OverallCorrelationFile.rds", sep= "")
saveRDS(categoricalCorrelation[[1]], correlationsOverallFilename)

for(i in 1:length(pairwiseTableNames)){
  pairwiseTableNames= gsub(" ", "", pairwiseTableNames)
  
  outputSubdirectoryNoslash = paste(outputFolderName, pairwiseTableNames[i], sep = "")
  if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
    dir.create(outputSubdirectoryNoslash)
  }
  outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
  
  correlationsPairFilename = paste(outputSubdirectory, filePrefix, pairwiseTableNames[i], "CorrelationFile",".rds", sep= "")
  saveRDS(categoricalCorrelation[[2]][[i]], correlationsPairFilename)
}





# -----------------------

report = mainTrees$report
view(report)

test = hist(rowSums(report))

length(which(rowSums(report)<400))

colnames(report) %in% names(phenotypeVector)
reportPruned = report[,colnames(report) %in% names(phenotypeVector)]

ncol(reportPruned)
colnames(reportPruned) %in% names(phenotypeVector)

test = hist(rowSums(reportPruned))
length(which(rowSums(reportPruned)<170))

?hist()



# -------------------------

mainTrees$masterTree$edge.length[1:length(mainTrees$masterTree$edge.length)] = 1

char2TreeCategoricalStates = function (tipvals, treesObj, useSpecies = NULL, model = "ER", 
          root_prior = "auto", plot = FALSE, anctrait = NULL) 
{
  mastertree = treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(mastertree$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      message(paste0("Species from master tree not present in useSpecies: ", 
                     paste(sp.miss, collapse = ",")))
    }
    useSpecies = intersect(mastertree$tip.label, useSpecies)
    mastertree = pruneTree(mastertree, useSpecies)
    mastertree = unroot(mastertree)
  }
  else {
    mastertree = pruneTree(mastertree, intersect(mastertree$tip.label, 
                                                 names(tipvals)))
    mastertree = unroot(mastertree)
  }
  if (is.null(anctrait)) {
    tipvals <- tipvals[mastertree$tip.label]
    intlabels <- map_to_state_space(tipvals)
    print("The integer labels corresponding to each category are:")
    print(intlabels$name2index)
    ancliks = getAncLiks(mastertree, intlabels$mapped_states, 
                         rate_model = model, root_prior = root_prior)
    states = rep(0, nrow(ancliks))
    for (i in 1:length(states)) {
      states[i] = which.max(ancliks[i, ])
    }
    states = c(intlabels$mapped_states, states)
    tree = mastertree
    tree$edge.length = states[tree$edge[, 2]]
    if (length(unique(tipvals)) == 2) {
      if (sum(!unique(tipvals) %in% c(TRUE, FALSE)) > 0) {
        message("Returning categorical tree for binary phenotype because phenotype values are not TRUE/FALSE")
      }
      else {
        tree$edge.length = ifelse(tree$edge.length == 
                                    2, 1, 0)
        print("There are only 2 categories: returning a binary phenotype tree.")
        if (plot) {
          plotTree(tree)
        }
        return(tree)
      }
    }
    if (plot) {
      plotTreeCategorical(tree, category_names = intlabels$state_names, 
                          master = mastertree, node_states = states)
    }
    return(states)
    return(tree)
  }
  else {
    if (length(unique(tipvals)) <= 2) {
      fgspecs <- names(tipvals)[tipvals != anctrait]
      res <- foreground2Tree(fgspecs, treesObj, plotTree = plot, 
                             clade = "terminal", useSpecies = useSpecies)
      print("There are only 2 categories: returning a binary phenotype tree.")
      if (plot) {
        plotTree(res)
      }
      return(res)
    }
    else {
      tipvals <- tipvals[mastertree$tip.label]
      intlabels <- map_to_state_space(tipvals)
      j <- which(intlabels$state_names == anctrait)
      if (length(j) < 1) {
        warning("The ancestral trait provided must match one of the traits in the phenotype vector.")
      }
      res = mastertree
      res$edge.length <- rep(j, length(res$edge.length))
      traits <- intlabels$state_names
      for (trait in traits) {
        if (trait == anctrait) {
          next
        }
        i <- which(intlabels$state_names == trait)
        res$edge.length[nameEdges(res) %in% names(tipvals)[tipvals == 
                                                             trait]] = i
      }
      names(res$edge.length) = nameEdges(res)
      if (plot) {
        states = res$edge.length[order(res$edge[, 2])]
        states = c(j, states)
        plotTreeCategorical(res, category_names = traits, 
                            master = treesObj$masterTree, node_states = states)
      }
      print("Category names are mapped to integers as follows:")
      print(intlabels$name2index)
      return(res)
    }
  }
}

commonStates = char2TreeCategoricalStates(commonPhenotypeVector, commonMainTrees, commonSpeciesFilter, model = modelType, anctrait = ancestralTrait, plot = T)
states = char2TreeCategoricalStates(phenotypeVector, mainTrees, speciesFilter, model = modelType, anctrait = ancestralTrait, plot = T) #use the phenotype vector to make a tree

categoricalTree
commonCategoricalTree

length(speciesFilter)
length(commonSpeciesFilter)

all.equal(categoricalTree$edge.length, commonCategoricalTree$edge.length)
all.equal(states, commonStates)

commonCategoricalTree$tip.label[which(duplicated(commonCategoricalTree$tip.label))]

manualAnnots$CommonName[which(duplicated(manualAnnots$CommonName))]

commonCategoricalTree = ZoonomTreeNameToCommon(categoricalTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
stableMaintrees = readRDS(mainTreesLocation)
stableCommonMainTrees = stableMaintrees
stableCommonMainTrees$masterTree = ZoonomTreeNameToCommon(stableCommonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)

?plotTreeCategorical
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableCommonMainTrees$masterTree)

plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableMaintrees$masterTree)



plotTreeCategorical(categoricalTree, c("Carnivore", "Herbivore", "Omnivore"), master = stableMaintrees$masterTree, node_states = states)

plotTreeCategorical(commonCategoricalTree, c("Carnivore", "Herbivore", "Omnivore"), master = stableCommonMainTrees$masterTree, node_states = states)



# --------------------------

phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #make a filename based on the prefix
phenotypeVector = readRDS(phenotypeVectorFilename)

lowCategoryGeneDropper = function(mainTrees, phenotypeVector){
  genesToDrop = vector()
  for(i in 1:length(mainTrees$trees)){
    currentTree = mainTrees$trees[[i]]
    currentTreeName = names(mainTrees$trees[i])
    currentTips = currentTree$tip.label
    
    phenotypedTips = currentTips[which(currentTips %in% names(phenotypeVector))]
    
    phenotypeValues = phenotypeVector[match(phenotypedTips, names(phenotypeVector))]
    phenotypeNumbers = table(phenotypeValues)
    if(any(phenotypeNumbers <3)){
      message(currentTreeName)
      print(phenotypeNumbers)
      genesToDrop = append(genesToDrop, currentTreeName)
    }
  }
  return(genesToDrop)
}  

lowCategoryGeneDropper(mainTrees, phenotypeVector)

# -------------------
ZoonomTreeNameToCommon(commonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
commonRERs = RERObject

colnames(commonRERs) = ZonomNameConvertVectorCommon(colnames(commonRERs), annotationLocation = spreadSheetLocation, tipCol = nameColumn)

plotRers(commonRERs, "ECI2", pathsObject)

?plotRers



# -----------------
testTree = readRDS("Output/CategoricalMobivoreTree/CategoricalMobivoreTreeCategoricalTree.rds")
table(testTree$edge.length)


REROne = readRDS("Output/CategoricalMobivoreTree/CategoricalMobivoreTreeRERFile-Other.rds")
RERTwo = readRDS("Output/CategoricalMobivoreTree/CategoricalMobivoreTreeRERFileOLd.rds")

all.equal(REROne, RERTwo)

demoTree = readRDS
mainTrees$masterTree$tip.label[mainTrees$masterTree$tip.label %in% "vs_OrnAna3"]

mainTrees$masterTree$edge.length[1:length(mainTrees$masterTree$edge.length)] = 1


testMain = readRDS("Data/zoonomiaAllMammalsTrees.rds")
testMain$master

grep("Ana", mainTrees$masterTree$tip.label)

phenotypeVector = readRDS(phenotypeVectorFilename)

categoricalTree$tip.label[!categoricalTree$tip.label %in% testTree$tip.label]
testTree$tip.label[!testTree$tip.label %in% categoricalTree$tip.label]
length(commonSpeciesFilter)

commonPhenotypeVector[names(commonPhenotypeVector) %in% "Platypus"]

prunedTree$edge.length

mainTrees2 = readRDS("Data/zoonomiaAllMammalsTrees.rds")
length(mainTrees2$masterTree$tip.label)

all.equal(mainTrees$masterTree, mainTrees2$masterTree)

length(mainTrees$masterTree$tip.label)
grep("ornAna", mainTrees$masterTree$tip.label)
mainTrees$masterTree$edge.length


testTrees = readRDS("data/RemadeTreesAllZoonomiaSpecies.rds")
testTrees$masterTree$edge.length

which(names(phenotypeVector) == "vs_HLlniGeo1")
phenotypeVector[94]

plotTreeCategorical2 = function (tree, category_names = NULL, master = NULL, node_states = NULL) 
{
  n = length(unique(tree$edge.length))
  if (n > length(palette())) {
    colors = colorRampPalette(palette())(n)
  }
  else {
    colors = c("yellowgreen", "lightgray", "yellow", "darkgreen", "darkblue", "lightblue", "gold", "black", "pink", "red")
  }
  edge_colors = tree$edge.length
  edge_colors = sapply(edge_colors, function(x) {
    colors[x]
  })
  par(mar = c(5, 4, 4, 10), xpd = TRUE)
  if (!is.null(master)) {
    cm = intersect(master$tip.label, tree$tip.label)
    master = pruneTree(master, cm)
    if (!is.null(node_states)) {
      node_colors = node_states
      node_colors = sapply(node_colors, function(x) {
        colors[x]
      })
      plot(master, cex = 0.25, edge.color = edge_colors, 
           node.color = node_colors)
    }
    else {
      plot(master, cex = 0.25, edge.color = edge_colors)
    }
  }
  else {
    if (!is.null(node_states)) {
      node_colors = node_states
      node_colors = sapply(node_colors, function(x) {
        colors[x]
      })
      plot(tree, cex = 0.25, edge.color = edge_colors, 
           use.edge.length = FALSE, node.depth = 2, node.color = node_colors)
    }
    else {
      plot(tree, cex = 0.25, edge.color = edge_colors, 
           use.edge.length = FALSE, node.depth = 2)
    }
  }
  if (!is.null(category_names)) {
    legend(x = "bottomright", inset = c(-0.25, 0), cex = 0.5, 
           legend = category_names, col = colors, lwd = 2)
  }
}

char2TreeCategorical2 = function (tipvals, treesObj, useSpecies = NULL, model = "ER", 
          root_prior = "auto", plot = FALSE, anctrait = NULL) 
{
  mastertree = treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(mastertree$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      message(paste0("Species from master tree not present in useSpecies: ", 
                     paste(sp.miss, collapse = ",")))
    }
    useSpecies = intersect(mastertree$tip.label, useSpecies)
    mastertree = pruneTree(mastertree, useSpecies)
    mastertree = unroot(mastertree)
  }
  else {
    mastertree = pruneTree(mastertree, intersect(mastertree$tip.label, 
                                                 names(tipvals)))
    mastertree = unroot(mastertree)
  }
  if (is.null(anctrait)) {
    tipvals <- tipvals[mastertree$tip.label]
    intlabels <- map_to_state_space(tipvals)
    print("The integer labels corresponding to each category are:")
    print(intlabels$name2index)
    ancliks = getAncLiks(mastertree, intlabels$mapped_states, 
                         rate_model = model, root_prior = root_prior)
    states = rep(0, nrow(ancliks))
    for (i in 1:length(states)) {
      states[i] = which.max(ancliks[i, ])
    }
    states = c(intlabels$mapped_states, states)
    tree = mastertree
    tree$edge.length = states[tree$edge[, 2]]
    if (length(unique(tipvals)) == 2) {
      if (sum(!unique(tipvals) %in% c(TRUE, FALSE)) > 0) {
        message("Returning categorical tree for binary phenotype because phenotype values are not TRUE/FALSE")
      }
      else {
        tree$edge.length = ifelse(tree$edge.length == 
                                    2, 1, 0)
        print("There are only 2 categories: returning a binary phenotype tree.")
        if (plot) {
          plotTree(tree)
        }
        return(tree)
      }
    }
    if (plot) {
      plotTreeCategorical2(tree, category_names = intlabels$state_names, 
                          master = mastertree, node_states = states)
    }
    return(tree)
  }
  else {
    if (length(unique(tipvals)) <= 2) {
      fgspecs <- names(tipvals)[tipvals != anctrait]
      res <- foreground2Tree(fgspecs, treesObj, plotTree = plot, 
                             clade = "terminal", useSpecies = useSpecies)
      print("There are only 2 categories: returning a binary phenotype tree.")
      if (plot) {
        plotTree(res)
      }
      return(res)
    }
    else {
      tipvals <- tipvals[mastertree$tip.label]
      intlabels <- map_to_state_space(tipvals)
      j <- which(intlabels$state_names == anctrait)
      if (length(j) < 1) {
        warning("The ancestral trait provided must match one of the traits in the phenotype vector.")
      }
      res = mastertree
      res$edge.length <- rep(j, length(res$edge.length))
      traits <- intlabels$state_names
      for (trait in traits) {
        if (trait == anctrait) {
          next
        }
        i <- which(intlabels$state_names == trait)
        res$edge.length[nameEdges(res) %in% names(tipvals)[tipvals == 
                                                             trait]] = i
      }
      names(res$edge.length) = nameEdges(res)
      if (plot) {
        states = res$edge.length[order(res$edge[, 2])]
        states = c(j, states)
        plotTreeCategorical(res, category_names = traits, 
                            master = treesObj$masterTree, node_states = states)
      }
      print("Category names are mapped to integers as follows:")
      print(intlabels$name2index)
      return(res)
    }
  }
}



treeImageFilename = paste(outputFolderName, filePrefix, "CategoricalTree.pdf", sep="") #make a filename based on the prefix
pdf(treeImageFilename, height = length(phenotypeVector)/18)                     #make a pdf to store the plot, sized based on tree size
char2TreeCategorical2(commonPhenotypeVector, commonMainTrees, commonSpeciesFilter, model = modelType, anctrait = ancestralTrait, plot = T)

categoricalTree = char2TreeCategorical2(phenotypeVector, mainTrees, speciesFilter, model = modelType, anctrait = ancestralTrait, plot = T) #use the phenotype vector to make a tree
dev.off()   



mainTrees3 = read.tree("Results/NewZoonomiaMasterTreePrunedToAlignmentSpecies.nwk")







# ------------------------------

fullTree = readRDS("Data/zoonomiaAllMammalsTrees.rds")
plotTree(commonMainTrees$masterTree)

fullTreeTrees = fullTree[[1]]

fullTree$masterTree

tipNumberList = numeric()
for(i in 1:length(fullTreeTrees)){
  currentTipNumber = length(fullTreeTrees[[i]]$tip.label)
  message(currentTipNumber)
  tipNumberList = append(tipNumberList, currentTipNumber)
}
length(tipNumberList)

max(tipNumberList)

tipList = list()
for(i in 1:length(fullTreeTrees)){
  currentTips = fullTreeTrees[[i]]$tip.label
  tipList = append(tipList, list(currentTips))
}

tipNumberList[order(tipNumberList, decreasing = T)]

plotTree(mainTrees$masterTree)
mainTrees$masterTree$tip.label
grep("ana", mainTrees$masterTree$tip.label)


biggestTree = fullTreeTrees[[8866]]$tip.label

singleMissing = fullTree$masterTree$tip.label[which(!fullTree$masterTree$tip.label %in% biggestTree)]

fullTreeTips = fullTree$masterTree$tip.label


tipTreeNumber = numeric()
for(i in 1:length(fullTreeTips)){
  currentTip = fullTreeTips[i]
  currentTipTreeNumber = length(which(sapply(tipList, function(x) currentTip %in% x)))
  names(currentTipTreeNumber) = currentTip
  message(currentTipTreeNumber)
  tipTreeNumber = append(tipTreeNumber, currentTipTreeNumber)
}

lowTipTrees = tipTreeNumber[order(tipTreeNumber)]



which(sapply(tipList, length) > 467)
highTipGenes = tipList[which(sapply(tipList, length) > 467)]

n=1
TestTip1 = names(lowTipTrees[1])
TestTip2 = names(lowTipTrees[2])


which(sapply(highTipGenes, function(x) !TestTip1 %in% x))
which(sapply(highTipGenes, function(x) !TestTip2 %in% x))

tipListTestDrop = tipList

lowTipTreesDropping = lowTipTrees[lowTipTrees < 10000]
tipsToDrop = names(lowTipTreesDropping) 
length(tipsToDrop)

tipListTestDrop = lapply(tipListTestDrop, function(x) Filter(function(y) !(y %in% tipsToDrop), x))


dropedLengths = sapply(tipListTestDrop, length)
dropedLengths[order(dropedLengths, decreasing = T)]

tressWith458 = which(dropedLengths == 458)

all10kspeciesTrees = fullTreeTrees[tressWith458]


all10kspeciesTreesSameTest = all10kspeciesTrees

all10kspeciesTreesSameTest$KAT7$edge.length = rep(1, length(all10kspeciesTreesSameTest$KAT7$edge.length))
all10kspeciesTreesSameTest$METTL1$edge.length = rep(1, length(all10kspeciesTreesSameTest$KAT7$edge.length))
all10kspeciesTreesSameTest$WNT2B$edge.length = rep(1, length(all10kspeciesTreesSameTest$KAT7$edge.length))
all10kspeciesTreesSameTest$TGFBI$edge.length = rep(1, length(all10kspeciesTreesSameTest$KAT7$edge.length))
all10kspeciesTreesSameTest$CFAP97D1$edge.length = rep(1, length(all10kspeciesTreesSameTest$KAT7$edge.length))


all.equal(all10kspeciesTreesSameTest$KAT7, all10kspeciesTreesSameTest$METTL1)
all.equal(all10kspeciesTreesSameTest$KAT7, all10kspeciesTreesSameTest$WNT2B)
all.equal(all10kspeciesTreesSameTest$KAT7, all10kspeciesTreesSameTest$TGFBI)
all.equal(all10kspeciesTreesSameTest$KAT7, all10kspeciesTreesSameTest$CFAP97D1)


report = fullTree$report
#write.csv(report, file= "Results/geneTreesReport.csv")


sum(report$vs_HLthyCyn1)

# ---

report = read.csv("Results/geneTreesReport.csv")
rownames(report) = report$X
report = report[,-1]

numSpecies = rowSums(report)[order(rowSums(report), decreasing = T)]
test = table(numSpecies)

topGeneNames = names(numSpecies[1:length(which(numSpecies >453))])
topGeneNames = names(numSpecies)

topGenes = report[which(row.names(report)%in% topGeneNames),]

speciesInTopTrees = colSums(topGenes)[order(colSums(topGenes))]
speciesMissingFromTopTrees = speciesInTopTrees[speciesInTopTrees < length(topGeneNames)]
length(speciesMissingFromTopTrees)


speciesToKeep = c("vs_HLornAna3", "vs_HLtacAcu1", "vs_HLgymLea1", "vs_HLpseCup1", "vs_ptePar1", "vs_HLpseCor1")


testDrop = topGenes
i=1
tipsToDrop = character()
while(T){
  currentLowestSpecies = names(speciesMissingFromTopTrees[i])
  
  message(" -------------- ")
  message(" i = ", i )
  message(currentLowestSpecies)
  if(!currentLowestSpecies %in% speciesToKeep){
  #if(T){
    tipsToDrop = append(tipsToDrop, currentLowestSpecies)
    dropCol = which(colnames(testDrop) == currentLowestSpecies)
    
    testDrop = testDrop[,-dropCol]
  }

  
  #check
  ncol(testDrop)
  print(rowSums(testDrop)[order(rowSums(testDrop))])
  
  numberOFFullTrees = length(which(rowSums(testDrop) == ncol(testDrop)))
  message(paste("number of matching trees =",numberOFFullTrees))
  
  if(numberOFFullTrees >9){
    message("Tips to drop:")
    print(tipsToDrop)
    break()
  }else(
    i = i+1
  )
}

tipsToDrop1 = tipsToDrop
tipsToDrop2 = tipsToDrop
tipsToDrop3 = tipsToDrop
tipsToDrop4 = tipsToDrop 

tipsToDropFinal = tipsToDrop


ZonomNameConvertVectorCommon(tipsToDrop, tipColumn = "Zoonomia")


tipsToDrop1 %in% tipsToDrop2
tipsToDrop4 %in% tipsToDrop1



tipsInNewMaster3 = colnames(testDrop)


tipsInNewMaster2 = colnames(testDrop)

tipsInNewMaster = colnames(testDrop)

all.equal(tipsInNewMaster3, tipsInNewMaster2)

saveRDS(tipsInNewMaster2, "Results/newZoMasterTips.rds")

tipsInUpdatedMaster = colnames(testDrop)


length(tipsInUpdatedMaster)
length(tipsInNewMaster)


togaTree = read.tree("Data/togaTree.nwk")
plot.phylo(togaTree)

tipsToDrop = togaTree$tip.label[!togaTree$tip.label %in% tipsInNewMaster2]

togaPruned = drop.tip(togaTree, tipsToDrop)
plot.phylo(togaPruned)

write.tree(togaPruned, "Results/NewZoonomiaMasterTreePrunedToAlignmentSpecies.nwk")

tipsToDrop = togaTree$tip.label[!togaTree$tip.label %in% tipsInUpdatedMaster]
togaPruned = drop.tip(togaTree, tipsToDrop)
plot.phylo(togaPruned)
write.tree(togaPruned, "Results/NewZoonomiaMasterTreePrunedToAlignmentSpecies.nwk")


#

test
?hist

dev.off()

which(is.na(commonMainTrees$masterTree$tip.label))
testPlot = plotTree(mainTrees$masterTree)

pdf(file = "test.pdf", width = 1000, height = 1000)
dev.off()
commonMainTrees$masterTree = commonMainTrees$masterTree

length(mainTrees$masterTree$edge.length)
which(is.na(mainTrees$masterTree$edge.length))

plot.phylo(mainTrees$masterTree)

mainTrees$masterTree$edge[which(is.na(mainTrees$masterTree$edge.length)),]

length(mainTrees$masterTree$tip.label)

NAFindTree = mainTrees$masterTree

NAFindTree$edge.length[which(is.na(NAFindTree$edge.length))] = 0.123456

plotTree(NAFindTree)

pdf("output/TestTree.pdf", height = length(NAFindTree$tip.label)/18)  

plotTreeHighlightBranches(NAFindTree, hlspecies = which(NAFindTree$edge.length == 0.123456), hlcols = "blue")

dev.off()

?plotTree


hillerMain = readRDS("data/NewHillerMainTrees.rds")
oldZoMain = readRDS("data/RemadeTreesAllZoonomiaSpecies.rds")

which(is.na(hillerMain$masterTree$edge.length))
which(is.na(oldZoMain$masterTree$edge.length))



#------ 

mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
masterTree = mainTrees$masterTree

masterTree$tip.label
