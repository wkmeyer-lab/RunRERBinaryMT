a = b #this is to prevent accidental full runs

library(RERconverge)

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
write.csv(report, file= "Results/geneTreesReport.csv")

numSpecies = rowSums(report)[order(rowSums(report), decreasing = T)]
test = table(numSpecies)

topGeneNames = names(numSpecies[1:length(which(numSpecies >453))])
topGeneNames = names(numSpecies)

topGenes = report[which(row.names(report)%in% topGeneNames),]

speciesInTopTrees = colSums(topGenes)[order(colSums(topGenes))]
speciesMissingFromTopTrees = speciesInTopTrees[speciesInTopTrees < length(topGeneNames)]
length(speciesMissingFromTopTrees)


# ---
testDrop = topGenes
i=1
tipsToDrop = character()

while(T){
  currentLowestSpecies = names(speciesMissingFromTopTrees[i])
  
  message(" -------------- ")
  message(" i = ", i )
  message(currentLowestSpecies)
  tipsToDrop = append(tipsToDrop, currentLowestSpecies)
  
  dropCol = which(colnames(testDrop) == currentLowestSpecies)
  
  testDrop = testDrop[,-dropCol]

  
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
tipsInNewMaster3 = colnames(testDrop)


tipsInNewMaster2 = colnames(testDrop)

tipsInNewMaster = colnames(testDrop)

all.equal(tipsInNewMaster3, tipsInNewMaster2)

saveRDS(tipsInNewMaster2, "Results/newZoMasterTips.rds")


length(tipsInNewMaster3)
length(tipsInNewMaster)


togaTree = read.tree("Data/togaTree.nwk")
plot.phylo(togaTree)

tipsToDrop = togaTree$tip.label[!togaTree$tip.label %in% tipsInNewMaster2]

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
