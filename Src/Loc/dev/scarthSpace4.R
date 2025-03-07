a = b #this is to prevent accidental full runs

palette(c("yellowgreen", "darkgray", "yellow", "darkgreen", "darkblue", "lightblue", "gold", "black", "pink", "red"))
palette(c("yellowgreen", "yellow", "darkgreen", "darkblue", "lightblue", "gold", "black", "pink", "red"))


palette(c("yellow", "darkgreen", "darkblue", "lightblue", "black", "pink", "red"))
palette(c( "darkgreen", "darkblue", "lightblue", "black", "red"))
palette(c( "darkgreen", "darkblue", "lightblue", "black", "pink", "red"))
palette(c( "darkgreen", "darkblue", "black", "red"))
palette(c(  "red", "darkgreen", "black"))


palette(c( "darkgreen", "black", "darkblue", "red"))
palette(c("black", "darkblue"))

palette(c( "darkgreen", "darkblue", "black", "red", "gray"))
palette(c( "darkgreen", "blue", "pink", "red"))
palette(c( "darkgreen", "blue", "purple", "red"))
palette(c( "gray", "darkblue", "darkgreen",  "black", "red"))

library(RERconverge)
source("Src/Reu/ZoonomTreeNameToCommon.R")
# ---------------------------------------
length(phenotypeVector)


?correlateWithCategoricalPhenotype
# -----------------------------------------
stableMaintrees = mainTrees
mainTrees = stableMaintrees
stableMaintrees = readRDS(mainTreesLocation)

mainTrees$masterTree$edge.length[1:length(mainTrees$masterTree$edge.length)] = 1

stableCommonMainTrees = stableMaintrees
stableCommonMainTrees$masterTree = ZoonomTreeNameToCommon(stableCommonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)

?plotTreeCategorical
pdf(treeImageFilename, height = length(phenotypeVector)/18, width = 10)     
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableCommonMainTrees$masterTree)

plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableMaintrees$masterTree)
dev.off()  


plotTreeCategorical(categoricalTree, c("Carnivore", "Herbivore", "Omnivore"), master = stableMaintrees$masterTree)

plotTreeCategorical(commonCategoricalTree, c("Carnivore", "Herbivore", "Omnivore"), master = stableCommonMainTrees$masterTree)

# --- make plots to demonstrate binary trees --- 

mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")
mainCategoricalTree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")

plotTreeCategorical(mainCategoricalTree, c("Carnivore", "Herbivore", "Omnivore"), master = stableMaintrees$masterTree)
stableCommonMainTrees = stableMaintrees
stableCommonMainTrees$masterTree = ZoonomTreeNameToCommon(stableCommonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)

commonCategoricalTree = ZoonomTreeNameToCommon(mainCategoricalTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
phenotypeVector = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPhenotypeVector.rds")
pdf("Results/BinaryDemoTrees.pdf", height = length(phenotypeVector)/18, width = 10)     
palette(c( "darkgreen", "gray", "gray", "gray"))
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableCommonMainTrees$masterTree)
palette(c( "gray", "darkblue", "gray", "gray"))
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableCommonMainTrees$masterTree)
palette(c( "gray", "gray", "black", "gray"))
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableCommonMainTrees$masterTree)
palette(c( "gray", "gray", "gray", "red"))
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableCommonMainTrees$masterTree)
dev.off()


# ------ plot RERs of CategoricalINsvertivore to check rho direction meaning 
library(RERconverge)
RERobject = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeRERFile.rds")
pathsObject = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPathsFile.rds")
mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
phenotypeSet = c("Herbivore", "Insectivore", "Omnivore", "Vertivore")
colorset = c( "darkgreen", "darkblue", "black", "red")

?plotRers
plotRers(RERobject, "BPIFB1", pathsObject)

source("Src/Reu/rerViolinPlot.R")
rerViolinPlot(mainTrees, RERobject, pathsObject, phenotypeSet , geneOfInterest = "BPIFB1", colorScale = colorset)
rerViolinPlot()

# -- extract genes from top results -- 


lines <- readLines("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/TopResults.txt")

# Filter lines that start with two tab characters
filtered_lines <- grep("^\\t\\t", lines, value = TRUE)

# Print or save the filtered lines
print(filtered_lines)

genesToRun = unique(filtered_lines)
genesToRun = gsub("\t", "", genesToRun)
HIgenesToRun = genesToRun

writeLines(genesToRun, "Results/hyphyGenesToRun.txt")

parseTopResultFileToGenes = function(file){
  lines <- readLines(file)
  
  # Filter lines that start with two tab characters
  filtered_lines <- grep("^\\t\\t", lines, value = TRUE)
  
  # Print or save the filtered lines
  #print(filtered_lines)
  
  genesToRun = unique(filtered_lines)
  genesToRun = gsub("\t", "", genesToRun)
  genesToRun
}

HItopResults = "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/TopResults.txt"
HVtopResults = "Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/TopResults.txt"
IVtopResults = "Output/CategoricalInsVertivoreTree/Insectivore-Vertivore/TopResults.txt"

testOut = parseTopResultFileToGenes(HItopResults)

HVgenesToRun = parseTopResultFileToGenes(HVtopResults)
IVgenesToRun = parseTopResultFileToGenes(IVtopResults)

VgenesToRun = append(HVgenesToRun, IVgenesToRun)
VgenesToRun[duplicated(VgenesToRun)]
VgenesToRun = unique(VgenesToRun)

VOonlyGenes = VgenesToRun[VgenesToRun %in% HIgenesToRun]
writeLines(VOonlyGenes, "Results/hyphyGenesToRunVO.txt")

HVOonlyGenes = VgenesToRun[!VgenesToRun %in% HIgenesToRun]
writeLines(HVOonlyGenes, "Results/hyphyGenesToRunHVO.txt")
#-----------------------------------

<<<<<<< HEAD
# -- exmaining the maturity results ---
??rer
rerTree = returnRersAsTree(mainTrees, RERObject, index = "PTCD1", phenv = pathsObject)
treePlotRers(mainTrees, RERObject, index = "PTCD1", phenv = pathsObject, type = "color")
=======
source("Src/Reu/treeColorPlots.R")

treeColorByLabel(phenMasterTree)
nodelabels(frame="none")


# --- making nexus trees of genes of interest ------

EHHADHTree = mainTrees$trees$EHHADH

write.nexus(EHHADHTree, "Results/EHHADHTree.nex")
write.tree(EHHADHTree, "Results/EHHADHTree.tree")
writeNexus(EHHADHTree, "Results/EHHADHTree.nex")
EHHADHTree$node.label = NULL

masterTree = mainTrees$masterTree
writeNexus(EHHADHTree, "Results/masterTree.nex")

?writeNexus


# ----- Trim fasta file to master tree ------

fasta = read_fasta("Results/ENST00000231887.EHHADH.filt.fa")

fastaTipHeaders = fasta$headers
fastaTipHeaders = sub("\\t.*", "", fastaTipHeaders)
fastaTipHeaders[fastaTipHeaders == "REFERENCE"] = "vs_hg38"

which(!fastaTipHeaders %in% masterTree$tip.label)

fasta$headers = fasta$headers[-which(!fastaTipHeaders %in% masterTree$tip.label)]
fasta$sequences = fasta$sequences[-which(!fastaTipHeaders %in% masterTree$tip.label)]

read.dna("Results/ENST00000231887.EHHADH.filt.fa")
# -------- Redo with new package ---- 


??fasta
library(seqinr)
fastaLocation = "Results/ENST00000231887.EHHADH.filt.fa"
mainTreesLocation = 'data/zoonomiaAllMammalsTrees.rds'

fasta = read.fasta(fastaLocation)

fastaTipHeaders = names(fasta)
fastaTipHeaders = sub("\\t.*", "", fastaTipHeaders)
fastaTipHeaders[fastaTipHeaders == "REFERENCE"] = "vs_hg38"

if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
masterTree = mainTrees$masterTree
noDataTips = masterTree$tip.label[!masterTree$tip.label %in% fastaTipHeaders]
masterTree = drop.tip(masterTree, noDataTips)
masterTree$node.label = NULL


fastaToDrop = which(!fastaTipHeaders %in% masterTree$tip.label)

fasta = fasta[-fastaToDrop]
fastaTipHeaders = fastaTipHeaders[-fastaToDrop]

names(fasta) = fastaTipHeaders

write.fasta(fasta, names = names(fasta), file.out = "Results/outFasta.fa")
writeNexus(masterTree, "Results/masterTreeOut.nex")
write.tree(masterTree, "Results/masterTreeOut.tree")


masterTree$tip.label[order(masterTree$tip.label)]
fastaTipHeaders[order(fastaTipHeaders)]

#
masterTree2 = mainTrees$masterTree
masterTree2$tip.label
masterTree2$node.label = NULL

write.tree(masterTree2, "Results/fullMasterTreeOut.tree")


phenotypeTree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")
masterTree2$tip.label
phenotypeTree$tip.label

phenMasterTree = masterTree2
phenMasterTree = drop.tip(phenMasterTree, phenMasterTree$tip.label[!phenMasterTree$tip.label %in% phenotypeTree$tip.label])

# ---- Remake file creation code from start cleanly -------
library(seqinr)
fastaLocation = "Results/ENST00000231887.EHHADH.filt.fa"
mainTreesLocation = 'data/zoonomiaAllMammalsTrees.rds'
foregroundCategory = "1"

useManualTree = F
filePrefix = "CategoricalInsVertivoreTree"
index = "EHHADH"
phenotypeTreeLocation = "Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds"




if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
masterTree = mainTrees$masterTree
masterTree$node.label = NULL

if(useManualTree){
  phenotypeTree = readRDS(phenotypeTreeLocation)
}else{
  source("Src/Reu/paths2Tree.R")
  outputFolderName = paste("Output/",filePrefix,"/", sep = "")
  pathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "") #make a filename based on the prefix
  
  pathsObject = readRDS(pathsFilename)
  pathsTree = paths2Tree(mainTrees, pathsObject, index)
  #paths tree actually currently being unused because of how the master tree phenotype matching works. Because it's relying on the trees beingthe same shape and therefore having matching node numbers, I can't use the paths -- or, at least, it's very messy to try, so I'm not.
  
  
  phenotypeTreeCategoricalLocation = paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep="") #make a filename based on the prefix
  phenotypeTreeBinaryLocation = paste(outputFolderName, filePrefix, "BinaryTree.rds", sep="") #make a filename based on the prefix
  if(file.exists(phenotypeTreeCategoricalLocation)){
    fullPhenotypeTree = readRDS(phenotypeTreeCategoricalLocation)
  }else if(file.exists(phenotypeTreeBinaryLocation)){
    fullPhenotypeTree = readRDS(phenotypeTreeBinaryLocation)
  }else{
    stop("The prefix has neither a categorical or binary phenotype tree")
  }
    
}
phenotypeTree = fullPhenotypeTree


phenMasterTree = masterTree
phenMasterTree = drop.tip(phenMasterTree, phenMasterTree$tip.label[!phenMasterTree$tip.label %in% fullPhenotypeTree$tip.label])
phenMasterTree = drop.tip(phenMasterTree, phenMasterTree$tip.label[!phenMasterTree$tip.label %in% phenotypeTree$tip.label])

#add category as label to nodes
allLabels = rep("", (length(phenMasterTree$tip.label)+phenMasterTree$Nnode))
for(i in 1:length(allLabels)){
  message(i)
  parentEdge = which(phenMasterTree$edge[,2]==i)
  if(!length(parentEdge)==0){
    allLabels[i] = parentEdge
    allLabels[i] = phenotypeTree$edge.length[parentEdge]
    }
}

allLabels[which(allLabels == foregroundCategory)] = "Foreground"

tipLabels = allLabels[c(1:length(phenotypeTree$tip.label))]
originalTipValues = phenMasterTree$tip.label
phenMasterTree$tip.label = paste0(phenMasterTree$tip.label, "{", tipLabels, "}")
internalLabels = allLabels[-c(1:length(phenotypeTree$tip.label))]
phenMasterTree$node.label = paste0("{", internalLabels, "}")



# - Read Fasta file - 
fasta = read.fasta(fastaLocation)

fastaTipHeaders = names(fasta)
fastaTipHeaders = sub("\\t.*", "", fastaTipHeaders)
fastaTipHeaders[fastaTipHeaders == "REFERENCE"] = "vs_hg38"

# trim files to match eachother 
noDataTips = phenMasterTree$tip.label[!originalTipValues %in% fastaTipHeaders]
phenMasterTree = drop.tip(phenMasterTree, noDataTips)

fastaToDrop = which(!fastaTipHeaders %in% originalTipValues)
fasta = fasta[-fastaToDrop]
fastaTipHeaders = fastaTipHeaders[-fastaToDrop]
names(fasta) = fastaTipHeaders

#Write output file
write.fasta(fasta, names = names(fasta), file.out = "Results/outFasta.fa")
write.tree(phenMasterTree, "Results/masterTreeOut.tree")

fastaOut = readLines("Results/outFasta.fa")
treeOut = readLines("Results/masterTreeOut.tree")
combinedContent = c(fastaOut, treeOut)
writeLines(combinedContent, "Results/hyphyOutput.fna")

# ------------------------------------------
treeColorPlot(phenMasterTree)
treeColorByLabel(phenMasterTree)
prunedPaths = drop.tip(pathsTree, pathsTree$tip.label[!pathsTree$tip.label %in% substr(phenMasterTree$tip.label, 1, nchar(phenMasterTree$tip.label)-3)])

treeColorByLabel(prunedPaths)
treeColorPlot(prunedPaths)
plot.phylo(prunedPaths)
prunedPaths$edge.length

fasta = fasta[1:74]

#




#
?correlateWithContinuousPhenotype

pathsVector = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPathsFile.rds")
??paths
tree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")
treesObj = mainTrees
categorical = T
tree2Paths
useSpecies = NULL
source("Src/Reu/RERConvergeFunctions.R")



?RERconverge
?plotTreeCategorical()
?treePlotNew
?plotTreeHighlightBranches()


pathsVector
length(pathsVector)
rerObject = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeRERFile.rds")
ncol(rerObject)

rerObjectTrimmed = (rerObject[1:2,])

test = rbind(pathsVector, pathsVector)

returnRersAsTree(mainTrees, test, 2)
returnRersAsTree(mainTrees, pathsVector, 1)

for(i in 1:10){
 outname = paste0("path", i)
 outPath = paths2Tree(mainTrees, pathsVector, i)
 assign(outname, outPath)
}
all.equal(path1, path2)

# -------------------------------

categoricalTree = readRDS()

??paths
??tree

rawMaturityVector = readRDS("Output/MaturityLogRaw/MaturityLogRawContinuousPhenotypeVector.rds")
percentMaturityVector = readRDS("Output/MaturityLifespanPercent/MaturityLifespanPercentContinuousPhenotypeVector.rds")
rawPaths = readRDS("Output/MaturityLogRaw/MaturityLogRawContinuousPathsFile.rds")
percentPaths = readRDS("Output/MaturityLifespanPercent/MaturityLifespanPercentContinuousPathsFile.rds")

percentRERs = readRDS("Output/MaturityLifespanPercent/MaturityLifespanPercentRERFile.rds")
rawMatRers = readRDS("Output/MaturityLogRaw/MaturityLogRawRERFile.rds")
all.equal(percentRERs, rawMatRers)



commonRERs = RERObject
colnames(commonRERs) = ZonomNameConvertVectorCommon(colnames(commonRERs), annotationLocation = spreadSheetLocation, tipCol = nameColumn)

targetGene = "PTCD1"
targetGene = "ACOT8"

targetGene = "RPS8"
targetGene = "RPS3"
targetGene = "RPS17"
targetGene = "RPL27"

targetGene = "ZIC1"
targetGene = "GTSF1"
targetGene = "OR51E1"

which(colnames(commonRERs) == "Chinese River Dolphin")


testMat = commonRERs[1:2, 1:3]

testDropRER = commonRERs[,-which(colnames(commonRERs) %in% c("Chinese River Dolphin", "Star-nosed mole", "Beaver"))]
testDropPaths = pathsObject[-which(colnames(commonRERs) %in% c("Chinese River Dolphin", "Star-nosed mole","Beaver"))]

testDropRER = commonRERs[,-which(colnames(commonRERs) %in% c("Chinese River Dolphin", "Star-nosed mole"))]
testDropPaths = pathsObject[-which(colnames(commonRERs) %in% c("Chinese River Dolphin", "Star-nosed mole"))]


commonRERs[which(rownames(commonRERs) == "GTSF1"),][order(commonRERs[which(rownames(commonRERs) == "GTSF1"),])]
{
x=pathsObject
y=commonRERs[targetGene,]
names(y)==namePathsWSpecies(mainTrees$masterTree)

plot(x,y, cex.axis=1, cex.lab=1, cex.main=1, xlab="Maturity Percentage Change",
     ylab="Evolutionary Rate", main=paste("Gene",targetGene,"Pearson Correlation"),
     pch=19, cex=1)
text(x,y, labels=names(y), pos=4)
abline(lm(y~x), col='red',lwd=3)
} # plot with percentage labels

{
  x=pathsObject
  y=commonRERs[targetGene,]
  names(y)==namePathsWSpecies(mainTrees$masterTree)
  
  plot(x,y, cex.axis=1, cex.lab=1, cex.main=1, xlab="Maturity Raw Change",
       ylab="Evolutionary Rate", main=paste("Gene",targetGene,"Pearson Correlation"),
       pch=19, cex=1, ylim=c(-1, 3))
  text(x,y, labels=names(y), pos=4)
  abline(lm(y~x), col='red',lwd=3)
} #plot with Raw labels
m <- lm(y ~ x)


{
  x=testDropPaths
  y=testDropRER[targetGene,]
  names(y)==namePathsWSpecies(mainTrees$masterTree)
  
  plot(x,y, cex.axis=1, cex.lab=1, cex.main=1, xlab="Maturity Raw Change",
       ylab="Evolutionary Rate", main=paste("Gene",targetGene,"Pearson Correlation"),
       pch=19, cex=1, ylim=c(-1, 3))
  text(x,y, labels=names(y), pos=4)
  abline(lm(y~x), col='red',lwd=3)
} #plot with Raw labels
m <- lm(y ~ x)

eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                 list(a = format(unname(coef(m)[1]), digits = 2),
                      b = format(unname(coef(m)[2]), digits = 2),
                      r2 = format(summary(m)$r.squared, digits = 3)))
as.character(as.expression(eq));

{
x=rawPaths
y=percentPaths
names(y)==namePathsWSpecies(mainTrees$masterTree)

plot(x,y, cex.axis=1, cex.lab=1, cex.main=1, xlab="Lifespan Percent Change",
     ylab="Log Raw change", main=paste("Comparission of percentage and raw change"),
     pch=19, cex=1)
text(x,y, labels=names(y), pos=4)
abline(lm(y~x), col='red',lwd=3)
}
# ---- Compare the ranking of the downsampled and non-downsampled results ---- 

mainIhCorrels = readRDS("Output/CategoricalInsvertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreCorrelationFile.rds")
mainVhCorrels = readRDS("Output/CategoricalInsvertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreCorrelationFile.rds")

downIhCorrels = readRDS("Output/CategoricalDownsampledInsvertTree/Herbivore-Insectivore/CategoricalDownsampledInsVertTreeHerbivore-InsectivoreCorrelationFile.rds")
downVhCorrels = readRDS("Output/CategoricalDownsampledInsvertTree/Herbivore-Vertivore/CategoricalDownsampledInsVertTreeHerbivore-VertivoreCorrelationFile.rds")

halfIhCorrels = readRDS("Output/CategoricalNoMegabranchInsvertTree/Herbivore-Insectivore/CategoricalNoMegabranchInsvertTreeHerbivore-InsectivoreCorrelationFile.rds")
halfVhCorrels = readRDS("Output/CategoricalNoMegabranchInsvertTree/Herbivore-Vertivore/CategoricalNoMegabranchInsvertTreeHerbivore-VertivoreCorrelationFile.rds")


mainIhCorrels$order = rep(1:nrow(mainIhCorrels))
mainVhCorrels$order = rep(1:nrow(mainVhCorrels))
downIhCorrels$order = rep(1:nrow(downIhCorrels))
downVhCorrels$order = rep(1:nrow(downVhCorrels))
halfIhCorrels$order = rep(1:nrow(halfIhCorrels))
halfVhCorrels$order = rep(1:nrow(halfVhCorrels))


mainIhCorrels = mainIhCorrels[order(mainIhCorrels$p.adj),]
mainVhCorrels = mainVhCorrels[order(mainVhCorrels$p.adj),]
downIhCorrels = downIhCorrels[order(downIhCorrels$p.adj),]
downVhCorrels = downVhCorrels[order(downVhCorrels$p.adj),]
halfIhCorrels = halfIhCorrels[order(halfIhCorrels$p.adj),]
halfVhCorrels = halfVhCorrels[order(halfVhCorrels$p.adj),]

mainIhCorrels$rank = rep(1:nrow(mainIhCorrels))
mainVhCorrels$rank = rep(1:nrow(mainVhCorrels))
downIhCorrels$rank = rep(1:nrow(downIhCorrels))
downVhCorrels$rank = rep(1:nrow(downVhCorrels))
halfIhCorrels$rank = rep(1:nrow(halfIhCorrels))
halfVhCorrels$rank = rep(1:nrow(halfVhCorrels))

mainIhCorrels = mainIhCorrels[order(mainIhCorrels$order),]
colnames(mainIhCorrels) = paste0("mih", colnames(mainIhCorrels))
mainVhCorrels = mainVhCorrels[order(mainVhCorrels$order),]
colnames(mainVhCorrels) = paste0("mvh", colnames(mainVhCorrels))
downIhCorrels = downIhCorrels[order(downIhCorrels$order),]
colnames(downIhCorrels) = paste0("dih", colnames(downIhCorrels))
downVhCorrels = downVhCorrels[order(downVhCorrels$order),]
colnames(downVhCorrels) = paste0("dvh", colnames(downVhCorrels))
halfIhCorrels = halfIhCorrels[order(halfIhCorrels$order),]
colnames(halfIhCorrels) = paste0("hih", colnames(halfIhCorrels))
halfVhCorrels = halfVhCorrels[order(halfVhCorrels$order),]
colnames(halfVhCorrels) = paste0("hvh", colnames(halfVhCorrels))

combinedData = cbind(mainIhCorrels, mainVhCorrels, downIhCorrels, downVhCorrels, halfIhCorrels, halfVhCorrels)


plot(combinedData$mihrank, combinedData$dihrank)
plot(combinedData$mvhrank, combinedData$dvhrank)
plot(combinedData$hvhrank, combinedData$dvhrank)
plot(combinedData$mvhrank, combinedData$hvhrank)

plot(combinedData$mihp.adj, combinedData$dihp.adj)
plot(combinedData$mvhp.adj, combinedData$dvhp.adj)
plot(combinedData$hvhp.adj, combinedData$dvhp.adj)


# --- Check downsampling RERS are equal ---- 

mainRERs = readRDS("Output/CategoricalInsvertivoreTree/CategoricalInsVertivoreTreeRERFile.rds")
halfDownRERs = readRDS("Output/CategoricalNoMegabranchInsvertTree/CategoricalNoMegabranchInsvertTreeRERFile.rds")
downRERs = readRDS("Output/CategoricalDownsampledInsvertTree/CategoricalDownsampledInsvertTreeRERFile.rds")

all.equal(mainRERs, halfDownRERs)
all.equal(downRERs, halfDownRERs)

mainFilter = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeSpeciesFilterJanThirty.rds")
halfdownFilter = readRDS("Output/CategoricalNoMegabranchInsvertTree/CategoricalNoMegabranchInsvertTreeSpeciesFilter.rds")
downFilter = readRDS("Output/CategoricalDownsampledInsvertTree/CategoricalDownsampledInsvertTreeSpeciesFilter.rds")

mainFilterOld = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeSpeciesFilter.rds")

mainFilter
halfdownFilter
downFilter
all.equal(halfdownFilter, downFilter)

all.equal(mainFilter, downFilter)
all.equal(mainFilter, mainFilterOld)

# ---- Make phylogeneticaly matching downsampled trees ----
insVertTree = readRDS("Output/CategoricalInsvertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")
insVertPhenv = readRDS("Output/CategoricalInsvertivoreTree/CategoricalInsVertivoreTreeCategoricalPhenotypeVector.rds")
table(insVertTree$edge.length)
table(insVertPhenv)

pdf(height = 18, width = 10)  
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableCommonMainTrees$masterTree)
edgelabels(col = "darkgreen", frame = "none")
dev.off()

noMegaBranchesTree = insVertTree
noMegaBranchesTree$edge.length[c(1,2,367)] = NA

categoricalTree$edge.length[c(1,2)] = 5
commonCategoricalTree$edge.length[c(1,2)] = 5

pathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "") #make a filename based on the prefix
paths = tree2Paths(categoricalTree, mainTrees, useSpecies = speciesFilter, categorical = T)
#char2PathsCategorical(phenotypeVector, mainTrees, speciesFilter, model = modelType, anctrait = ancestralTrait) #make a path based on the phenotype vector
saveRDS(paths, file = pathsFilename)

?tree2Paths
categoricalTree$edge.length[c(389,388, 370, 366, 365, 359, 351, 352, 350, 347, 339, 337, 110, 336, 335, 334)] = 5
categoricalTree$edge.length[c(1,2, 364, 360, 358, 342, 343, 344, 345, 346, 348, 329, 330, 331, 332, 333, 338, 145, 3, 147)] = 5



commonCategoricalTree$edge.length[c(389,388, 370, 366, 365, 359, 351, 352, 350, 347, 339, 337, 110, 336, 335, 334)] = 5
commonCategoricalTree$edge.length[c(1,2, 364, 360, 358, 342, 343, 344, 345, 346, 348, 329, 330, 331, 332, 333, 338, 145, 3, 147)] = 5



#-------------------

maturityRER = readRDS("Output/MaturityLifespanPercent/MaturityLifespanPercentRERFile.rds")
maturityPath = readRDS("Output/MaturityLifespanPercent/MaturityLifespanPercentContinuousPathsFile.rds")
maturityMaintrees= readRDS("data/newHillerMainTrees.rds")

commonMaturityRER = maturityRER
colnames(commonMaturityRER)
colnames(commonMaturityRER) = ZonomNameConvertVectorCommon(colnames(commonMaturityRER), tipColumn = "ZoonomiaName")

returnRersAsTree(maturityMaintrees, maturityRER, "NDRG4", maturityPath)
plotRers(commonMaturityRER, "NDRG4", phenv = maturityPath)



#-----------------------------------

InsectivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-GO_Biological_Process_2023.rds")
write.csv(InsectivoreGoData, "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-GO_Biological_Process_2023.csv")

?correlateWithBinaryPhenotype
install.packages()

#----------------------------------------

InsectivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-GO_Biological_Process_2023.rds")
VertivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreEnrichment-GO_Biological_Process_2023.rds")
CarnivoreGoData = readRDS("Output/CategoricalPrunedCarnivoryTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreEnrichment-GO_Biological_Process_2023.rds")


#----------------------------------------
library(ggvenn)
InsectivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-GO_Biological_Process_2023.rds")[[1]]
VertivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreEnrichment-GO_Biological_Process_2023.rds")[[1]]
CarnivoreGoData = readRDS("Output/CategoricalPrunedCarnivoreTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreEnrichment-GO_Biological_Process_2023.rds")[[1]]

InsectivoreGoData = InsectivoreGoData[order(InsectivoreGoData$p.adj),]
VertivoreGoData = VertivoreGoData[order(VertivoreGoData$p.adj),]
CarnivoreGoData = CarnivoreGoData[order(CarnivoreGoData$p.adj),]


signficiantInsectivore = InsectivoreGoData[which(InsectivoreGoData$pval <0.05),]
signficiantVertivore = VertivoreGoData[which(VertivoreGoData$pval <0.05),]
signficiantCarnivore = CarnivoreGoData[which(CarnivoreGoData$pval <0.05),]

vennData = list(
  Insectivore = rownames(signficiantInsectivore),
  Vertivore = rownames(signficiantVertivore),
  Carnivore = rownames(signficiantCarnivore)
)

signficianterInsectivore = InsectivoreGoData[which(InsectivoreGoData$p.adj <0.05),]
signficianterVertivore = VertivoreGoData[which(VertivoreGoData$p.adj <0.05),]
signficianterCarnivore = CarnivoreGoData[which(CarnivoreGoData$p.adj <0.05),]

topInsectivore = InsectivoreGoData[1:100,]
topVertivore = VertivoreGoData[1:100,]
topCarnivore = CarnivoreGoData[1:100,]

vennData = list(
  Insectivore = rownames(topInsectivore),
  Vertivore = rownames(topVertivore),
  Carnivore = rownames(topCarnivore)
)
ggvenn(vennData, fill_color = c("blue", "red", "pink"))


# --------

InsectivoreGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreCorrelationFile.rds")
VertivoreGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreCorrelationFile.rds")
CarnivoreGeneData = readRDS("Output/CategoricalPrunedCarnivoreTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreCorrelationFile.rds")

InsectivoreGeneData = InsectivoreGeneData[order(InsectivoreGeneData$p.adj),]
VertivoreGeneData = VertivoreGeneData[order(VertivoreGeneData$p.adj),]
CarnivoreGeneData = CarnivoreGeneData[order(CarnivoreGeneData$p.adj),]

sigInsectGenes = InsectivoreGeneData[which(InsectivoreGeneData$p.adj <0.05),]
sigVertGenes = VertivoreGeneData[which(VertivoreGeneData$p.adj <0.05),]
sigCarnGenes = CarnivoreGeneData[which(CarnivoreGeneData$p.adj <0.05),]

geneVennData = list(
  Insectivore = rownames(sigInsectGenes),
  Vertivore = rownames(sigVertGenes),
  Carnivore = rownames(sigCarnGenes)
)
ggvenn(geneVennData, fill_color = c("blue", "red", "pink"))

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
stableMaintrees = mainTrees
mainTrees$masterTree$edge.length[1:length(mainTrees$masterTree$edge.length)] = 1
stableMaintrees = readRDS(mainTreesLocation)
stableCommonMainTrees = stableMaintrees
stableCommonMainTrees$masterTree = ZoonomTreeNameToCommon(stableCommonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)

?plotTreeCategorical
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableCommonMainTrees$masterTree)

plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableMaintrees$masterTree)



plotTreeCategorical(categoricalTree, c("Carnivore", "Herbivore", "Omnivore"), master = stableMaintrees$masterTree)

plotTreeCategorical(commonCategoricalTree, c("Carnivore", "Herbivore", "Omnivore"), master = stableCommonMainTrees$masterTree)



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
