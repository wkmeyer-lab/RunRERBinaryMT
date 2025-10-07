a = b #prevent full runs
library(RERconverge)
library(tools)


# ------------------------------------------------------------------
# --- Looking into differecen between liam and non-liam results  ----- 
# ------------------------------------------------------------------

which(!rownames(GOResults$HI) %in% rownames(GOResults$HV))
which(!rownames(GOResults$HO) %in% rownames(GOResults$HV))
which(!rownames(GOResults$IO) %in% rownames(GOResults$HV))
rownames(GOResults$HI)[279]
rownames(GOResults$HO)[279]

GOResults$HI = GOResults$HI[-279,]
GOResults$HO = GOResults$HO[-279,]
GOResults$IO = GOResults$IO[-279,]


nonLiamReuslts = readRDS("Output/CategoricalInsVertivoreTree/")




# ------------------------------------------------------------------
# --- Looking into tree sizes  ----- 
# ------------------------------------------------------------------

mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")
mainTrees$masterTree

unprunedTree = readRDS("Output/Categorical4CategoryUnprunedTree/Categorical4CategoryUnprunedTreeCategoricalTree.rds")
prunedTree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")

prunedTree

?fastwilcoxGMTall 

function (vals, annotList, alternative = "two.sided", ...) 
{
  reslist = list()
  for (n in names(annotList)) {
    reslist[[n]] = fastwilcoxGMT(vals, annotList[[n]], alternative = alternative, 
                                 ...)
    message(paste0(nrow(reslist[[n]]), " results for annotation set ", 
                   n))
  }
  reslist
}
<bytecode: 0x000001ee6f9a39b0>
  <environment: namespace:RERconverge>
# ------------------------------------------------------------------
# --- using liam's zero-length-added-tip ancestral infrence method  ----- 
# ------------------------------------------------------------------

# --------------
# -- section which is actually run during phenotype generation --- 

masterTree = mainTrees$masterTree

nodesToAdd = c(455, 457, 471, 650, 492)
names(nodesToAdd) = c("Mammalia", "Marsupalia", "Placentalia", "Chiroptera", "Primates")

masterTreeAdded = masterTree
for(i in 1:length(nodesToAdd)){
  M<-matchNodes(masterTree,masterTreeAdded,method="distances")
  masterTreeAdded<-bind.tip(masterTreeAdded,names(nodesToAdd)[i],edge.length=0,
                              where=M[which(M[,1]==as.numeric(nodesToAdd[i])),2])
}
mainTrees$masterTree = masterTreeAdded

phenToAdd = c("Insectivore", "Insectivore", "Insectivore", "Insectivore", "Omnivore")
names(phenToAdd) = names(nodesToAdd)

phenotypeVector = append(phenotypeVector, phenToAdd)
speciesFilter = append(speciesFilter, names(phenToAdd))


# - Make common name versions of objects (used in visualization) - 
commonMainTrees = mainTrees
commonMainTrees$masterTree = ZoonomTreeNameToCommon(commonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
commonPhenotypeVector = phenotypeVector
names(commonPhenotypeVector) = ZonomNameConvertVectorCommon(names(commonPhenotypeVector), annotationLocation = spreadSheetLocation, tipColumn = nameColumn)
commonSpeciesFilter = ZonomNameConvertVectorCommon(speciesFilter, annotationLocation = spreadSheetLocation, tipColumn = nameColumn)

# - Categorical Tree - 
treeImageFilename = paste(outputFolderName, filePrefix, "CategoricalTree.pdf", sep="") #make a filename based on the prefix
palette(c( "darkgreen", "darkblue","black", "red"))

pdf(treeImageFilename, height = length(phenotypeVector)/18, width = 10)                     #make a pdf to store the plot, sized based on tree size
commonCategoricalTree = char2TreeCategorical(commonPhenotypeVector, commonMainTrees, commonSpeciesFilter, model = modelType, anctrait = ancestralTrait, plot = F)
categoricalTree = char2TreeCategorical(phenotypeVector, mainTrees, speciesFilter, model = modelType, anctrait = ancestralTrait, plot = F) #use the phenotype vector to make a tree

commonCategoricalTreeExtraTip = commonCategoricalTree
categoricalTreeExtraTip = categoricalTree

commonCategoricalTree = drop.tip(commonCategoricalTree, names(nodesToAdd))
categoricalTree = drop.tip(categoricalTree, names(nodesToAdd))
mainTrees$masterTree = drop.tip(mainTrees$masterTree, names(nodesToAdd))
commonMasterAdded = commonMainTrees$masterTree
commonMainTrees$masterTree = drop.tip(commonMainTrees$masterTree, names(nodesToAdd))

plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMainTrees$masterTree)
plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)

plotTreeCategorical(commonCategoricalTreeExtraTip, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = commonMasterAdded)
plotTreeCategorical(categoricalTreeExtraTip, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = masterTreeAdded)
dev.off()  
# ----------------------------
# comaprision of this with the primary analysis

primaryCategoricalTree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")


length(primaryCategoricalTree$edge)
length(categoricalTree$edge)
all.equal(primaryCategoricalTree$edge,categoricalTree$edge)

length(which(!primaryCategoricalTree$edge == categoricalTree$edge))

all.equal(primaryCategoricalTree$Nnode,categoricalTree$Nnode)
all.equal(primaryCategoricalTree$node.label,categoricalTree$node.label)

all.equal(primaryCategoricalTree$tip.label,categoricalTree$tip.label)

primaryCategoricalTree$tip.label %in% categoricalTree$tip.label
# What this is saying is that the tips are the same, but are stored in a different order. 
# Suggested that save-loading it will reload ordering potentially casued by tip changes) 
write.tree(categoricalTree, "Output/CategoricalInsVertivoreTreeLiamInference/CategoricalTreeRaw.tree")
categoricalTreeNew = read.tree("Output/CategoricalInsVertivoreTreeLiamInference/CategoricalTreeRaw.tree")
all.equal(categoricalTreeNew$tip.label,categoricalTree$tip.label)
#no, that didn't do it. 

categoricalTreeNew = reorder.phylo(categoricalTree)
primaryTreeNew = reorder.phylo(primaryCategoricalTree)
#also no effect

match(categoricalTree$tip.label, primaryCategoricalTree$tip.label)

all.equal(primaryCategoricalTree$edge.length,categoricalTree$edge.length)

library(gridExtra)

liamTree = plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)
primaryTree = plotTreeCategorical(primaryCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)

sideBySide = grid.arrange(liamTree, primaryTree)

pdf(height = length(phenotypeVector)/18, width = 20) 
par(mfrow=c(1,2))
liamTree = plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)
primaryTree = plotTreeCategorical(primaryCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = mainTrees$masterTree)

dev.off()
# ----------
#Code writteen to learn how to do the above 

#Make a plot of the master tree with node labels to find fing nodes
pdf("Output/CategoricalInsVertivoreTree/MasterTreeWithNodeLables.pdf", width = 12, height = 24,)
plotTree(commonMaster,offset=1,direction="rightwards",
         lwd=1)

nodelabels(commonMaster$node.label, cex=0.6, col= "red", frame="none")
nodelabels(cex=0.6, col= "green", frame="none", adj = c(0, 2))

dev.off()

#Nodes to set 
nodesToAdd = c(455, 457, 471, 650, 492)
names(nodesToAdd) = c("Mammalia", "Marsupalia", "Placentalia", "Chiroptera", "Primates")

commonMasterAdded = commonMaster
for(i in 1:length(nodesToAdd)){
  M<-matchNodes(commonMaster,commonMasterAdded,method="distances")
  commonMasterAdded<-bind.tip(commonMasterAdded,names(nodesToAdd)[i],edge.length=0,
                   where=M[which(M[,1]==as.numeric(nodesToAdd[i])),2])
}

pdf("Output/CategoricalInsVertivoreTreeLiamInference/MasterTreeWithAddedNodes.pdf", width = 24, height = 24,)
par(mfrow=c(1,2))
plotTree(commonMaster,fsize=0.8,lwd=1)
nodelabels(cex=0.6)
plotTree(commonMasterAdded,fsize=0.8,lwd=1)
dev.off()

# Liam's code

## load phytools
library(phytools)
## simulate a small tree
tree<-pbtree(n=26,tip.label=LETTERS,scale=1)
## set a value of Q
Q<-matrix(c(-1,1,0,1,-2,1,0,1,-1),3,3,
          dimnames=list(letters[1:3],letters[1:3]))
## simulate data
xx<-sim.Mk(tree,Q,internal=TRUE)
xx

plotTree(tree,offset=1,direction="upwards",
         lwd=1)
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cols<-setNames(palette.colors(3,"Polychrome 36"),
               letters[1:3])
points(pp$xx,pp$yy,pch=16,col=cols[xx],cex=1.5)

nn<-xx[sort(sample(1:tree$Nnode+Ntip(tree),10))]
nn

x<-xx[tree$tip.label]
x

nntree<-tree
for(i in 1:length(nn)){
  M<-matchNodes(tree,nntree,method="distances")
  nntree<-bind.tip(nntree,names(nn)[i],edge.length=0,
                   where=M[which(M[,1]==as.numeric(names(nn)[i])),2])
}

par(mfrow=c(1,2))
plotTree(tree,fsize=0.8,lwd=1)
nodelabels(cex=0.6)
plotTree(nntree,fsize=0.8,lwd=1)

# ------------------------------------------------------------------
# --- looking into data for methods section  ----- 
# ------------------------------------------------------------------
maintrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
masterTree = maintrees$masterTree

otherCategoicalTree = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalTree.rds")

categoricalTree$tip.label[which(!categoricalTree$tip.label %in% otherCategoicalTree$tip.label)]

all.equal(categoricalTree, otherCategoicalTree)

length(droppedTips)

# ------------------------------------------------------------------
# --- Getting named paths for emily ----- 
# ------------------------------------------------------------------

paths = readRDS(pathsFilename)
paths = paths-1
all.equal(binaryPaths, paths)

paths != binaryPaths

# ------------------------------------------------------------------
# --- Looking into I-V signifiacnt results ----- 
# ------------------------------------------------------------------

geneSet = "KeggReactome"

combinedDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResults.rds")
combinedResults = readRDS(combinedDataFilename)

combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResults-", geneSet, ".rds")
GoCombinedResults = readRDS(combinedGODataFilename)

length(which(combinedResults$`IV-significant`))

IVSigGenes = rownames(combinedResults)[which(combinedResults$`IV-significant`)]

length(which(GoCombinedResults$`IV-significant`))

GoCombinedResults[which(GoCombinedResults$`IV-significant`),c(29,30)]


inverseDNARepairGenes = which(GoCombinedResults$`HI-HV-Overlap` & !GoCombinedResults$`HI-HV-CH-Overlap`)
inverseOlfactionGene = which(rownames(GoCombinedResults) == "REACTOME_OLFACTORY_SIGNALING_PATHWAY")
inverseGenes = append(inverseDNARepairGenes, inverseOlfactionGene)
rownames(GoCombinedResults)[inverseGenes]

GoCombinedResults[inverseGenes,c(29,30)]

length(which(GoCombinedResults$`OV-significant`))
GoCombinedResults[which(GoCombinedResults$`OV-significant`),c(35,36)]

vertCandidates = which(rownames(GoCombinedResults) %in% c("REACTOME_GLUCOCORTICOID_BIOSYNTHESIS", "REACTOME_METABOLISM_OF_STEROID_HORMONES"))

sigcols = grep("significant", colnames(GoCombinedResults))

GoCombinedResults[vertCandidates,sigcols]

# ------------------------------------------------------------------
# --- Making RERPlots fro the oxidation genes using the new maintrees ----- 
# ------------------------------------------------------------------
library(data.table)
source("Src/Reu/ZoonomTreeNameToCommon.R")
mainTrees = readRDS("data/categoricalInsVertivoreMaintrees.rds")
mainTrees$masterTree

filePrefix = "CategoricalSlimMainInsVertivoreTree"
outputFolderName = "Output/CategoricalSlimMainInsVertivoreTree/"
cat4phenotypeSet = c("Herbivore", "Insectivore",  "Omnivore", "Vertivore")
cat4colorset = c( "darkgreen", "darkblue","black", "red")

RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
cat4RERObject = readRDS(RERFileName)                                              #Use the existing ones
PathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "")  
pathObject = readRDS(PathsFilename)

commonRERs = cat4RERObject
colnames(commonRERs) = ZonomNameConvertVectorCommon(colnames(commonRERs), tipColumn = "ZoonomiaTip")


palette(c( "darkgreen", "darkblue","black", "red"))


genesOfNote = c("EHHADH", "ECI2", "ACADM", "ACOT12", "ACAD11", "ACOT13", "ACAA2")
k=1
plotRers(commonRERs, genesOfNote[k], pathObject)

pdf("results/temp.pdf", 20, 20)
treePlotRers(mainTrees, cat4RERObject, genesOfNote[k], phenv = pathObject)
dev.off()


mainTrees$trees$EHHADH$tip.label



categoricalTree = readRDS(paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep= ""))

trimmedMainTrees = mainTrees
trimmedMasterTree = drop.tip(mainTrees$masterTree, which(!mainTrees$masterTree$tip.label %in% categoricalTree$tip.label))

trimmedMainTrees$masterTree = trimmedMasterTree
i = 1
for(i in 1:length(trimmedMainTrees$trees)){
  currentTree = trimmedMainTrees$trees[[i]]
  trimmedCurrentTree = drop.tip(currentTree, which(!currentTree$tip.label %in% categoricalTree$tip.label))
  trimmedMainTrees$trees[[i]] = trimmedCurrentTree
}


RERTree = returnRersAsTree(trimmedMainTrees, cat4RERObject, genesOfNote[k])


# -- Switching to making a new maintrees object fro returnRERs that's pruned to the right size. 

write.tree(trimmedMasterTree, "Results/InsVertMasterTree.tree")


# ------------------------------------------------------------------
# ---  -Making scatterplots of the DNA repair genes---- 
# ------------------------------------------------------------------
#Run alongside the MakeOVerlapFigure main script 
combinedDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResults.rds")
combinedResults = readRDS(combinedDataFilename)

combinedGODataFilename = paste0(outputFolderName, filePrefix, "combinedGOResults-", geneSet, ".rds")
GoCombinedResults = readRDS(combinedGODataFilename)

DNARepairGenesFilename = "Results/dnaRepairGenes.rds"
genesets = readRDS(DNARepairGenesFilename)


rhoValuesSpecific = rhoValues
specificGenes = which(rownames(rhoValues) %in% genesets$combined)

rhoValuesSpecific = rhoValues[specificGenes, ]

library(ggrepel)

for(i in 1:length(rhoValuesSpecific)){
  xName = names(rhoValuesSpecific)[i]
  if(i+1 <= length(rhoValuesSpecific)){
    for(j in (i+1):length(rhoValuesSpecific)){
      yName = names(rhoValuesSpecific)[j]
      
      rhoCorrellPlot = ggplot(rhoValuesSpecific, aes(x = .data[[xName]], y = .data[[yName]])) + 
        geom_point() + geom_pointdensity() + scale_color_viridis()
      
      
      denstiyScaleValue = ggplot_build(rhoCorrellPlot)$plot$scales$scales[[1]]$get_limits()[2]
      densityScaleSet = append(densityScaleSet, denstiyScaleValue)
      rm(rhoCorrellPlot)
    }
  }
}
densityScale = c(1, max(densityScaleSet))


rhoPlotSet = list()
netIndex= 0
for(i in 1:length(rhoValuesSpecific)){
  xName = names(rhoValuesSpecific)[i]
  if(i <= length(rhoValuesSpecific)){
    if(bothAxis){jStart = 1}else{jStart = i+1}
    for(j in (jStart):length(rhoValuesSpecific)){
      yName = names(rhoValuesSpecific)[j]
      yLabel =  paste0(replacePrefixWithName(addDashes(gsub("-Rho", "", yName))), " Dunn Z Statistic")
      xLabel =  paste0(replacePrefixWithName(addDashes(gsub("-Rho", "", xName))), " Dunn Z Statistic")
      
      rhoCorrellPlot = ggplot(rhoValuesSpecific, aes(x = .data[[xName]], y = .data[[yName]])) + 
        geom_point() + geom_pointdensity() + scale_color_viridis(name = "Density of genes", limits = densityScale) + 
        stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),formula = y ~ x,parse = TRUE, size = 6) +
        theme_classic()+
        geom_text_repel(aes(label = rownames(rhoValuesSpecific)), size = 3)+
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") +  # vertical line at x=0
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  # horizontal line at y=0
        xlab(xLabel) + ylab(yLabel)+
        theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16))
      
      netIndex = netIndex +1
      rhoPlotSet[[netIndex]] = rhoCorrellPlot
      names(rhoPlotSet)[netIndex] = paste(xName, yName, sep="-")
      rm(rhoCorrellPlot)
    }
  }
}

rhoPlotSet$`HI-Rho-HV-Rho`

png("Results/DNARepairGenes.png", 1200,1200)
rhoPlotSet$`HI-Rho-HV-Rho`
dev.off()
# ------------------------------------------------------------------
# --- Making RERPlots fro the oxidation genes ----- 
# ------------------------------------------------------------------
library(data.table)
source("Src/Reu/ZoonomTreeNameToCommon.R")
mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")


filePrefix = "CategoricalInsVertivoreTree"
outputFolderName = "Output/CategoricalInsVertivoreTree/"
cat4phenotypeSet = c("Herbivore", "Insectivore",  "Omnivore", "Vertivore")
cat4colorset = c( "darkgreen", "darkblue","black", "red")

RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
cat4RERObject = readRDS(RERFileName)                                              #Use the existing ones
PathsFilename = paste(outputFolderName, filePrefix, "CategoricalPathsFile.rds", sep= "")  
pathObject = readRDS(PathsFilename)

commonRERs = cat4RERObject
colnames(commonRERs) = ZonomNameConvertVectorCommon(colnames(commonRERs), tipColumn = "ZoonomiaTip")


palette(c( "darkgreen", "darkblue","black", "red"))
plotRers(commonRERs, "ACAA2", pathObject)

genesOfNote = c("EHHADH", "ECI2", "ACADM", "ACOT12", "ACAD11", "ACOT13", "ACAA2")

k=1

pdf("results/temp.pdf", 20, 20)
treePlotRers(mainTrees, cat4RERObject, genesOfNote[k], phenv = pathObject)
dev.off()


categoricalTree = readRDS(paste(outputFolderName, filePrefix, "CategoricalTree.rds", sep= ""))

trimmedMainTrees = mainTrees
trimmedMasterTree = drop.tip(mainTrees$masterTree, which(!mainTrees$masterTree$tip.label %in% categoricalTree$tip.label))

trimmedMainTrees$masterTree = trimmedMasterTree
i = 1
for(i in 1:length(trimmedMainTrees$trees)){
  currentTree = trimmedMainTrees$trees[[i]]
  trimmedCurrentTree = drop.tip(currentTree, which(!currentTree$tip.label %in% categoricalTree$tip.label))
  trimmedMainTrees$trees[[i]] = trimmedCurrentTree
}


RERTree = returnRersAsTree(trimmedMainTrees, cat4RERObject, genesOfNote[k])


# -- Switching to making a new maintrees object fro returnRERs that's pruned to the right size. 

write.tree(trimmedMasterTree, "Results/InsVertMasterTree.tree")

# ------------------------------------------------------------------
# ---  Getting gene list of the pathways in opposite directions  ----- 
# ------------------------------------------------------------------
library(RERconverge)
library("tools")
source("Src/Reu/cmdArgImport.R")
library(xlsx)

pathwayNames = c("REACTOME_HDR_THROUGH_SINGLE_STRAND_ANNEALING_SSA", "REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR", "REACTOME_DNA_DOUBLE_STRAND_BREAK_REPAIR", "REACTOME_HOMOLOGY_DIRECTED_REPAIR", "REACTOME_DISEASES_OF_DNA_REPAIR")

gmtFileName = "Data/KeggReactome.gmt"
gmtFile = read.gmt(gmtFileName)

which(gmtFile$geneset.names %in% pathwayNames)

genesets = gmtFile$genesets[which(gmtFile$geneset.names %in% pathwayNames)]
names(genesets) = gmtFile$geneset.names[which(gmtFile$geneset.names %in% pathwayNames)]

unlist(genesets[1:5])

genesets$combined = list(unlist(genesets[1:5]))
genesets$combined = unlist(genesets[1:5])
genesets$combined = unique(genesets$combined)

DNARepairGenesFilename = "Results/dnaRepairGenes.rds"
saveRDS(genesets, DNARepairGenesFilename)

lines <- sapply(genesets, function(x) paste(x, collapse = "\t"))

writeLines(lines, "Results/dnaRepairGenes.txt")

names(genesets)

which(is.numeric(names(droppedTips)))
which(!is.na(as.numeric(names(droppedTips))))
length(which(!is.na(as.numeric(names(droppedTips)))))
# ------------------------------------------------------------------
# ---  Make plots using the branchlength removed mastertrees  ----- 
# ------------------------------------------------------------------
stableMaintrees = mainTrees 
mainTrees = stableMaintrees
stableMaintrees = readRDS(mainTreesLocation)

stableCommonMainTrees = stableMaintrees
stableCommonMainTrees$masterTree = ZoonomTreeNameToCommon(stableCommonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)


mainTrees$masterTree$edge.length[1:length(mainTrees$masterTree$edge.length)] = 1

pdf(treeImageFilename, height = length(phenotypeVector)/14, width = 10)     
plotTreeCategorical(commonCategoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableCommonMainTrees$masterTree)

plotTreeCategorical(categoricalTree, c("Herbivore", "Insectivore", "Omnivore", "Vertivore"), master = stableMaintrees$masterTree)
dev.off()  

categoricalCommonTreeFilename = paste(outputFolderName, filePrefix, "CategoricalCommonTree.rds", sep="") #make a filename based on the prefix
saveRDS(commonCategoricalTree, categoricalCommonTreeFilename)

plotTreeCategorical(categoricalTree, c("Carnivore", "Herbivore", "Omnivore"), master = stableMaintrees$masterTree)

plotTreeCategorical(commonCategoricalTree, c("Carnivore", "Herbivore", "Omnivore"), master = stableCommonMainTrees$masterTree)

pdf(treeImageFilename, height = length(phenotypeVector)/14, width = 10)     
plotTreeCategorical(commonCategoricalTree, c("Background", "Carnivore"), master = stableCommonMainTrees$masterTree)

plotTreeCategorical(categoricalTree, c("Background", "Carnivore"), master = stableMaintrees$masterTree)
dev.off()  



# ------------------------------------------------------------------
# ---  Make subset tree for tyler ----- 
# ------------------------------------------------------------------

tylerSpecies = read.csv("Results/toMichael.csv")


tipsToKeep = tylerSpecies$fa


commonCategoricalTree

tipsToKeep = ZonomNameConvertVectorCommon(tipsToKeep, tipColumn = "ZoonomiaTip")
tipsToKeep = commonCategoricalTree$tip.label[which(1:64 %% 2 ==0)]
tipsToDrop = commonCategoricalTree$tip.label[!commonCategoricalTree$tip.label %in% tipsToKeep]
commonCategoricalTreePruned = drop.tip(commonCategoricalTree, tipsToDrop)

commonCategoricalTree = commonCategoricalTreePruned

ggTreeOut = ggtree(commonCategoricalTree) +scale_color_manual(values=palette()) 
ggTreeOut = ggTreeOut %<+% edge + aes(color=CategorylengthChar)
ggTreeOut = ggTreeOut %<+% tip_data 
ggTreeOut$data$label = paste(ggTreeOut$data$label, "-", ggTreeOut$data$node, sep="")
ggTreeOut = ggTreeOut + geom_tiplab()
#ggTreeOut = ggTreeOut + geom_tiplab(geom = "phylopic", aes(image = uuid))
#ggTreeOut + geom_phylopic(aes(uuid = uuid), color = "black", alpha = 1, size = 0.08)
ggTreeOut




# ------------------------------------------------------------------
# ---  Confirm Rho meaning for Categorical Diet ----- 
# ------------------------------------------------------------------


results = readRDS("OUtput/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCombinedCategoricalCorrelationFile.rds")
RERObject = readRDS("Output/CategoricalBinaryCarnivoreTree/CategoricalBinaryCarnivoreTreeRERFile.rds")
phenotypeVector = readRDS("Output/CategoricalBinaryCarnivoreTree/CategoricalBinaryCarnivoreTreeCategoricalPhenotypeVector.rds")
pathsObject = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPathsFile.rds")
mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
phenotypeSet = c("Herbivore", "Insectivore", "Omnivore", "Vertivore")
colorset = c( "darkgreen", "darkblue", "black", "red")

plotRers(RERObject, "CPB1" , pathsObject,)
source("Src/Reu/rerViolinPlot.R")
rerViolinPlot(mainTrees, RERObject, pathsObject, phenotypeSet , geneOfInterest = "CPB1", colorScale = colorset)

pairwiseResults = results[[2]]

genelist = c("CPB1")
rowindex = which(rownames(pairwiseResults$`1 - 2`)%in% genelist)

pairwiseResults$`1 - 2`[rowindex,]

# ------------------------------------------------------------------
# ---  Making Demos of the plots or other main-repo addable functions ----- 
# ------------------------------------------------------------------

#MakeMasterAndGeneTreePlot 

source("Src/Reu/makeMasterAndGeneTreePlots.R")

mainTrees = readRDS("Data/zoonomiaAllMammalsTrees.rds")
RERObject = readRDS("Output/CategoricalBinaryCarnivoreTree/CategoricalBinaryCarnivoreTreeRERFile.rds")
phenotypeVector = readRDS("Output/CategoricalBinaryCarnivoreTree/CategoricalBinaryCarnivoreTreeCategoricalPhenotypeVector.rds")
geneInQuestion = "CFTR"
dropTips = T
convertNames = T
tipColumn = "ZoonomiaTip"

png("Results/MasterAndGeneExmaple.png", width = 1200, height = 1200)
makeMasterAndGeneTreePlots(mainTrees, "M6PR", RERObject, tipColumn = "ZoonomiaTip", phenotypeVector = phenotypeVector, fgcols = "orange", bgcolor = "darkgreen")
dev.off()


#RER violin plot 
source("Src/Reu/rerViolinPlot.R")

RerData = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeRERFile.rds")
PathsData = readRDS("Output/CategoricalInsVertivoreTree/CategoricalInsVertivoreTreeCategoricalPathsFile.rds")
PhenotypeSet = c("Herbivore", "Insectivore",  "Omnivore", "Vertivore")
Colorset = c( "darkgreen", "darkblue","black", "red")
geneOfInterest = "CYP1A1"

png("Results/RERViolinPlotExample.png", 600, 600)
rerViolinPlot(mainTrees, RerData, pathsObject = PathsData, phenotypeSet = PhenotypeSet, colorScale = Colorset, geneOfInterest = geneOfInterest)
dev.off()

# ------------------------------------------------------------------
# --- Checking if mergedata has all hiller speices   ----- 
# ------------------------------------------------------------------
mergeData = read.csv("Data/mergedData.csv")
mainTrees = readRDS('data/zoonomiaAllMammalsTrees.rds')

mainTrees$masterTree$tip.label[which(!mainTrees$masterTree$tip.label %in% mergeData$ZoonomiaTip)]

which(mergeData)


# ------------------------------------------------------------------
# --- Making a column in mergeData that matchs the InsVertivore classification   ----- 
# ------------------------------------------------------------------

substitutions = list(
  c("C-Invertebrate-eater", "Insectivore"), c("C-InsVertivore-Insectivore", "Insectivore"),
  c("C-Herpetivore", "Vertivore"),
  c("C-Piscivore", "Vertivore"), c("C-InsVertivore-Piscivore", "Vertivore"),
  c("C-Endotherm-Carnivore", "Vertivore"), c("C-Scavenger", "Vertivore"), c("C-Nonspecific-Vertebrate-eater", "Vertivore"),
  c("C-Terrestrial-vertebrates-eater", "Vertivore"), c("C-All-vertebrate-eater", "Vertivore"), c("C-InsVertivore-Carnivore", "Vertivore"),
  c("C-InsVertivore-Mixed", "Omnivore"), 
  c("O-For Examination", "Omnivore"), c("O-Scavenger", "Omnivore"),
  c("H-Frugivore", "Herbivore"), 
  c("H-Nectarivore", "Herbivore"), 
  c("H-High-sugar-plants-Eater", "Herbivore"),
  c("H-Granivore", "Herbivore"), c("H-Nonspecific-Herbivore", "Herbivore"), 
  c("H-Low-sugar-plants-Eater", "Herbivore"), c("H-All-plants-Eater", "Herbivore"),
  c("O-Generalist", "Omnivore")
)
mergeData = read.csv("Data/mergedData.csv")

Dietvalues = mergeData$DerekDietClassification90InsVertivoreSorting

  for( i in 1:length(substitutions)){
    substitutePhenotypes = substitutions[[i]]
    message(paste("replacing", substitutePhenotypes[1], "with", substitutePhenotypes[2]))
    Dietvalues = gsub(substitutePhenotypes[1], substitutePhenotypes[2], Dietvalues)
  }

mergeData$insVertivoreDiet = Dietvalues

write.csv(mergeData, "Data/mergedData.csv", row.names = F)
# ------------------------------------------------------------------
# --- Making Updated RERConverge Explanation slide  ----- 
# ------------------------------------------------------------------
library(RERconverge)
source("Src/Reu/ZonomNameConvertMatrixCommon.R")
mainTrees = readRDS("Data/RemadeTreesAllZoonomiaSpecies.rds")
mainTrees = readRDS("data/zoonomiaAllMammalsTrees.rds")

CVHRERs = readRDS("Output/Old/CVHRemake/CVHRemakeRERFile.rds")
foregroundSpecies = readRDS("Output/Old/CVHRemake/CVHRemakeBinaryTreeForegroundSpecies.rds")
CVHPaths = readRDS("Output/Old/CVHRemake/CVHRemakePathsFile.rds")
commonRERs = ZonomNameConvertMatrixCommon(CVHRERs)


source("Src/Reu/makeMasterAndGeneTreePlots.R")

# this one is correctly sized for the most part, but the tip labels are too small (especially given the lighter orange). I don't know how to fix that -- changing the obvious values in the plotting funciton had no effect. 
png("Results/tempMasterandGenePlot.png", 575, 575)
makeMasterAndGeneTreePlots(mainTrees,"IQANK1", CVHRERs,  foregroundSpecies, correlationPlot = F, tipColumn = "manualAnnotations_FaName", fgcols = "orange", bgcolor = "darkgreen")
dev.off()

png("Results/tempCorrelationPlot.png", 420, 420)
makeMasterAndGeneTreePlots(mainTrees,"IQANK1", CVHRERs,  foregroundSpecies, correlationPlot = T, tipColumn = "manualAnnotations_FaName", fgcols = "orange", bgcolor = "darkgreen")
dev.off()


plotRers(commonRERs, "IQANK1", CVHPaths, sort = F)
plotRersNew = function (rermat = NULL, index = NULL, phenv = NULL, rers = NULL, method = "k", xlims = NULL, plot = 1, xextend = 0.2, sortrers = F) {
  {
    if (!is.null(phenv) && length(unique(phenv[!is.na(phenv)])) > 
        2) {
      categorical = TRUE
      if (method != "aov") {
        method = "kw"
      }
    }
    else {
      categorical = FALSE
    }
    if (is.null(rers)) {
      e1 = rermat[index, ][!is.na(rermat[index, ])]
      colids = !is.na(rermat[index, ])
      e1plot <- e1
      if (exists("speciesNames")) {
        names(e1plot) <- speciesNames[names(e1), ]
      }
      if (is.numeric(index)) {
        gen = rownames(rermat)[index]
      }
      else {
        gen = index
      }
    }
    else {
      e1plot = rers
      gen = "rates"
    }
    names(e1plot)[is.na(names(e1plot))] = ""
    if (!is.null(phenv)) {
      phenvid = phenv[colids]
      if (categorical) {
        fgdcor = getAllCor(rermat[index, , drop = F], phenv, 
                           method = method)[[1]]
      }
      else {
        fgdcor = getAllCor(rermat[index, , drop = F], phenv, 
                           method = method)
      }
      plottitle = paste0(gen, ": rho = ", round(fgdcor$Rho, 
                                                4), ", p = ", round(fgdcor$P, 4))
      if (categorical) {
        n = length(unique(phenvid))
        if (n > length(palette())) {
          pal = colorRampPalette(palette())(n)
        }
        else {
          pal = palette()[1:n]
        }
      }
      if (categorical) {
        df <- data.frame(species = names(e1plot), rer = e1plot, 
                         stringsAsFactors = FALSE) %>% mutate(mole = as.factor(phenvid))
      }
      else {
        df <- data.frame(species = names(e1plot), rer = e1plot, 
                         stringsAsFactors = FALSE) %>% mutate(mole = as.factor(ifelse(phenvid > 
                                                                                        0, 2, 1)))
      }
    }
    else {
      plottitle = gen
      df <- data.frame(species = names(e1plot), rer = e1plot, 
                       stringsAsFactors = FALSE) %>% mutate(mole = as.factor(ifelse(0, 
                                                                                    2, 1)))
    }
    if (sortrers) {
      df = filter(df, species != "") %>% arrange(desc(rer))
    }
    if (is.null(xlims)) {
      ll = c(min(df$rer) * 1.1, max(df$rer) + xextend)
    }
    else {
      ll = xlims
    }
  }
  if (categorical) {
    g <- ggplot(df, aes(x = rer, y = factor(species, levels = unique(ifelse(rep(sortrers, 
                                                                                nrow(df)), species[order(rer)], sort(unique(species))))), 
                        col = mole, label = species)) + scale_size_manual(values = c(1, 
                                                                                     1, 1, 1)) + geom_point(aes(size = mole)) + scale_color_manual(values = pal) + 
      scale_x_continuous(limits = ll) + geom_text(hjust = 1, 
                                                  size = 2) + ylab("Branches") + xlab("relative rate") + 
      ggtitle(plottitle) + geom_vline(xintercept = 0, linetype = "dotted") + 
      theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), 
            legend.position = "none", panel.background = element_blank(), 
            axis.text = element_text(size = 18, face = "bold", 
                                     colour = "black"), axis.title = element_text(size = 24, 
                                                                                  face = "bold"), plot.title = element_text(size = 24, 
                                                                                                                            face = "bold")) + theme(axis.line = element_line(colour = "black", 
                                                                                                                                                                             size = 1)) + theme(axis.line.y = element_blank())
  }
  else {
    g <- ggplot(df, 
                aes(x = rer, 
                    y = factor(species, levels = unique(ifelse(rep(sortrers, nrow(df)), species[order(rer)], sort(unique(species))))), 
                    col = mole, 
                    label = species
                )
    ) + 
      scale_size_manual(values = c(1, 1, 1, 1)) + 
      geom_point(aes(size = mole)) + 
      scale_color_manual(values = c("black", "blue")) + 
      scale_x_continuous(limits = ll) + 
      geom_text(hjust = "center", size = 3) + 
      ylab("Branches") + 
      xlab("relative rate") + 
      ggtitle(plottitle) + 
      geom_vline(xintercept = 0, linetype = "dotted") + 
      theme(
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.position = "none", 
        panel.background = element_blank(), 
        axis.text = element_text(size = 18, face = "bold", colour = "black"), 
        axis.title = element_text(size = 24, face = "bold"), 
        plot.title = element_text(size = 24, face = "bold")) + 
      theme(axis.line = element_line(colour = "black", size = 1)) + 
      theme(axis.line.y = element_blank())
  }
  if (plot) {
    print(g)
  }
  else {
    g
  }
}



mainTrees



# ------------------------------------------------------------------
# --- Getting trees for jack  ----- 
# ------------------------------------------------------------------
categoricalTreeFilename

categoricalTestTree = readRDS(categoricalTreeFilename)

masterTreeWithBranchLengths = stableMaintrees$masterTree
treeInMasterWithoutPhenotype = masterTreeWithBranchLengths$tip.label[!masterTreeWithBranchLengths$tip.label %in% categoricalTestTree$tip.label]
masterTreeWithBranchLengthsPruned = drop.tip(masterTreeWithBranchLengths, treeInMasterWithoutPhenotype)


togaTree = read.newick("Data/TogaTree.nwk")
treeInHillerWithoutPhenotype = togaTree$tip.label[!togaTree$tip.label %in% categoricalTestTree$tip.label]
hillerTreePruned = drop.tip(togaTree, treeInHillerWithoutPhenotype)

hillerTreePruned$edge.length
masterTreeWithBranchLengthsPruned$edge.length


write.tree(categoricalTestTree, paste0(outputFolderName, "ZoonomiaMaximalTreeCategoryBranchLengths.nwk"))
write.tree(masterTreeWithBranchLengthsPruned, paste0(outputFolderName, "ZoonomiaMaximalTreeAlignmentBranchLengths.nwk"))
write.tree(hillerTreePruned, paste0(outputFolderName, "ZoonomiaMaximalTreeHillerBranchLengths.nwk"))
saveRDS(categoricalTestTree, paste0(outputFolderName, "ZoonomiaMaximalTreeCategoryBranchLengths.rds"))
saveRDS(masterTreeWithBranchLengthsPruned, paste0(outputFolderName, "ZoonomiaMaximalTreeAlignmentBranchLengths.rds"))
saveRDS(hillerTreePruned, paste0(outputFolderName, "ZoonomiaMaximalTreeHillerBranchLengths.rds"))


# ------------------------------------------------------------------
# --- Making violin plots from the categorical data for presentaiton  ----- 
# ------------------------------------------------------------------

filePrefix = "CategoricalInsVertivoreTree"
outputFolderName = "Output/CategoricalInsVertivoreTree/"
cat4phenotypeSet = c("Herbivore", "Insectivore",  "Omnivore", "Vertivore")
cat4colorset = c( "darkgreen", "darkblue","black", "red")

RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
cat4RERObject = readRDS(RERFileName)                                              #Use the existing ones

pathsFileName = paste(outputFolderName, filePrefix, phenotypeStyle, "PathsFile.rds", sep= "") #Set a filename for the pathss based on the prefix and style
cat4pathsObject = readRDS(pathsFileName)                                          #If the file already exists, use the existing one.

combinedDataFilename = paste0(outputFolderName, filePrefix, "combinedGeneResults.rds")
combinedResults = readRDS(combinedDataFilename)

overlap = combinedResults[combinedResults$`HI-HV-CH-Overlap` & !is.na(combinedResults$`HI-HV-CH-Overlap`) & !combinedResults$`HO-significant`,]

overlapGenes = rownames(overlap[order(overlap$`HI-p.adj`),])

filePrefix = "CategoricalPrunedCarnivoreTree"
outputFolderName = "Output/CategoricalPrunedCarnivoreTree/"
cat3phenotypeSet = c("Carnivore", "Herbivore", "Omnivore")
cat3colorset = c( "orange", "darkgreen", "black")

RERFileName = paste(outputFolderName, filePrefix, "RERFile.rds", sep= "")       #Set a filename for the RERs based on the prefix
cat3RERObject = readRDS(RERFileName)                                              #Use the existing ones

pathsFileName = paste(outputFolderName, filePrefix, phenotypeStyle, "PathsFile.rds", sep= "") #Set a filename for the pathss based on the prefix and style
cat3pathsObject = readRDS(pathsFileName)                                          #If the file already exists, use the existing one.



if(file_ext(mainTreesLocation) == "rds"){
  if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
}else{
  if(!exists("mainTrees")){mainTrees = readTrees(mainTreesLocation)} 
}





source("Src/Reu/rerViolinPlot.R")
library(gridExtra)

i = 1
{
currentGene = overlapGenes[i]
cat4Plot = rerViolinPlot(mainTrees, cat4RERObject, cat4pathsObject, cat4phenotypeSet , geneOfInterest = currentGene, colorScale = cat4colorset)
cat4Plot = cat4Plot + xlab(c("Herbivore", "Invertivore", "Omnivore", "Vertivore"))
cat3Plot = rerViolinPlot(mainTrees, cat3RERObject, cat3pathsObject, cat3phenotypeSet , geneOfInterest = currentGene, colorScale = cat3colorset)
comboPlot = grid.arrange(cat4Plot, cat3Plot, ncol = 2)
comboPlot
i = i+1
}
goodGenes = c(1, 8, 10)
# ------------------------------------------------------------------
# --- Work on making a tree figure for lalitha  ----- 
# ------------------------------------------------------------------
lalithaData = read.csv("C:/Users/mit221/Downloads/SpeciesNamesAndPhenosForDE72.csv")


laltihaSpecies = lalithaData[,1]

#laltihaSpecies = unlist(laltihaSpecies)
#laltihaSpecies = laltihaSpecies[-1]
#laltihaSpecies = gsub("[0-9]+$", "", laltihaSpecies)
#laltihaSpecies = unique(laltihaSpecies)

mergeData = read.csv("Data/mergedData.csv")

for(i in 1:length(laltihaSpecies)){
laltihaSpecies[i] = gsub(" ", "_", laltihaSpecies[i]) #replace spaces with underscores 
laltihaSpecies[i] = sub('^([^_]+_[^_]+).*', '\\1', laltihaSpecies[i]) #remove anything  after a second underscore
laltihaSpecies[i] = tolower(laltihaSpecies[i])
}

lalithaData[,1] = laltihaSpecies
write.csv(lalithaData, "C:/Users/mit221/Downloads/SpeciesNamesAndPhenosForDE72.csv")

sum(!laltihaSpecies %in% mergeData$Scientific_Binomial)


ggTreeOut = ggtree(commonCategoricalTree) +scale_color_manual(values=palette()) 




# ------------------------------------------------------------------
# --- Assess driving diet of unqieu carnivory results  ----- 
# ------------------------------------------------------------------

combinedData = readRDS(paste0(combinedDataFilename, ".rds"))

carnivoryUNique = combinedData[which(combinedData$`CH-significant` & !combinedData$`HI-significant` & !combinedData$`HV-significant`),]


table(carnivoryUNique$`CH-Driver`)
table(sign(carnivoryUNique$`CH-Rho`))
# Right. I don't current have driver information because I haven't run the driver analysis for the 
# I need to see about making that. 

# quick fix for the driver analysis to retarget because this is not technically a part of the main analysis 

# ------------------------------------------------------------------
# --- Get set of genes most different between HI and HV  ----- 
# ------------------------------------------------------------------

combinedResults # from MakeOverlapFigure
combinedResultsModified = combinedResults


getComparisionDifference = function(dataframe, colOne, colTwo){
  colOneIndex = names(dataframe)[which(names(dataframe) == colOne)]
  colTwoIndex = names(dataframe)[which(names(dataframe) == colTwo)]
  
  distanceFromEqual = abs(dataframe[colOneIndex] - dataframe[colTwoIndex]) / sqrt(2)
  distanceFromEqual
}



getComparisionDifference(combinedResultsModified, "HI-Rho", "HV-Rho")
combinedResults$`HI-HV-Delta` = getComparisionDifference(combinedResultsModified, "HI-Rho", "HV-Rho")


test = combinedResultsModified[order(combinedResultsModified$`HI-HV-Delta`, decreasing = T),]


df$distance_from_y_eq_x <- abs(df$x - df$y) / sqrt(2)




# ------------------------------------------------------------------
# --- examining the GO sets in specific overlap sections ----- 
# ------------------------------------------------------------------
GoSignificanceResults # from MakeOverlapFigure

View(GoSignificanceResults)

# get pathways that show up in compionents but not carnivory 
IVnCpathways = rownames(GoSignificanceResults[which(GoSignificanceResults$`HI-significant` & GoSignificanceResults$`HV-significant` & !GoSignificanceResults$`CH-significant`),])

which(rownames(GoCombinedResults) %in% IVnCpathways)
IVnCResults = GoCombinedResults[which(rownames(GoCombinedResults) %in% IVnCpathways), ]
View(IVnCResults)
cat(IVnCpathways)

# get the pathways which are unique to vertivory 
VonlyPathways = rownames(GoSignificanceResults[which(!GoSignificanceResults$`HI-significant` & GoSignificanceResults$`HV-significant` & !GoSignificanceResults$`CH-significant`),])
IonlyPathways = rownames(GoSignificanceResults[which(GoSignificanceResults$`HI-significant` & !GoSignificanceResults$`HV-significant` & !GoSignificanceResults$`CH-significant`),])


# ------------------------------------------------------------------
# --- getting permualtion RERResult data fro rho plot creation ----- 
# ------------------------------------------------------------------

permulationIntermediate = readRDS(paste0(outputFolderName, "CategoricalInsVertivoreTreeCategoricalPermulationsIntermediates101.rds"))


permulationHIStatsValues = permulationIntermediate$Peffsize$`1 - 3`
permulationHVStatsValues = permulationIntermediate$Peffsize$`1 - 4`

permulationHIStatsCol = permulationHIStatsValues[,1]
permulationHVStatsCol = permulationHVStatsValues[,1]


rhoValuesPerm = data.frame(permulationHIStatsCol, permulationHVStatsCol)

for(i in 1:length(rhoValuesPerm)){
  xName = names(rhoValuesPerm)[i]
  if(i+1 <= length(rhoValuesPerm)){
    for(j in (i+1):length(rhoValuesPerm)){
      yName = names(rhoValuesPerm)[j]
      yLabel =  paste0(replacePrefixWithName(addDashes(gsub("-Rho", "", yName))), " Stat")
      xLabel =  paste0(replacePrefixWithName(addDashes(gsub("-Rho", "", xName))), " Stat")
      
      rhoCorrellPlot = ggplot(rhoValuesPerm, aes(x = .data[[xName]], y = .data[[yName]])) + 
        geom_point() + geom_pointdensity() + scale_color_viridis(name = "Number of nearby genes", limits = densityScale) + 
        stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),formula = y ~ x,parse = TRUE) +
        theme_classic()+
        xlab(xLabel) + ylab(yLabel)
    
      
      netIndex = netIndex +1
      rhoPlotSet[[netIndex]] = rhoCorrellPlot
      names(rhoPlotSet)[netIndex] = paste(xName, yName, sep="-")
    }
  }
}

rhoDfList = list()
permulationRhoDataframes = for(i in 1:2500){
  rhoDf = data.frame(permulationHIStatsValues[,i], permulationHVStatsValues[,i])
  rhoDfList[[i]] = rhoDf
}

rSquaredSet = NULL
for(k in 1:length(rhoDfList)){
  rhoDf = rhoDfList[[k]]
  names(rhoDf) = c("HI", "HV")
  model <- lm(HV ~ HI, data = rhoDf)
  rSquared = summary(model)$r.squared
  rSquaredSet = append(rSquaredSet, rSquared)
  #cat(rSquared, "\n")
}
pdf()
hist(rSquaredSet)
dev.off()

mean(rSquaredSet)
