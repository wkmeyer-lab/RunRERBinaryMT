a = b #This is to prevent accidental ful runs 

permdata = readRDS("Output/CategoricalDiet3Phen/CategoricalDiet3PhenPermulationsData1.rds")

permdtaaTrim = permdata

permdtaaTrim[[2]] = permdtaaTrim[[2]][1:3]
permdtaaTrim[[3]] = permdtaaTrim[[3]][1:3]

permDataTrim2 = permdata
permDataTrim2[[2]] = permDataTrim2[[2]][4:6]
permDataTrim2[[3]] = permDataTrim2[[3]][4:6]

permDataTrim3 = permdata
permDataTrim3[[2]] = permDataTrim3[[2]][7:9]
permDataTrim3[[3]] = permDataTrim3[[3]][7:9]

permDataTrim4 = permdata
permDataTrim4[[2]] = permDataTrim4[[2]][10:12]
permDataTrim4[[3]] = permDataTrim4[[3]][10:12]

saveRDS(permdtaaTrim, "Output/CategoricalDiet3Phen/CategoricalDiet3PhenPermulationsDataDev.rds")
rm(permdata)
rm(permdtaaTrim)

library(data.table)
source("Src/Reu/CategoricalPermulationsParallelFunctions.R")

permCorrelations = getPermPvalsCategorical3(correlationsObject, permulationsData$trees, phenotypeVector, mainTrees, RERObject)
permCorrelations = getPermPvalsCategorical2(correlationsObject, permulationsData$trees, phenotypeVector, mainTrees, RERObject)

permulationsDataSaving = permulationsData

permulationsData = permulationsDataSaving
permCorrelationsNew = CategoricalPermulationGetCor(correlationsObject, permulationsData$trees, phenotypeVector, mainTrees, RERObject, report=T)


permCorrelationsNew2 = CategoricalPermulationGetCor(correlationsObject, permDataTrim2$trees, phenotypeVector, mainTrees, RERObject, report=T)
permCorrelationsCombine = permCorrelationsNew


intermediate1 = permCorrelationsNew
intermediate2 = permCorrelationsNew2
i=1

combineCategoricalPermulationIntermediates = function(intermediate1, intermediate2){
  intermediatesCombined = intermediate1
  intermediatesCombined[[1]] = cbind(intermediatesCombined[[1]], intermediate2[[1]])
  for(i in 1:length(intermediatesCombined[[2]])){
    intermediatesCombined[[2]][[i]] = cbind(intermediatesCombined[[2]][[i]], intermediate2[[2]][[i]])
  }
  intermediatesCombined[[3]] = cbind(intermediatesCombined[[3]], intermediate2[[3]])
  for(i in 1:length(intermediatesCombined[[4]])){
    intermediatesCombined[[4]][[i]] = cbind(intermediatesCombined[[4]][[i]], intermediate2[[4]][[i]])
  }
  return(intermediatesCombined)
}

paste(basePermulationsFilename, 1, ".rds", sep= "")

permCorrelationsNew3 = CategoricalPermulationGetCor(correlationsObject, permDataTrim3$trees, phenotypeVector, mainTrees, RERObject, report=T)
permCorrelationsNew4 = CategoricalPermulationGetCor(correlationsObject, permDataTrim4$trees, phenotypeVector, mainTrees, RERObject, report=T)


saveRDS(permCorrelationsNew, paste(basePermulationsFilename, 1, ".rds", sep= ""))
saveRDS(permCorrelationsNew2, paste(basePermulationsFilename, 2, ".rds", sep= ""))
saveRDS(permCorrelationsNew3, paste(basePermulationsFilename, 3, ".rds", sep= ""))
saveRDS(permCorrelationsNew4, paste(basePermulationsFilename, 4, ".rds", sep= ""))


combinedPerms = combineCategoricalPermulationIntermediates(intermediate1, intermediate2)



ist1 = permCorrelationsNew
list2 = permCorrelationsNew2
listIndex = 1
matrix1 = list1[[listIndex]]
matrix2 = list2[[listIndex]]

matrixCombine = cbind(matrix1, matrix2)
matrix1Length = ncol(matrix1)
matrix2Length = ncol(matrix2)

matrixCombine[,4] = matrix2[,1]


permCorrelationsCombine[[1]][,matrix1Length+1:matrix1Length+matrix2Length] = permCorrelationsNew2[[1]][,1:matrix2Length]

permCorrelationsCombine[[1]][,matrix1Length+1:matrix1Length+matrix2Length]


time1 = Sys.time()
time2 = Sys.time()
timeDif = time2 - time1
str(timeDif)
attr(timeDif, "units")


permPval = CategoricalCalculatePermulationPValues(correlationsObject, permCorrelationsNew, report=F)
permPval = CategoricalCalculatePermulationPValues(correlationsObject, combinedPerms, report=F)



?correlateWithCategoricalPhenotype

all.equal(paths, pathsObject)
modelType ="ER"

#categoricalCorrelationOld = categoricalCorrelation

all.equal(categoricalCorrelation[[2]][1], categoricalCorrelationOld[[2]][1])


pairwiseTableNames = names(correlationsObject[[2]])                               #Prepare to repalce the number-number titles with phenotype-phenotype titles
for(i in 1:length(categoryNames)){                                            #for each phenotype
  pairwiseTableNames= gsub(i, names(categoryNames)[i], pairwiseTableNames)                        #replace the number with the phenotype name  
}
names(correlationsObject[[2]]) = pairwiseTableNames                               #update the dataframe titles

?install.packages()


source("Src/Reu/ZoonomTreeNameToCommon.R")
source("Src/Reu/ZonomNameConvertVector.R")
masterTree = mainTrees$masterTree

pdf(height = 25)
ZoonomTreeNameToCommon(readTest, isForegroundTree = T)


dev.off()

mainTips = masterTree$tip.label
commonTips = ZonomNameConvertVectorCommon(mainTips)
scietificTips = ZonomNameConvertVectorCommon(mainTips, common = F)
bothTips = append(scietificTips, commonTips)
write.csv(bothTips, file="Results/SpeciesNames.csv")


CVHRERs = readRDS("Output/CVHRemake/CVHRemakeRERFile.rds")
CVHCorrelation = readRDS("Output/CVHRemake/CVHRemakeCorrelationFile.rds")
CVHCorrelation[order(abs(CVHCorrelation$Rho), decreasing = F),]
source("Src/Reu/RERViolinPlot.R")
mainTrees = readRDS("Data/RemadeTreesAllZoonomiaSpecies.rds")
phenotypeTree = readRDS("Output/CVHRemake/CVHRemakeBinaryForegroundTRee.rds")
fgSpecies = readRDS("Output/CVHRemake/CVHRemakeBinaryTreeForegroundSpecies.rds")
rerViolinPlot()
quickViolin = function(gene){
  rerViolinPlot(mainTrees = mainTrees, CVHRERs, phenotypeTree = phenotypeTree, foregroundSpecies = foregroundSpecies, geneOfInterest = gene, foregroundName = "Carnivore", backgroundName = "Herbivore", backgroundColor = "darkgreen")
}
quickViolin("SLC47A1")
quickViolin("SPEGNB")
quickViolin("SPECC1L")
quickViolin("CRB1")


quickViolin("PROM2")


source("Src/Reu/ZoonomTreeNameToCommon.R")
pdf(height = 25)
ZoonomTreeNameToCommon(binaryForegroundTreeOutput)
dev.off()

source("Src/Reu/RERConvergeFunctions.R")

realCors = correlationsObject
intermediateList = permCorrelations
gene = 1
testOut = list(res = realCors, pvals = list(corsMatPvals, Ppvals), effsize = list(corsMatEffSize, Peffsize))
testOutTrim = testOut$res


categoricalPerms = readRDS("Output/CategoricalDiet3Phen/CategoricalDiet3PhenPermulationsPValueCorrelations.rds")
categoricalNoPerms = readRDS("Output/CategoricalDiet3Phen/CategoricalDiet3PhenCombinedCategoricalCorrelationFile.rds")
categoricalPairwise = readRDS("Output/CategoricalDiet3Phen/CategoricalDiet3PhenPairwiseCorrelationFile.rds")


                            #select the group of pairwise comparisons

phenotypeVectorFilename = paste(outputFolderName, filePrefix, "CategoricalPhenotypeVector.rds",sep="") #select the phenotype vector based on prefix
phenotypeVector = readRDS(phenotypeVectorFilename)                            #load in the phenotype vector 
categories = map_to_state_space(phenotypeVector)                              #and use it to connect branch lengths to phenotype name
categoryNames = categories$name2index                                         #store the length-phenotype connection

pairwiseTableNames = names(categoricalPerms[[2]])                               #Prepare to replace the number-number titles with phenotype-phenotype titles
for(i in 1:length(categoryNames)){                                            #for each phenotype
  pairwiseTableNames= gsub(i, names(categoryNames)[i], pairwiseTableNames)    #replace the number with the phenotype name  
}
names(categoricalPerms[[2]]) = pairwiseTableNames                               #update the dataframe titles

permulationsPValuesFilename = paste(outputFolderName, filePrefix, "PermulationsPValueCorrelations.rds", sep= "")
saveRDS(categoricalPerms, permulationsPValuesFilename)

permulationsPValuesOutput = categoricalPerms

permulationsPValuesOverallFilename = paste(outputFolderName, filePrefix, "PermulationsOverallCorrelations.rds", sep= "")
saveRDS(permulationsPValuesOutput[[1]], permulationsPValuesOverallFilename)

for(i in 1:length(pairwiseTableNames)){
  pairwiseTableNames= gsub(" ", "", pairwiseTableNames)
  permulationsPValuesPairFilename = paste(outputFolderName, filePrefix, "PermulationsCorrelations", pairwiseTableNames[i],".rds", sep= "")
  saveRDS(permulationsPValuesOutput[[2]][i], permulationsPValuesPairFilename)
}

#---

outputSubdirectoryNoslash = paste(outputFolderName, "Overall", sep = "")
if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
  dir.create(outputSubdirectoryNoslash)
}
outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")

permulationsPValuesOverallFilename = paste(outputSubdirectory, filePrefix, "PermulationsOverallCorrelations.rds", sep= "")
saveRDS(permulationsPValuesOutput[[1]], permulationsPValuesOverallFilename)

for(i in 1:length(pairwiseTableNames)){
  pairwiseTableNames= gsub(" ", "", pairwiseTableNames)
  
  outputSubdirectoryNoslash = paste(outputFolderName, pairwiseTableNames[i], sep = "")
  if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
    dir.create(outputSubdirectoryNoslash)
  }
  outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
  
  permulationsPValuesPairFilename = paste(outputSubdirectory, filePrefix, "PermulationsCorrelations", pairwiseTableNames[i],".rds", sep= "")
  saveRDS(permulationsPValuesOutput[[2]][i], permulationsPValuesPairFilename)
}


nonPermCorrelations = readRDS("Output/CategoricalDiet3Phen/CategoricalDiet3PhenCombinedCategoricalCorrelationFile.rds")


nonPermCorrelations[[1]]$Rho

permulationsPValuesOutput=categoricalPerms

i=1

for(i in 1:length(pairwiseTableNames)){
  pairwiseTableNames= gsub(" ", "", pairwiseTableNames)
  
  outputSubdirectoryNoslash = paste(outputFolderName, pairwiseTableNames[i], sep = "")
  if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
    dir.create(outputSubdirectoryNoslash)
  }
  outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
  
  permulationsPValuesPairFilename = paste(outputSubdirectory, filePrefix, pairwiseTableNames[i], "PermulationsCorrelationFile",".rds", sep= "")
  saveRDS(permulationsPValuesOutput[[2]][[i]], permulationsPValuesPairFilename)
}

test1 = permulationsPValuesOutput[[2]][[i]]


correlData$p.adj  

nonPermCorrelations[[2]]$`1 - 3`$p.adj
#saveRDS(mainTrees, "Data/FirstExpressionTrees.rds")
write.csv(mainTrees$masterTree$tip.label, "test.csv")

which(!is.na(RERObject[,1]))

mainTrees$masterTree$tip.label %in% speciesFilter


RERObjectNew = getAllResiduals(mainTrees, useSpecies = speciesFilter, plot = F, min.sp = 3, transform = "none", maxT = 5)

all(is.na(RERObjectNew))



?getAllResiduals()

function (treesObj, cutoff = NULL, transform = "sqrt", weighted = T, 
          useSpecies = NULL, min.sp = 10, scale = T, doOnly = NULL, 
          maxT = NULL, scaleForPproj = F, mean.trim = 0.05, plot = T) 
{
  if (is.null(cutoff)) {
    cutoff = quantile(treesObj$paths, 0.05, na.rm = T)
    message(paste("cutoff is set to", cutoff))
  }
  if (weighted) {
    weights = computeWeightsAllVar(treesObj$paths, transform = transform, 
                                   plot = plot)
    residfunc = fastLmResidMatWeighted
  }
  else {
    residfunc = fastLmResidMat
  }
  if (is.null(useSpecies)) {
    useSpecies = treesObj$masterTree$tip.label
  }
  if (is.null(maxT)) {
    maxT = treesObj$numTrees
  }
  if (transform != "none") {
    transform = match.arg(transform, c("sqrt", "log"))
    transform = get(transform)
  }
  else {
    transform = NULL
  }
  cm = intersect(treesObj$masterTree$tip.label, useSpecies)
  sp.miss = setdiff(treesObj$masterTree$tip.label, useSpecies)
  if (length(sp.miss) > 0) {
    message(paste0("Species from master tree not present in useSpecies: ", 
                   paste(sp.miss, collapse = ",")))
  }
  rr = matrix(nrow = nrow(treesObj$paths), ncol = ncol(treesObj$paths))
  maxn = rowSums(treesObj$report[, cm])
  if (is.null(doOnly)) {
    doOnly = 1
  }
  else {
    maxT = 1
  }
  skipped = double(nrow(rr))
  skipped[] = 0
  for (i in doOnly:(doOnly + maxT - 1)) {
    if (sum(!is.na(rr[i, ])) == 0 && !skipped[i] == 1) {
      tree1 = treesObj$trees[[i]]
      both = intersect(tree1$tip.label, cm)
      if (length(both) < min.sp) {
        next
      }
      tree1 = unroot(pruneTree(tree1, both))
      allreport = treesObj$report[, both]
      ss = rowSums(allreport)
      iiboth = which(ss == length(both))
      if (length(iiboth) < 2) {
        message(paste("Skipping i =", i, "(no other genes with same species set)"))
        next
      }
      nb = length(both)
      ai = which(maxn[iiboth] == nb)
      message(paste("i=", i))
      if (T) {
        ee = edgeIndexRelativeMaster(tree1, treesObj$masterTree)
        ii = treesObj$matIndex[ee[, c(2, 1)]]
        allbranch = treesObj$paths[iiboth, ii]
        if (is.null(dim(allbranch))) {
          message(paste("Issue with gettiing paths for genes with same species as tree", 
                        i))
          return(list(iiboth = iiboth, ii = ii))
        }
        if (weighted) {
          allbranchw = weights[iiboth, ii]
        }
        if (scaleForPproj) {
          nv = apply(scaleMatMean(allbranch), 2, mean, 
                     na.rm = T, trim = mean.trim)
        }
        else {
          nv = apply(allbranch, 2, mean, na.rm = T, trim = mean.trim)
        }
        iibad = which(allbranch < cutoff)
        if (!is.null(transform)) {
          nv = transform(nv)
          allbranch = transform(allbranch)
        }
        allbranch[iibad] = NA
        if (!scale) {
          if (!weighted) {
            proj = residfunc(allbranch[ai, , drop = F], 
                             model.matrix(~1 + nv))
          }
          else {
            proj = residfunc(allbranch[ai, , drop = F], 
                             model.matrix(~1 + nv), allbranchw[ai, , 
                                                               drop = F])
          }
        }
        else {
          if (!weighted) {
            proj = residfunc(allbranch[, , drop = F], 
                             model.matrix(~1 + nv))
          }
          else {
            proj = residfunc(allbranch[, , drop = F], 
                             model.matrix(~1 + nv), allbranchw)
          }
          proj = scale(proj, center = F)[ai, , drop = F]
        }
        rr[iiboth[ai], ii] = proj
      }
    }
  }
  message("Naming rows and columns of RER matrix")
  rownames(rr) = names(treesObj$trees)
  colnames(rr) = namePathsWSpecies(treesObj$masterTree)
  rr
}

correlData[order(correlData$p.adj),]

upham = ReadTntTree("Data/PrunedUphamTree65LiverTPMSpecies.tre")
??Tnt

library(TreeTools)

upham = read.tree("Data/PrunedUphamTree65LiverTPMSpecies.tre")

masterTree = mainTrees$masterTree
plotTree(upham)
plotTree(masterTree)

unrootupham = UnrootTree(upham)

#saveRDS(mainTrees, "Data/NoSignFirstExpressionTrees.rds")

mainTrees$masterTree = unrootupham
#saveRDS(mainTrees, "Data/NoSignFirstExpressionTreesNewMaster.rds")

oldRERs = readRDS("Output/CVHRemake/CVHRemakeRERFile.rds")


cat4Crrels = readRDS("Output/CategoricalDiet4Phen/CategoricalDiet4PhenPairwiseCorrelationFile.rds")
pairwiseTableNames = names(cat4Crrels)

cat5Crrels = readRDS("Output/CategoricalDiet5Phen/CategoricalDiet5PhenPairwiseCorrelationFile.rds")
pairwiseTableNames = names(cat5Crrels)

cat4Crrels[[1]]

for(i in 1:length(pairwiseTableNames)){
  pairwiseTableNames= gsub(" ", "", pairwiseTableNames)
  
  outputSubdirectoryNoslash = paste(outputFolderName, pairwiseTableNames[i], sep = "")
  if(!dir.exists(outputSubdirectoryNoslash)){                       #create that directory if it does not exist
    dir.create(outputSubdirectoryNoslash)
  }
  outputSubdirectory = paste(outputSubdirectoryNoslash, "/", sep="")
  
  permulationsPValuesPairFilename = paste(outputSubdirectory, filePrefix, pairwiseTableNames[i], "CorrelationFile",".rds", sep= "")
  saveRDS(cat5Crrels[[i]], permulationsPValuesPairFilename)
}
j=1

correlation[order(correlation$p.adj),]

saveRDS(mainTrees, "Data/NoSignExpressionTreesRound3.rds")

signmulitmat = readRDS("Data/SignMultiplierExpressionTreesRound3.rds")

unsignedRER = RERObject
signedRER = unsignedRER * signmulitmat

RERObject = signedRER

all.equal(signedRER, unsignedRER)
all.equal(signedRER, RERObject)


# ----------# 


foreground = permulatedForeground; treesObj = mainTrees; plotTree=F; clade="all"; transition="bidirectional"; useSpecies=speciesFilter; weighted = F; 

source("Src/Reu/RERCOnvergeFunctions.R")
foreground2TreeDebug = function (foreground, treesObj, plotTree = T, clade = c("ancestral", "terminal", "all"), weighted = F, transition = "unidirectional", useSpecies = NULL) 
{
  res = treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(res$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      message(paste0("Species from master tree not present in useSpecies: ", 
                     paste(sp.miss, collapse = ",")))
    }
    useSpecies = intersect(useSpecies, res$tip.label)
    res = pruneTree(res, useSpecies)
  }
  foreground = intersect(foreground, useSpecies)
  res$edge.length <- rep(0, length(res$edge.length))
  {
    if (transition == "bidirectional") {
      res <- inferBidirectionalForegroundCladesDebug(res, foreground, 
                                                ancestralOnly = F)
    }
  }
  res
}


treeinput = res; foreground = foreground; ancestralOnly = F; i=1

inferBidirectionalForegroundCladesDebug = function(treeinput, foreground = NULL, ancestralOnly = F){
  tree <- treeinput
  tip.vals=rep(0, length(tree$tip.label))
  names(tip.vals)=tree$tip.label
  tip.vals[foreground]=1
  tmp=cbind(as.character(tip.vals))
  rownames(tmp)=names(tip.vals)
  tip.vals=tmp
  #Add option to function for "type" within ancestral.pars
  ancres=ancestral.pars(tree, df<-as.phyDat(tip.vals, type="USER", levels=unique(as.character(tip.vals))),type="ACCTRAN" )
  ancres=unlist(lapply(ancres, function(x){x[2]}))
  internalVals=ancres
  #evals=matrix(nrow=nrow(treesObj$masterTree$edge), ncol=2)
  evals=matrix(nrow=nrow(tree$edge), ncol=2)
  eres=ancres
  #evals[,1]=eres[treesObj$masterTree$edge[,1]]
  evals[,1]=eres[tree$edge[,1]]
  #evals[,2]=eres[treesObj$masterTree$edge[,2]]
  evals[,2]=eres[tree$edge[,2]]
  tree$edge.length=evals[,2]-evals[,1]
  #res$edge.length[res$edge.length<1]=0
  if(!ancestralOnly){
    edgeIndex=which(tree$edge.length>0)
    edgeIndexNeg=which(tree$edge.length<0)
    edgeIndexAll = c(edgeIndex,edgeIndexNeg)
    edgeDirection = c(rep(1, length(edgeIndex)),rep(-1, length(edgeIndexNeg)))
    edgedf = data.frame(edgeIndexAll,edgeDirection)
    edgedf = edgedf[order(edgedf$edgeIndexAll),]
    clade.edges=NA
    clade.lengths=NA
    cladedf = data.frame(clade.edges,clade.lengths)
    for(i in 1:nrow(edgedf)) { #Does this go from ancestral to terminal?
      #save the clade until the edges no longer overlap
      clade.edges=getAllCladeEdgesDebug(tree, edgedf$edgeIndexAll[i])
      clade.edges=unique(c(edgedf$edgeIndexAll[i], clade.edges))
      if (any(clade.edges %in% cladedf$clade.edges)==F) {
        tree$edge.length[cladedf$clade.edges[which(cladedf$clade.lengths==1)]]=1
        if (edgedf$edgeDirection[i] == 1) {
          clade.lengths = c(rep(1,length(clade.edges)))
        } else {
          clade.lengths = c(rep(0,length(clade.edges)))
        }
        cladedf = data.frame(clade.edges,clade.lengths)
      } else {
        #update df lengths
        if (edgedf$edgeDirection[i] == 1) {
          cladedf$clade.lengths[which(cladedf$clade.edges %in% clade.edges)] = 1
        } else {
          cladedf$clade.lengths[which(cladedf$clade.edges %in% clade.edges)] = 0
        }
      }
    }
    #update edge lengths from the final clade
    tree$edge.length[cladedf$clade.edges[which(cladedf$clade.lengths==1)]]=1
    tree$edge.length[tree$edge.length<0]=0
  }
  tree$edge.length[tree$edge.length<0]=0
  tree
}

tree= tree; AncEdge = edgedf$edgeIndexAll[i]

getAllCladeEdgesDebug = function(tree, AncEdge){
  node=tree$edge[AncEdge,2]
  #get descendants
  iid=getDescendantsDebug(tree, node)
  #find their edges
  iim=match(iid, tree$edge[,2])
  iim
}

tree= tree; node = node; curr = NULL

getDescendantsDebug = function (tree, node, curr = NULL) 
{
  
  daughters <- tree$edge[which(tree$edge[, 1] == node), 2]
  curr <- daughters
  if(is.na(length(curr)) || is.null(length(curr)) || is.na(Ntip(tree)) || is.null(Ntip(tree)) || is.na(node) ||is.null(node)){
    print(node)
    print(curr)
    print(tree)
  }
  if (length(curr) == 0 && node <= Ntip(tree)) 
    curr <- node
  w <- which(daughters > Ntip(tree))
  if (length(w) > 0) 
    for (i in 1:length(w)) curr <- getDescendants(tree, daughters[w[i]], 
                                                  curr)
  return(curr)
}


readRDS("Output/LiverExpression3/LiverExpression3CombinedPrunedFastAppendedPermulationsPValue.rds")

permVal = readRDS("Output/LiverExpression3/LiverExpression3CombinedPrunedFastAppendedPermulationsPValue.rds")
permValNew = readRDS("Output/LiverExpression3/LiverExpression3CombinedPrunedFastAppendedPermulationsPValue.rds")

which(is.null(names(permVal)))
tail(names(permVal))
cleanedPermVal = a
  
which(!names(permVal) %in% unique(names(permVal)))

rownames(correlData)
which(!names(permVal) %in% rownames(correlData))

length(rownames(correlData))
length(names(permVal))

which(!rownames(correlData) %in% names(permVal))

cleanedPermVals = unique(permVal)

length(unique(names(permVal)))

which(duplicated(names(permVal)))
length(which(duplicated(names(permVal))))
4250-3401

permSub1 = permVal[1:850]
permSub2 = permVal[3401:4250]

names(permSub2)

all(names(permSub1) %in% names(permSub2))

all.equal(permSub1, permSub2)

permValClean = permVal[-c(3401:4250)]

all.equal(rownames(correlData), names(permValClean))
match(rownames(correlData), names(permValClean), nomatch = T)

names = data.frame(rownames(correlData)[c(1:17000)], names(permValNew))

unMatched = subset(names, names[[1]] != names[[2]])

all.equal(rownames(correlData), names(permValNew))

tail(rownames(unMatched))


pValset1= readRDS("Output/LiverExpression3/LiverExpression3CombinedPrunedFast1-3400PermulationsPValue.rds")
pValset2= readRDS("Output/LiverExpression3/LiverExpression3CombinedPrunedFast3401-6800PermulationsPValue.rds")
pValset3= readRDS("Output/LiverExpression3/LiverExpression3CombinedPrunedFast6801-10200PermulationsPValue.rds")
pValset4= readRDS("Output/LiverExpression3/LiverExpression3CombinedPrunedFast10201-13600PermulationsPValue.rds")
pValset5= readRDS("Output/LiverExpression3/LiverExpression3CombinedPrunedFast13601-17000PermulationsPValue.rds")


pvalCombined = append(pValset1, pValset2)
pvalCombined = append(pvalCombined, pValset3)
pvalCombined = append(pvalCombined, pValset4)
pvalCombined = append(pvalCombined, pValset5)

names = data.frame(rownames(correlData)[c(1:17000)], names(pvalCombined))
all.equal(rownames(correlData), names(pvalCombined))
unMatched = subset(names, names[[1]] != names[[2]])

saveRDS(pvalCombined, file = "Output/LiverExpression3/LiverExpression3CombinedPrunedFastAllPermulationsPValue.rds")
correlData = correlData[1:17000,]

#correlDataBackup = correlData
correlData = correlDataBackup

correlData = correlData[-which(correlData$permPValue ==0),]
correlDataPositive = correlDataPositive[which(correlDataPositive$permPValue==0)]


saveRDS(correlData, "Output/LiverExpression3/LiverExpression3PermulatedCorrelations.rds")
write.csv(correlData, "Output/LiverExpression3/LiverExpression3PermulatedCorrelations.CSV")
