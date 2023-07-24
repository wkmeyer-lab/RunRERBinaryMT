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
