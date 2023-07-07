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
ZoonomTreeNameToCommon(masterTree)

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
