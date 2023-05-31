a = b #This is to prevent accidental ful runs 

permdata = readRDS("Output/CategoricalDiet3Phen/CategoricalDiet3PhenPermulationsData1.rds")

permdtaaTrim = permdata

permdtaaTrim[[2]] = permdtaaTrim[[2]][1:3]
permdtaaTrim[[3]] = permdtaaTrim[[3]][1:3]


saveRDS(permdtaaTrim, "Output/CategoricalDiet3Phen/CategoricalDiet3PhenPermulationsDataDev.rds")
rm(permdata)
rm(permdtaaTrim)

library(data.table)
source("Src/Reu/CategoricalPermulationsParallelFunctions.R")

permCorrelations = getPermPvalsCategorical3(correlationsObject, permulationsData$trees, phenotypeVector, mainTrees, RERObject)
permCorrelations = getPermPvalsCategorical2(correlationsObject, permulationsData$trees, phenotypeVector, mainTrees, RERObject)

permCorrelations = CategoricalPermulationGetCor(correlationsObject, permulationsData$trees, phenotypeVector, mainTrees, RERObject, report=T)

time1 = Sys.time()
time2 = Sys.time()
timeDif = time2 - time1
str(timeDif)
attr(timeDif, "units")


permPval = CategoricalCalculatePermulationPValues(correlationsObject, permCorrelations)

?correlateWithCategoricalPhenotype

all.equal(paths, pathsObject)
modelType ="ER"

#categoricalCorrelationOld = categoricalCorrelation

all.equal(categoricalCorrelation[[2]][1], categoricalCorrelationOld[[2]][1])
