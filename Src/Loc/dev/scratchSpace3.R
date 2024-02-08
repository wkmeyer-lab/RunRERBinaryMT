phenotypeVector = readRDS("Output/IPCRelaxTest/HillerIPCPhenotype.rds")
speciesFilter = names(phenotypeVector) 
speciesFilter = speciesFilter[-which(speciesFilter == "ornAna2")]


phenotypeVector = readRDS("Output/HMGRelaxTest/HillerHGMPhenotype.rds")
speciesFilter = names(phenotypeVector) 


library(RERconverge)
mainTrees = readRDS("Data/RemadeTreesAllZoonomiaSpecies.rds")
masterTree = mainTrees$masterTree
write.tree(masterTree, "Results/ZoonomiaMasterTree.tree")


list.files("data")

masterTree2 =  ZoonomTreeNameToCommon(masterTree, scientific = T)
write.tree(masterTree2, "Results/zoonomiaMasterTreeScientific.tree")


source("src/reu/ZoonomTreeNameToCommon.R")

Phen5Tree = readRDS("Output/CategoricalDiet5Phen/CategoricalDiet5PhenCategoricalTree.rds")

ZoonomTreeNameToCommon(Phen5Tree)


hillerConversionTable = read.csv("Data/HIllerZoonomPhenotypeTable.csv")

relevantSpecies

hillerNames = match(speciesNames, hillerConversionTable$Zoonomia)
?match()
