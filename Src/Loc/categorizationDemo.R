library(RERconverge)
source("Src/Reu/ZonomNameConvertMatrixCommon.R")
source("Src/Reu/ZonomNameConvertVector.R")

allInsectivoryData = read.csv("Data/allInsectivoryCorrelationFile.csv")

RERs = readRDS("Data/allInsectivoryRERFile.rds")
Paths = readRDS("Data/allInsectivoryPathsFile.rds")
plotRers(RERs,"KIAA0825", Paths )

CNRers = ZonomNameConvertMatrixCommon(RERs)
CNNames = CNRers[1,]
names(CNNames)
bats = grep("bat", names(CNNames))

binaryTree = readRDS("Data/allInsectivoryBinaryForegroundTree.rds")
plotTree(binaryTree)

zonomMaster = readRDS("Data/RemadeTreesAllZoonomiaSpecies.rds")
functionPaths = categorizePaths(binaryTree, zonomMaster, "functionPath", overwrite = T)
plotRers(CNRers,"KIAA0825", functionPaths, sortrers = T)

plotRers(CNRers,"ZNF292", manualPaths, sortrers = T)
plotRers(CNRers,"ZNF292", Paths, sortrers = T)
