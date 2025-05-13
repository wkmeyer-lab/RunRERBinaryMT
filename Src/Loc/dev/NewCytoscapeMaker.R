#------ Run driver assessment code if not already run  --- 
mainPrefix = "CategoricalInsVertivoreTree"
subdirectory = "Herbivore-Vertivore"
BinaryTreeOne = "CategoricalBinaryHerbivoreTree"
BinaryPhenotype = "Herbivore"
BinaryTreeTwo = "CategoricalBinaryInsectivoreTree"
BinaryPhenotypeTwo = "Insectivore"
GOGroup = "KeggReactome"



driverTable = AssessRERDirection(mainPrefix,subdirectory , BinaryTreeOne, BinaryPhenotype,  BinaryTreeTwo, BinaryPhenotypeTwo)
driverTableFilename = paste0("Output/", mainPrefix, "/", subdirectory, "/", mainPrefix, subdirectory, "DriverTable")
write.csv(driverTable, paste0(driverTableFilename, ".csv"))
saveRDS(driverTable, paste0(driverTableFilename, ".rds"))


GODriver = AssessGoCategoryDirection(mainPrefix, subdirectory, GOGroup, 0.1, 1, F)
View(GODriver)

goDriverTableFilename = paste0("Output/", mainPrefix, "/", subdirectory, "/", mainPrefix, subdirectory, "GoDriverTable")
write.csv(GODriver, paste0(goDriverTableFilename, ".csv"))


#--- split GO Driver into positive and negative
goDriverTableFilename = paste0("Output/", mainPrefix, "/", subdirectory, "/", mainPrefix, subdirectory, "GoDriverTable")
GODriver = read.csv(paste0(goDriverTableFilename, ".csv"))


GoDriverPositive = GODriver[which(GODriver$stat > 0),]
GoDriverNegative = GODriver[which(GODriver$stat < 0),]

GODriverPositiveColored = GODriver
GODriverPositiveColored$p.adj[GODriverPositiveColored$stat < 0] = 0.11
GODriverPositiveColored$pval[GODriverPositiveColored$stat < 0] = 0.11


cytoscapeDirectory = paste0("Output/", mainPrefix, "/", subdirectory, "/", "Cytoscape")
if(!dir.exists(cytoscapeDirectory)){                       #create that directory if it does not exist
  dir.create(cytoscapeDirectory)
}


goPositiveGoDriverTableFilename = paste0("Output/", mainPrefix, "/", subdirectory, "/", "Cytoscape/", mainPrefix, subdirectory, "GoDriverPositiveTable")
goNegativeGoDriverTableFilename = paste0("Output/", mainPrefix, "/", subdirectory, "/", "Cytoscape/", mainPrefix, subdirectory, "GoDriverNegativeTable")
goColoredGoDriverTableFilename = paste0("Output/", mainPrefix, "/", subdirectory, "/", "Cytoscape/", mainPrefix, subdirectory, "GoDriverPositiveColoredTable")

write.table(GODriverPositiveColored, paste0(goColoredGoDriverTableFilename, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.csv(GoDriverPositive, paste0(goPositiveGoDriverTableFilename, ".csv"), row.names = F)
write.csv(GoDriverNegative, paste0(goNegativeGoDriverTableFilename, ".csv"), row.names = F)
write.csv(GODriverPositiveColored, paste0(goColoredGoDriverTableFilename, ".csv"), row.names = F)



# ---- Convert a Driver table to cytoscape format 

GODriver

GOCytoscape = GODriver[,c(1,7,3,4,2,6)]
colnames(GOCytoscape) = c("GO.ID", "Description", "p.adj", "DriverPackagedAsQval", "Phenotype", "Gene.vals")
GOCytoscape$Phenotype = sign(GOCytoscape$Phenotype)
GOCytoscape$Phenotype[GOCytoscape$Phenotype == 1] = "+1"
write.table(GOCytoscape, paste0(cytoscapeDirectory, "/CytoscapeInput.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

gmtfile = readRDS(paste0("Data/", GOGroup, ".rds"))
saveRDS(gmtfile, paste0(cytoscapeDirectory, "/", GOGroup, ".rds"))                

readLines("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/cytoscape/KeggReactome.gmt")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GSEABase")
library(GSEABase)
gmtData = getGmt("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/cytoscape/KeggReactome.gmt")

which(names(gmtData) %in% GoDriverPositive$X)
gmtList = as.list(gmtData)

gmtData[which(names(gmtData) %in% GoDriverPositive$X)]


gmtNames = names(gmtData)
gmtDriver = rep(NA, length(gmtNames))
gmtUpdate = data.frame(gmtNames, gmtDriver)
gmtDirections = GODriver$Driver[match(gmtNames, GODriver$X)]
write.csv(gmtDirections, "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/cytoscape/gmtDirectionColumn.csv")
