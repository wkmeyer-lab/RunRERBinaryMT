phenotypeVector = readRDS("Output/IPCRelaxTest/HillerIPCPhenotype.rds")
speciesFilter = names(phenotypeVector) 
speciesFilter = speciesFilter[-which(speciesFilter == "ornAna2")]


phenotypeVector = readRDS("Output/HMGRelaxTest/HillerHGMPhenotype.rds")
speciesFilter = names(phenotypeVector) 
