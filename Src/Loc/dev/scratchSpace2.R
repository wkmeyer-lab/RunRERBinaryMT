a = b #This is to prevent accidental ful runs 

permdata = readRDS("Output/CategoricalDiet3Phen/CategoricalDiet3PhenPermulationsData1.rds")

permdtaaTrim = permdata

permdtaaTrim[[2]] = permdtaaTrim[[2]][1:10]
permdtaaTrim[[3]] = permdtaaTrim[[3]][1:10]


saveRDS(permdtaaTrim, "Output/CategoricalDiet3Phen/CategoricalDiet3PhenPermulationsDataDev.rds")
rm(permdata)
rm(permdtaaTrim)
