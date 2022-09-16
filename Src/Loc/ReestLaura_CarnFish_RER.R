#Use existing gene trees, RER, and phenotype tree
#from whole Boreoeutheria set to run RERconverge
#for just Laurasiatheria subset
#Redo with re-estimated branch lengths using gene trees
#Trees are here: Phylotrees_FishCarn_BranchReest.rds
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries
library(RERconverge) #load RERconverge package
setwd("/share/ceph/wym219group/shared/projects/MammalDiet/Zoonomia/CarnFishvHerbs/") #go to directory

#Load existing data
#fcrer <- readRDS("mamRERCMU_FishCarn.rds") #load existing RER
ttree <- read.tree("Culled_0511_FishCarnvHerb.txt") #load existing phenotype tree

#Use spreadsheet to find Laurasiatheria species within tree
ma <- read.csv("MA_Laura.csv") #read in spreadsheet
subma <- ma[intersect(which(ma$Laurasiatheria==1),
		      which(ma$CarnFish_Herbs %in% c(0,1))),] 
#subset data frame to just phenotyped Laurasiatheria species
laurasp <- subma$FaName[subma$FaName %in% ttree$tip.label] #get overlap of Laurasiatheria and species in phenotype tree

#Prune tree if needed
ttree <- pruneTree(ttree, laurasp)

#Re-compute RERs
FishTree <- readRDS('Phylotrees_FishCarn_BranchReest.rds') #Is this just Laurasiatheria/
if (!file.exists("mamRERCMU_FishCarn_Laura_reest.rds")) {
        fcrer <- getAllResiduals(FishTree, useSpecies = laurasp, plot = FALSE)
	saveRDS(fcrer, file="Phylotrees_FishCarn_BranchReest_Laura.rds")
} else {
	fcrer <- readRDS("mamRERCMU_FishCarn_Laura_reest.rds")
}

#Create paths from existing tree, restricting to Laurasiatheria with useSpecies
if (!file.exists("Paths_CarnFish_Laura.rds")) {
	tpaths <- tree2Paths(ttree, FishTree, binarize=T, useSpecies=laurasp)
	saveRDS(tpaths, file="Paths_CarnFish_Laura.rds")
} else {
	tpaths <- readRDS("Paths_CarnFish_Laura.rds")
}
#Correlate phenotype paths with gene RERs
tcor <- correlateWithBinaryPhenotype(fcrer, tpaths, min.sp=35) #run correlation, using min.sp
#(1/2 of all species with phenotypes)

write.csv(tcor, file="Reest_CorrelationFishCarnHerbLaurasiatheria_wminsp.csv", 
	  row.names=T, quote=F) #write out the results
saveRDS(tcor, file= "Reest_Laura_GOinput_wminsp.rds") 
