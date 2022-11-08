library(RERconverge)

carnVHerbTree = read.tree("Data/Culled_0511_FishCarnvHerb.txt")
plotTreeHighlightBranches(carnVHerbTree, hlspecies = carnVHerbTree$edge.length == 1)

relevantSpeciesNames = carnVHerbTree$tip.label
binaryForegroundTreeOutput = carnVHerbTree
binaryTreeFilename = paste(outputFolderName, filePrefix, "BinaryForegroundTree.rds", sep="")
saveRDS(binaryForegroundTreeOutput, file = binaryTreeFilename)
saveRDS()