knitr::opts_chunk$set(echo = TRUE)
library(RERconverge)







# --- Walkthrough ----
rerpath = find.package('RERconverge')
#Get main tree
toytreefile = "subsetMammalGeneTrees.txt" 
toyTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)

#-Generate binary path tree-
#Generate binary tree
marineFg = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee")
sisters_marine = list("clade1"=c("Killer_whale", "Dolphin"))
marineFgTree = foreground2TreeClades(marineFg,sisters_marine,toyTrees,plotTree=F)
marineplot1 = plotTreeHighlightBranches(marineFgTree,
                                        hlspecies=which(marineFgTree$edge.length==1),
                                        hlcols="blue", main="Marine mammals trait tree")
#convert to path
pathvec = tree2PathsClades(marineFgTree, toyTrees)
# Calculating paths from the foreground tree
pathvec = tree2PathsClades(marineFgTree, toyTrees)

# Calculate RERs
mamRERw = getAllResiduals(toyTrees, transform="sqrt", weighted=T, scale=T)

# Calculate correlation 
res = correlateWithBinaryPhenotype(mamRERw, pathvec, min.sp=10, min.pos=2,
                                   weighted="auto")
