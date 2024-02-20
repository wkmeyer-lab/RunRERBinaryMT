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

hillerConversionTableLocation = "Data/HillerZoonomiaPhenotypeTable.csv"
hillerConversionTable = read.csv("Data/HIllerZoonomPhenotypeTable.csv")

relevantSpecies

hillerNames = match(speciesNames, hillerConversionTable$Zoonomia)
?match()

newHillerMasterTree = read.tree("Data/newHillerMasterTree.txt")
plotTree(newHillerMasterTree)
nodelabels(col = "red", adj = c(0, -0), frame = "none")
tiplabels(col = "blue", frame = "none")
edgelabels(col = "darkgreen", frame = "none")


is.binary(newHillerMasterTree)

newHillerTrees = readTrees('Data/newHillerMainTrees2.txt')


testTree = newHillerTrees$trees[[1]]

is.binary(newHillerMasterTree)

is.binary(testTree)
plotTree(testTree)

demoTree= testTree
demoTree$edge.length = rep(1, length = length(demoTree$edge))

plotTree(demoTree)
ZoonomTreeNameToCommon(demoTree, manualAnnotLocation = hillerConversionTableLocation)


lapply(newHillerTrees$trees[[]], is.binary())


oldHilTrees = oldHillerTrees$trees


all(sapply(oldHilTrees, is.binary))


?readTrees

saveRDS(newHillerTrees, "Data/newHillerMainTrees.rds")

oldHillerTrees = readRDS("Data/mam120aa_trees.rds")

newHillerTrees = names(newHillerTrees$trees)
oldHillerTrees = names(oldHillerTrees$trees)

missingGenes = oldHillerTrees[which(oldHillerTrees %in% newHillerTrees)]





# --------------------------

length(phenotypeVector)
length(mainTrees$masterTree$tip.label)

removedSpecies = mainTrees$masterTree$tip.label[which(!mainTrees$masterTree$tip.label %in% speciesNames)]


droppedspecies = manualAnnots[[annotColumn]][which(manualAnnots$FaName %in% removedSpecies)]
names(droppedspecies) = manualAnnots$FaName[which(manualAnnots$FaName %in% removedSpecies)]
droppedspecies

which(manualAnnots$FaName %in% removedSpecies)



correlation[order(correlation$p.adj),]


function (tipvals, treesObj, useSpecies = NULL, model = "ER", 
          root_prior = "auto", plot = FALSE, anctrait = NULL) 
{
  mastertree = treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(mastertree$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      message(paste0("Species from master tree not present in useSpecies: ", 
                     paste(sp.miss, collapse = ",")))
    }
    useSpecies = intersect(mastertree$tip.label, useSpecies)
    mastertree = pruneTree(mastertree, useSpecies)
    mastertree = unroot(mastertree)
  }
  else {
    mastertree = pruneTree(mastertree, intersect(mastertree$tip.label, 
                                                 names(tipvals)))
    mastertree = unroot(mastertree)
  }
  if (is.null(anctrait)) {
    tipvals <- tipvals[mastertree$tip.label]
    intlabels <- map_to_state_space(tipvals)
    print("The integer labels corresponding to each category are:")
    print(intlabels$name2index)
    ancliks = getAncLiks(mastertree, intlabels$mapped_states, 
                         rate_model = model, root_prior = root_prior)
    states = rep(0, nrow(ancliks))
    for (i in 1:length(states)) {
      states[i] = which.max(ancliks[i, ])
    }
    states = c(intlabels$mapped_states, states)
    tree = mastertree
    tree$edge.length = states[tree$edge[, 2]]
    if (length(unique(tipvals)) == 2) {
      if (sum(!unique(tipvals) %in% c(TRUE, FALSE)) > 0) {
        message("Returning categorical tree for binary phenotype because phenotype values are not TRUE/FALSE")
      }
      else {
        tree$edge.length = ifelse(tree$edge.length == 
                                    2, 1, 0)
        print("There are only 2 categories: returning a binary phenotype tree.")
        if (plot) {
          plotTree(tree)
        }
        return(tree)
      }
    }
    if (plot) {
      plotTreeCategorical(tree, category_names = intlabels$state_names, 
                          master = mastertree, node_states = states)
    }
    return(tree)
  }
  else {
    if (length(unique(tipvals)) <= 2) {
      fgspecs <- names(tipvals)[tipvals != anctrait]
      res <- foreground2Tree(fgspecs, treesObj, plotTree = plot, 
                             clade = "terminal", useSpecies = useSpecies)
      print("There are only 2 categories: returning a binary phenotype tree.")
      if (plot) {
        plotTree(res)
      }
      return(res)
    }
    else {
      tipvals <- tipvals[mastertree$tip.label]
      intlabels <- map_to_state_space(tipvals)
      j <- which(intlabels$state_names == anctrait)
      if (length(j) < 1) {
        warning("The ancestral trait provided must match one of the traits in the phenotype vector.")
      }
      res = mastertree
      res$edge.length <- rep(j, length(res$edge.length))
      traits <- intlabels$state_names
      for (trait in traits) {
        if (trait == anctrait) {
          next
        }
        i <- which(intlabels$state_names == trait)
        res$edge.length[nameEdges(res) %in% names(tipvals)[tipvals == 
                                                             trait]] = i
      }
      names(res$edge.length) = nameEdges(res)
      if (plot) {
        states = res$edge.length[order(res$edge[, 2])]
        states = c(j, states)
        plotTreeCategorical(res, category_names = traits, 
                            master = treesObj$masterTree, node_states = states)
      }
      print("Category names are mapped to integers as follows:")
      print(intlabels$name2index)
      return(res)
    }
  }
}
