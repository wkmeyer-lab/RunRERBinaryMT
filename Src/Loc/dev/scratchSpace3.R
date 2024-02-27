

plotRers(rermat = RERObject, index = "OR10J5", phenv = pathsObject)





rownames(enrichmentReorder)

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


testSubfolderData = readRDS("Output/CategoricalDiet3Phen/Overall/CategoricalDiet3PhenOverallPermulationsCorrelations.rds")


?plotRers

plotRers(RERObject, "ALOX15", phenv = pathsObject)

# --------------------------

length(phenotypeVector)
length(mainTrees$masterTree$tip.label)

removedSpecies = mainTrees$masterTree$tip.label[which(!mainTrees$masterTree$tip.label %in% speciesNames)]


droppedspecies = manualAnnots[[annotColumn]][which(manualAnnots$FaName %in% removedSpecies)]
names(droppedspecies) = manualAnnots$FaName[which(manualAnnots$FaName %in% removedSpecies)]
droppedspecies

which(manualAnnots$FaName %in% removedSpecies)



correlation[order(correlation$p.adj),]


?getPermsBinary
getPermsBinary
function (trees, mastertree, root_sp, fg_vec, sisters_list = NULL, 
          pathvec, plotTreeBool = F) 
{
  tip.labels = mastertree$tip.label
  res = getForegroundInfoClades(fg_vec, sisters_list, trees, 
                                plotTree = F, useSpecies = tip.labels)
  fg_tree = res$tree
  fg.table = res$fg.sisters.table
  t = root.phylo(trees$masterTree, root_sp, resolve.root = T)
  rm = ratematrix(t, pathvec)
  if (!is.null(sisters_list)) {
    fg_tree_info = getBinaryPermulationInputsFromTree(fg_tree)
    num_tip_sisters_true = unlist(fg_tree_info$sisters_list)
    num_tip_sisters_true = num_tip_sisters_true[which(num_tip_sisters_true %in% 
                                                        tip.labels)]
    num_tip_sisters_true = length(num_tip_sisters_true)
    fg_tree_depth_order = getDepthOrder(fg_tree)
  }
  else {
    fg_tree_depth_order = NULL
  }
  fgnum = length(which(fg_tree$edge.length == 1))
  if (!is.null(sisters_list)) {
    internal = nrow(fg.table)
  }
  else {
    internal = 0
  }
  tips = fgnum - internal
  testcondition = FALSE
  while (!testcondition) {
    blsum = 0
    while (blsum != fgnum) {
      sims = sim.char(t, rm, nsim = 1)
      nam = rownames(sims)
      s = as.data.frame(sims)
      simulatedvec = s[, 1]
      names(simulatedvec) = nam
      top = names(sort(simulatedvec, decreasing = TRUE))[1:tips]
      t_iter = foreground2Tree(top, trees, clade = "all", 
                               plotTree = F)
      blsum = sum(t_iter$edge.length)
    }
    t_info = getBinaryPermulationInputsFromTree(t_iter)
    if (!is.null(sisters_list)) {
      num_tip_sisters_fake = unlist(t_info$sisters_list)
      num_tip_sisters_fake = num_tip_sisters_fake[which(num_tip_sisters_fake %in% 
                                                          tip.labels)]
      num_tip_sisters_fake = length(num_tip_sisters_fake)
      t_depth_order = getDepthOrder(t_iter)
      testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order)) && 
        (num_tip_sisters_fake == num_tip_sisters_true)
    }
    else {
      t_depth_order = getDepthOrder(t_iter)
      testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order))
    }
  }
  if (plotTreeBool) {
    plot(t_iter)
  }
  return(t_iter)
}

function (tree, numperms, trees, root_sp, fg_vec, sisters_list, 
          pathvec, permmode = "cc") 
{
  if (permmode == "cc") {
    tree_rep = lapply(1:numperms, rep_tree, tree = trees)
    permulated.binphens = lapply(tree_rep, simBinPhenoCC, 
                                 mastertree = trees$masterTree, root_sp = root_sp, 
                                 fg_vec = fg_vec, sisters_list = sisters_list, pathvec = pathvec, 
                                 plotTreeBool = F)
  }
  else if (permmode == "ssm") {
    tree_rep = lapply(1:numperms, rep_tree, tree = tree)
    permulated.binphens = lapply(tree_rep, simBinPhenoSSM, 
                                 trees = trees, root_sp = root_sp, fg_vec = fg_vec, 
                                 sisters_list = sisters_list, pathvec = pathvec)
  }
  else {
    stop("Invalid binary permulation mode.")
  }
  output.list <- list()
  output.list[[1]] <- permulated.binphens
  return(output.list)
}


function (numperms, fg_vec, sisters_list, root_sp, RERmat, trees, 
          mastertree, permmode = "cc", method = "k", min.pos = 2, trees_list = NULL, 
          calculateenrich = F, annotlist = NULL) 
{
  pathvec = foreground2Paths(fg_vec, trees, clade = "all", 
                             plotTree = F)
  col_labels = colnames(trees$paths)
  names(pathvec) = col_labels
  if (permmode == "cc") {
    print("Running CC permulation")
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhen(trees$masterTree, 
                                                    numperms, trees, root_sp, fg_vec, sisters_list, pathvec, 
                                                    permmode = "cc")
    permulated.fg = mapply(getForegroundsFromBinaryTree, 
                           permulated.binphens[[1]])
    permulated.fg.list = as.list(data.frame(permulated.fg))
    phenvec.table = mapply(foreground2Paths, permulated.fg.list, 
                           MoreArgs = list(treesObj = trees, clade = "all"))
    phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[, 
                                                                                  i])
    print("Calculating correlations")
    corMatList = lapply(phenvec.list, correlateWithBinaryPhenotype, 
                        RERmat = RERmat)
    permPvals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permPvals) = rownames(RERmat)
    permRhovals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permRhovals) = rownames(RERmat)
    permStatvals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permStatvals) = rownames(RERmat)
    for (i in 1:length(corMatList)) {
      permPvals[, i] = corMatList[[i]]$P
      permRhovals[, i] = corMatList[[i]]$Rho
      permStatvals[, i] = sign(corMatList[[i]]$Rho) * -log10(corMatList[[i]]$P)
    }
  }
  else if (permmode == "ssm") {
    print("Running SSM permulation")
    if (is.null(trees_list)) {
      trees_list = trees$trees
    }
    RERmat = RERmat[match(names(trees_list), rownames(RERmat)), 
    ]
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhenSSMBatched(trees_list, 
                                                              numperms, trees, root_sp, fg_vec, sisters_list, pathvec)
    df.list = lapply(trees_list, getSpeciesMembershipStats, 
                     masterTree = mastertree, foregrounds = fg_vec)
    df.converted = data.frame(matrix(unlist(df.list), nrow = length(df.list), 
                                     byrow = T), stringsAsFactors = FALSE)
    attr = attributes(df.list[[1]])
    col_names = attr$names
    attr2 = attributes(df.list)
    row_names = attr2$names
    colnames(df.converted) = col_names
    rownames(df.converted) = row_names
    df.converted$num.fg = as.integer(df.converted$num.fg)
    df.converted$num.spec = as.integer(df.converted$num.spec)
    spec.members = df.converted$spec.members
    grouped.trees = groupTrees(spec.members)
    ind.unique.trees = grouped.trees$ind.unique.trees
    ind.unique.trees = unlist(ind.unique.trees)
    ind.tree.groups = grouped.trees$ind.tree.groups
    unique.trees = trees_list[ind.unique.trees]
    unique.map.list = mapply(matchAllNodesClades, unique.trees, 
                             MoreArgs = list(treesObj = trees))
    unique.permulated.binphens = permulated.binphens[ind.unique.trees]
    unique.permulated.paths = calculatePermulatedPaths_apply(unique.permulated.binphens, 
                                                             unique.map.list, trees)
    permulated.paths = vector("list", length = length(trees_list))
    for (j in 1:length(permulated.paths)) {
      permulated.paths[[j]] = vector("list", length = numperms)
    }
    for (i in 1:length(unique.permulated.paths)) {
      ind.unique.tree = ind.unique.trees[i]
      ind.tree.group = ind.tree.groups[[i]]
      unique.path = unique.permulated.paths[[i]]
      for (k in 1:length(ind.tree.group)) {
        permulated.paths[[ind.tree.group[k]]] = unique.path
      }
    }
    attributes(permulated.paths)$names = row_names
    print("Calculating correlations")
    RERmat.list = lapply(seq_len(nrow(RERmat[])), function(i) RERmat[i, 
    ])
    corMatList = mapply(calculateCorPermuted, permulated.paths, 
                        RERmat.list)
    permPvals = extractCorResults(corMatList, numperms, mode = "P")
    rownames(permPvals) = names(trees_list)
    permRhovals = extractCorResults(corMatList, numperms, 
                                    mode = "Rho")
    rownames(permRhovals) = names(trees_list)
    permStatvals = sign(permRhovals) * -log10(permPvals)
    rownames(permStatvals) = names(trees_list)
  }
  else {
    stop("Invalid binary permulation mode.")
  }
  if (calculateenrich) {
    realFgtree = foreground2TreeClades(fg_vec, sisters_list, 
                                       trees, plotTree = F)
    realpaths = tree2PathsClades(realFgtree, trees)
    realresults = getAllCor(RERmat, realpaths, method = method, 
                            min.pos = min.pos)
    realstat = sign(realresults$Rho) * -log10(realresults$P)
    names(realstat) = rownames(RERmat)
    realenrich = fastwilcoxGMTall(na.omit(realstat), annotlist, 
                                  outputGeneVals = F)
    groups = length(realenrich)
    c = 1
    while (c <= groups) {
      current = realenrich[[c]]
      realenrich[[c]] = current[order(rownames(current)), 
      ]
      c = c + 1
    }
    permenrichP = vector("list", length(realenrich))
    permenrichStat = vector("list", length(realenrich))
    c = 1
    while (c <= length(realenrich)) {
      newdf = data.frame(matrix(ncol = numperms, nrow = nrow(realenrich[[c]])))
      rownames(newdf) = rownames(realenrich[[c]])
      permenrichP[[c]] = newdf
      permenrichStat[[c]] = newdf
      c = c + 1
    }
    counter = 1
    while (counter <= numperms) {
      stat = permStatvals[, counter]
      names(stat) = rownames(RERmat)
      enrich = fastwilcoxGMTall(na.omit(stat), annotlist, 
                                outputGeneVals = F)
      groups = length(enrich)
      c = 1
      while (c <= groups) {
        current = enrich[[c]]
        enrich[[c]] = current[order(rownames(current)), 
        ]
        enrich[[c]] = enrich[[c]][match(rownames(permenrichP[[c]]), 
                                        rownames(enrich[[c]])), ]
        permenrichP[[c]][, counter] = enrich[[c]]$pval
        permenrichStat[[c]][, counter] = enrich[[c]]$stat
        c = c + 1
      }
      counter = counter + 1
    }
  }
  if (calculateenrich) {
    data = vector("list", 5)
    data[[1]] = permPvals
    data[[2]] = permRhovals
    data[[3]] = permStatvals
    data[[4]] = permenrichP
    data[[5]] = permenrichStat
    names(data) = c("corP", "corRho", "corStat", "enrichP", 
                    "enrichStat")
  }
  else {
    data = vector("list", 3)
    data[[1]] = permPvals
    data[[2]] = permRhovals
    data[[3]] = permStatvals
    names(data) = c("corP", "corRho", "corStat")
  }
  data
}
<bytecode: 0x00000174ff39f338>
  <environment: namespace:RERconverge>
  > generatePermulatedBinPhen
function (tree, numperms, trees, root_sp, fg_vec, sisters_list, 
          pathvec, permmode = "cc") 
{
  if (permmode == "cc") {
    tree_rep = lapply(1:numperms, rep_tree, tree = trees)
    permulated.binphens = lapply(tree_rep, simBinPhenoCC, 
                                 mastertree = trees$masterTree, root_sp = root_sp, 
                                 fg_vec = fg_vec, sisters_list = sisters_list, pathvec = pathvec, 
                                 plotTreeBool = F)
  }
  else if (permmode == "ssm") {
    tree_rep = lapply(1:numperms, rep_tree, tree = tree)
    permulated.binphens = lapply(tree_rep, simBinPhenoSSM, 
                                 trees = trees, root_sp = root_sp, fg_vec = fg_vec, 
                                 sisters_list = sisters_list, pathvec = pathvec)
  }
  else {
    stop("Invalid binary permulation mode.")
  }
  output.list <- list()
  output.list[[1]] <- permulated.binphens
  return(output.list)
}


function (numperms, fg_vec, sisters_list, root_sp, RERmat, trees, 
          mastertree, permmode = "cc", method = "k", min.pos = 2, trees_list = NULL, 
          calculateenrich = F, annotlist = NULL) 
{
  pathvec = foreground2Paths(fg_vec, trees, clade = "all", 
                             plotTree = F)
  col_labels = colnames(trees$paths)
  names(pathvec) = col_labels
  if (permmode == "cc") {
    print("Running CC permulation")
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhen(trees$masterTree, 
                                                    numperms, trees, root_sp, fg_vec, sisters_list, pathvec, 
                                                    permmode = "cc")
    permulated.fg = mapply(getForegroundsFromBinaryTree, 
                           permulated.binphens[[1]])
    permulated.fg.list = as.list(data.frame(permulated.fg))
    phenvec.table = mapply(foreground2Paths, permulated.fg.list, 
                           MoreArgs = list(treesObj = trees, clade = "all"))
    phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[, 
                                                                                  i])
    print("Calculating correlations")
    corMatList = lapply(phenvec.list, correlateWithBinaryPhenotype, 
                        RERmat = RERmat)
    permPvals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permPvals) = rownames(RERmat)
    permRhovals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permRhovals) = rownames(RERmat)
    permStatvals = data.frame(matrix(ncol = numperms, nrow = nrow(RERmat)))
    rownames(permStatvals) = rownames(RERmat)
    for (i in 1:length(corMatList)) {
      permPvals[, i] = corMatList[[i]]$P
      permRhovals[, i] = corMatList[[i]]$Rho
      permStatvals[, i] = sign(corMatList[[i]]$Rho) * -log10(corMatList[[i]]$P)
    }
  }
  else if (permmode == "ssm") {
    print("Running SSM permulation")
    if (is.null(trees_list)) {
      trees_list = trees$trees
    }
    RERmat = RERmat[match(names(trees_list), rownames(RERmat)), 
    ]
    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhenSSMBatched(trees_list, 
                                                              numperms, trees, root_sp, fg_vec, sisters_list, pathvec)
    df.list = lapply(trees_list, getSpeciesMembershipStats, 
                     masterTree = mastertree, foregrounds = fg_vec)
    df.converted = data.frame(matrix(unlist(df.list), nrow = length(df.list), 
                                     byrow = T), stringsAsFactors = FALSE)
    attr = attributes(df.list[[1]])
    col_names = attr$names
    attr2 = attributes(df.list)
    row_names = attr2$names
    colnames(df.converted) = col_names
    rownames(df.converted) = row_names
    df.converted$num.fg = as.integer(df.converted$num.fg)
    df.converted$num.spec = as.integer(df.converted$num.spec)
    spec.members = df.converted$spec.members
    grouped.trees = groupTrees(spec.members)
    ind.unique.trees = grouped.trees$ind.unique.trees
    ind.unique.trees = unlist(ind.unique.trees)
    ind.tree.groups = grouped.trees$ind.tree.groups
    unique.trees = trees_list[ind.unique.trees]
    unique.map.list = mapply(matchAllNodesClades, unique.trees, 
                             MoreArgs = list(treesObj = trees))
    unique.permulated.binphens = permulated.binphens[ind.unique.trees]
    unique.permulated.paths = calculatePermulatedPaths_apply(unique.permulated.binphens, 
                                                             unique.map.list, trees)
    permulated.paths = vector("list", length = length(trees_list))
    for (j in 1:length(permulated.paths)) {
      permulated.paths[[j]] = vector("list", length = numperms)
    }
    for (i in 1:length(unique.permulated.paths)) {
      ind.unique.tree = ind.unique.trees[i]
      ind.tree.group = ind.tree.groups[[i]]
      unique.path = unique.permulated.paths[[i]]
      for (k in 1:length(ind.tree.group)) {
        permulated.paths[[ind.tree.group[k]]] = unique.path
      }
    }
    attributes(permulated.paths)$names = row_names
    print("Calculating correlations")
    RERmat.list = lapply(seq_len(nrow(RERmat[])), function(i) RERmat[i, 
    ])
    corMatList = mapply(calculateCorPermuted, permulated.paths, 
                        RERmat.list)
    permPvals = extractCorResults(corMatList, numperms, mode = "P")
    rownames(permPvals) = names(trees_list)
    permRhovals = extractCorResults(corMatList, numperms, 
                                    mode = "Rho")
    rownames(permRhovals) = names(trees_list)
    permStatvals = sign(permRhovals) * -log10(permPvals)
    rownames(permStatvals) = names(trees_list)
  }
  else {
    stop("Invalid binary permulation mode.")
  }
  if (calculateenrich) {
    realFgtree = foreground2TreeClades(fg_vec, sisters_list, 
                                       trees, plotTree = F)
    realpaths = tree2PathsClades(realFgtree, trees)
    realresults = getAllCor(RERmat, realpaths, method = method, 
                            min.pos = min.pos)
    realstat = sign(realresults$Rho) * -log10(realresults$P)
    names(realstat) = rownames(RERmat)
    realenrich = fastwilcoxGMTall(na.omit(realstat), annotlist, 
                                  outputGeneVals = F)
    groups = length(realenrich)
    c = 1
    while (c <= groups) {
      current = realenrich[[c]]
      realenrich[[c]] = current[order(rownames(current)), 
      ]
      c = c + 1
    }
    permenrichP = vector("list", length(realenrich))
    permenrichStat = vector("list", length(realenrich))
    c = 1
    while (c <= length(realenrich)) {
      newdf = data.frame(matrix(ncol = numperms, nrow = nrow(realenrich[[c]])))
      rownames(newdf) = rownames(realenrich[[c]])
      permenrichP[[c]] = newdf
      permenrichStat[[c]] = newdf
      c = c + 1
    }
    counter = 1
    while (counter <= numperms) {
      stat = permStatvals[, counter]
      names(stat) = rownames(RERmat)
      enrich = fastwilcoxGMTall(na.omit(stat), annotlist, 
                                outputGeneVals = F)
      groups = length(enrich)
      c = 1
      while (c <= groups) {
        current = enrich[[c]]
        enrich[[c]] = current[order(rownames(current)), 
        ]
        enrich[[c]] = enrich[[c]][match(rownames(permenrichP[[c]]), 
                                        rownames(enrich[[c]])), ]
        permenrichP[[c]][, counter] = enrich[[c]]$pval
        permenrichStat[[c]][, counter] = enrich[[c]]$stat
        c = c + 1
      }
      counter = counter + 1
    }
  }
  if (calculateenrich) {
    data = vector("list", 5)
    data[[1]] = permPvals
    data[[2]] = permRhovals
    data[[3]] = permStatvals
    data[[4]] = permenrichP
    data[[5]] = permenrichStat
    names(data) = c("corP", "corRho", "corStat", "enrichP", 
                    "enrichStat")
  }
  else {
    data = vector("list", 3)
    data[[1]] = permPvals
    data[[2]] = permRhovals
    data[[3]] = permStatvals
    names(data) = c("corP", "corRho", "corStat")
  }
  data
}


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
