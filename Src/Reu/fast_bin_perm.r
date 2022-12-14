#Fast code for binary permulations
#Daniel Schaffer, Pfenning Lab

library(geiger)


#' Sets a leaf to be a descendent of a node in a bitmap
#' @param bits An array of bitmaps represented by a 2d integer array where each column corrsponds to one node
#' @param internal An internal node node such that 1 <= internal <= dim(bits)[2]
#' @param leaf A value representing a leaf node such that 1 <= leaf <= dim(bits)[1]*32
#' @param nInts dim(bits)[1]
#' @export
Rcpp::cppFunction("void setLeafDesc(IntegerVector bits, int internal, int leaf, int nInts) {
    int idx = nInts * (internal - 1) + (leaf - 1) / 32;
    int bit = (leaf - 1) % 32;
    bits[idx] |= 1 << bit;
    //return bits;
} 
")


#' Performs recursive DFS traversal of a tree to populate an array of bitmaps for ancestor-leaf relationships
#' @param desc array of bitmaps represented by a 2d integer array where each column corrsponds to one node and each bit corresponds to one leaf
#' @param node Current node in the tree
#' @param parents Array of the ancestors of node
#' @param parents Adjacency-list representation of the internal nodes of the tree.
#' @param nLeaves the number of leaves of the tree. Leaves are numbered 1-nLeaves; internal nodes are numbered (nLeaves+1)-(2nLeaves-1)
#' @param nInts Number of ints used for each bitmap, should be ceiling(nLeaves / 32)
#' @return Nothing, desc is modified-in-place
#' @export
leafDFS=function(desc, node, parents, adj, nLeaves, nInts){ #I'm not sure what desc is exactly. 
  if (node <= nLeaves){                                  #This checks if the node is a species, or internal. It does this by seeing if the current node is less than the number of leaves (species)
    for (i in 1:length(parents)) {                       #If so, it records the parents of the node. 
      setLeafDesc(desc, parents[i], node, nInts)
    }
  } else{                                               #If not, it's an internal node.                                     
    internalIndex = node - nLeaves                      # The internal index is the number of internal node that you are counting. This is the node's number, minus the number of terminal nodes. So if you have 10 terminal nodes, and you are on node 11, you are on internal node 1 
    parents[length(parents) + 1] = internalIndex
    l = adj[internalIndex,1]
    r = adj[internalIndex,2]
    leafDFS(desc, l, parents, adj, nLeaves, nInts)
    leafDFS(desc, r, parents, adj, nLeaves, nInts)
  }
}

#' Generates an array of bitmaps for ancestor-leaf relationships
#' @param t A fully dichotomous binary tree of class phylo
#' @return array of bitmaps represented by a 2d integer array where each column corrsponds to one node and each bit corresponds to one leaf
#' @export
makeLeafMap=function(t) {
  #initialize
  nInternal = t$Nnode
  nLeaves = nInternal + 1
  nInts = ceiling(nLeaves / 32)
  desc = array(dim=c(nInts, nInternal))
  adj = array(dim=c(nInternal, 2))
  for (i in 1:(nInternal)) {
    #node = toString(i+nLeaves)
    desc[,i] = integer(nInts)
    adj[i,] = integer(2)
  }
  #Compute adjanceny list
  edges = t$edge
  nEdges = length(edges[,1])
  for (i in 1:nEdges){
    par = edges[i,1] - nLeaves
    chi = edges[i,2]
    if (adj[par,1] == 0) {
      adj[par,1] = chi
    } else {
      adj[par,2] = chi
    }
  }
  #Because we cheat and use C to modify, we still get aliasing
  leafDFS(desc, nLeaves + 1, integer(0), adj, nLeaves, nInts)
  desc
}


#' Creates a bitmap representing an array of integers
#' @param indices Some positive integers
#' @param len The largest integer to represent in the bitmap; len >= max(indices)
#' @return An array of integerrs representing a bitmap such that bit (i-1) = 1 iff i is in indices
#' @export
Rcpp::cppFunction("IntegerVector makeBitMap(IntegerVector indices, int len) {
    int nInts = len / 32 + ((len % 32 == 0) ? 0 : 1);
    IntegerVector bits (nInts);
    for(int i=0; i<indices.length(); ++i){
        int leaf = indices[i] - 1;  //Zero-indexed
        int idx = leaf / 32;
        int bit = leaf % 32;
        bits[idx] |= 1 << bit;
    }
    return bits;
} 
")

#' Compares two bitmaps represented by arrays of integers
#' @param bits A bitmap represented by an array of integers
#' @param fgBits A bitmap represented by an array of integers with length(fgBits) >= length(bits) 
#' @return True if bits is a subset of fgBits, i.e. bits[i] = 1 => fgBits[i] = 1
#' @export
Rcpp::cppFunction("bool compare(IntegerVector bits, IntegerVector fgBits) {
    bool res = true;
    for(int i=0; i<bits.length(); ++i){
      res &= ((bits[i] & fgBits[i]) == bits[i]);
    }
    return res;
} 
")

#' Counts the foreground internal nodes of a tree
#' @param t A tree of class phylo, with tip labels
#' @param bitMaps The result of makeLeafMap(t) 
#' @param fg A subset of tip labels of t, the foreground species
#' @return The number of internal nodes in t with all leaf descendents in the foreground
#' @export
countInternal=function(t, bitMaps, fg){
  nInternal = t$Nnode
  fgNums = which(t$tip.label %in% fg)
  fgBitMap = makeBitMap(fgNums, nInternal + 1)
  count = 0
  for (i in 1:nInternal) {
    if (compare(bitMaps[,i], fgBitMap)) {
      count = count +1
    }
  }
  count
}


#This function is modified from the simBinPhenoVec function in RERConverge
#'Generates a permulated phenotype vector whose phylogeny matches a desired structure.  
#'User may specify the number of foreground branches that are internal branches.
#' @param tree A tree of type phylo
#' @param phenvec Named vector of 1's and 0's representing phenotype values for each species in tree
#' @param internal Total number of foreground internal nodes desired
#' @param tips Number of foreground species desired. Inferred from pfenvec if phenvec is specified
#' @param rm Precomputed rate matrix for T. Computed if not specified and phenvec is given.
#' @param leafBitMaps Precomputed result of makeLeafMap(tree)
#' @return A vector of permulated foreground species
fastSimBinPhenoVec=function(tree, tips=0, internal, phenvec = NULL, rm=NULL, leafBitMaps=NULL){
  insum=0
  top = character(0)
  if (!is.null(phenvec)) {
    tips = sum(phenvec)
    if(is.null(rm)) {
      rm = ratematrix(tree, phenvec)
    }
  } else if (is.null(rm)){
    print("No rate matrix supplied & cannot compute one")
    return()
  }
  
  if(is.null(leafBitMaps)){
    leafBitMaps = makeLeafMap(tree)
  }
  
  repeat { #implements a do-while loop
    t=tree 
    #root.phylo(tree, root, resolve.root = T) #For now, let's use the ancestral root
    sims=sim.char(t, rm, nsim = 1)
    nam=rownames(sims)
    s=as.data.frame(sims)
    simulatedvec=s[,1]
    names(simulatedvec)=nam
    top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
    insum=countInternal(tree, leafBitMaps, top)
    if (insum==internal) {break}
  }
  # plot(t)
  return(top)
}

fastSimBinPhenoVecReport=function(tree, tips=0, internal, phenvec = NULL, rm=NULL, leafBitMaps=NULL){
  insum=0
  top = character(0)
  if (!is.null(phenvec)) {
    tips = sum(phenvec)
    if(is.null(rm)) {
      rm = ratematrix(tree, phenvec)
      message(rm)
    }
  } else if (is.null(rm)){
    print("No rate matrix supplied & cannot compute one")
    return()
  }
  
  if(is.null(leafBitMaps)){
    leafBitMaps = makeLeafMap(tree)
  }
  loop=0
  insumList = vector()
  repeat { #implements a do-while loop
    t=tree 
    #root.phylo(tree, root, resolve.root = T) #For now, let's use the ancestral root
    sims=sim.char(t, rm, nsim = 1)
    nam=rownames(sims)
    s=as.data.frame(sims)
    simulatedvec=s[,1]
    names(simulatedvec)=nam
    top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
    insum=countInternal(tree, leafBitMaps, top)
    loop = loop+1
    message("loop: ", loop)
    
    message("insum: ",insum)
    #insumList = append(insumList, insum)
    
    if (insum==internal) {break}
  }
  # plot(t)
  return(top)
}


