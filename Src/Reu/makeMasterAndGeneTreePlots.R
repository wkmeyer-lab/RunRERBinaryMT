makeMasterVsGeneTreePlots = function(mainTrees, RERObject, geneInQuestion, foregroundVector, colors = "blue", correlationPlot = F){
  source("Src/Reu/ZonomNameConvertVector.R")
  source("Src/Reu/ZoonomTreeNameToCommon.R")
  source("Src/Reu/GetForegroundEdges.R")
  foregroundVector = ZonomNameConvertVectorCommon(foregroundVector)
  RERObject[geneInQuestion,]
  
  namesToKeep = names(RERObject[geneInQuestion,][!is.na(RERObject[geneInQuestion,])])
  namesToKeep = namesToKeep[!is.na(namesToKeep)]
  
  masterTree = mainTrees$masterTree
  prunedMaster = drop.tip(masterTree, which(!masterTree$tip.label %in% namesToKeep))
  geneTree = mainTrees$trees[[geneInQuestion]]
  prunedGeneTree = drop.tip(geneTree, which(!geneTree$tip.label %in% namesToKeep))
  
  commonMaster = ZoonomTreeNameToCommon(prunedMaster, plot =F)
  commonGene = ZoonomTreeNameToCommon(prunedGeneTree, plot =F)
  
  par(mfrow = c(1,2), mai = c(0.5, 0.1, 0.2, 0.1))
  masterFGEdges = getForegroundEdges(commonMaster, foregroundVector)
  plotTreeHighlightBranches(commonMaster, hlspecies = masterFGEdges, main = "Overall Genome Average", hlcols = colors)
  
  geneFGEdges = getForegroundEdges(commonGene, foregroundVector)
  plotTreeHighlightBranches(commonGene, main = paste("Gene:", geneInQuestion), hlspecies = geneFGEdges, hlcols = colors)
  
  if(correlationPlot){
    if(all.equal(masterFGEdges, geneFGEdges)){ #Only run this code if the trees are the same shape
      
      edgeGround = rep(NA, nrow(edgeCorrealtions))
      edgeGround[masterFGEdges] = 1
      edgeGround[!1:length(edgeGround) %in% masterFGEdges] = 0
      
      endNodes = commonMaster$edge[,2]
      edgeNames = commonMaster$tip.label[endNodes]
      
      edgeCorrealtions = data.frame(commonMaster$edge.length, commonGene$edge.length, as.factor(edgeGround), edgeNames)
      names(edgeCorrealtions) = c("Master", "Gene", "Ground", "Names")
      
      
      
      g <- ggplot(edgeCorrealtions, 
                  aes(x = Master, 
                      y = Gene, 
                      col = Ground, 
                      label = Names)
                  ) + 
                  scale_size_manual(values = c(1, 1, 1, 1)) + 
                  geom_point(aes(size = Ground)) + 
                  scale_color_manual(values = c("black", colors)) + 
                  coord_fixed() + 
                  expand_limits(x = max(edgeCorrealtions$Gene), y = max(edgeCorrealtions$Master))+
                  geom_text(hjust = "bottom", size = 4, check_overlap = T) + 
                  ylab("Gene-Specific Branch Length") + 
                  xlab("Genome Average Branch Length") + 
                  #ggtitle(plottitle) + 
                  geom_abline(intercept = 0, slope =1, linetype = "dotted") + 
                  theme(
                      #axis.ticks.y = element_blank(), 
                      #axis.text.y = element_blank(), 
                      legend.position = "none", 
                      panel.background = element_blank(), 
                      axis.text = element_text(size = 12, colour = "black"), 
                      axis.title = element_text(size = 20, face = "bold"), 
                      plot.title = element_text(size = 24, face = "bold")
                      ) +
                  theme(axis.line = element_line(colour = "black", size = 0.2))  
                  #theme(axis.line.y = element_blank())
      
    }
  }
  
}

