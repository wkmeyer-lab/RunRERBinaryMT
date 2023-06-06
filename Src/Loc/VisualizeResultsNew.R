# -- Libraries 
.libPaths("/share/ceph/wym219group/shared/libraries/R4") #add path to custom libraries to searched locations
library(RERconverge)
library(ggplot2)
library(ggpmisc)
library(cowplot)
library(qvalue)
library(insight)
source("Src/Reu/cmdArgImport.R")

# -- Usage:
# This script generates several visualizations of the RER outputs, and then organizes them into a single PDF file. 

# -- Command arguments list
# r = filePrefix                                     This is a prefix used to organize and separate files by analysis run. Always required. 
# v = <T or F>                                       This prefix is used to force the regeneration of the script's output, even if the files already exist. Not required, not always used.
# g = <T OR F>                                       This sets if Gene Enrichment analysis is included. Default TRUE. 
# l = <"enrichmentListName" OR c("name1","name2")>   This defines the enrichment list(s) used in the GO analysis the script should display

# p = <T or F or B>                                  This sets if the code should use permulated or unpermualted values, or both.
# f = "permulationPvalueFileLocation.rds"            This is a manual override to specify the script use a specific Permulation p-value file. 
      #If using any permulation p-value file other than "MainCombined" with no run instance number, it must be specified manually.
              


#----------------
args = c('r=CVHRemake', 'p=B') #This is a debug argument set. It is used to set arguments locally, when not running the code through a bash script.

# --- Standard start-up code ---
args = commandArgs(trailingOnly = TRUE)
{  # Bracket used for collapsing purposes
  #File Prefix
  if(!is.na(cmdArgImport('r'))){                                                #This cmdArgImport script is a way to import arguments from the command line. 
    filePrefix = cmdArgImport('r')
  }else{
    stop("THIS IS AN ISSUE MESSAGE; SPECIFY FILE PREFIX")
  }
  
  #  Output Directory 
  if(!dir.exists("Output")){                                      #Make output directory if it does not exist
    dir.create("Output")
  }
  outputFolderNameNoSlash = paste("Output/",filePrefix, sep = "") #Set the prefix sub directory
  if(!dir.exists(outputFolderNameNoSlash)){                       #create that directory if it does not exist
    dir.create(outputFolderNameNoSlash)
  }
  outputFolderName = paste("Output/",filePrefix,"/", sep = "")
  
  #  Force update argument
  forceUpdate = FALSE
  if(!is.na(cmdArgImport('v'))){                                 #Import if update being forced with argument 
    forceUpdate = cmdArgImport('v')
    forceUpdate = as.logical(forceUpdate)
  }else{
    paste("Force update not specified, not forcing update")
  }
}

# --- Argument Imports ---
# Defaults
permulationDefaultFilename = "MainCombinedPermulationsPValue.rds" #This argument isn't actually updated from the command line, it is used to easily swap what the default permulation filename the script expects.
permulationDefaultFilename= "CombinedPrunedFastAllPermulationsPValue.rds"

usePermulations = TRUE
permulationPValOverride = NULL 
useGeneEnrichment = TRUE
useBoth = FALSE

{ # Bracket used for collapsing purposes
  
  #Permulation use
  if(!is.na(cmdArgImport('p'))){
    usePermulations = cmdArgImport('p')
    if(usePermulations == "B"){
      useBoth = TRUE
      usePermulations = TRUE
    }else{
      usePermulations = as.logical(usePermulations)
      if(is.na(usePermulations)){
        usePermulations = TRUE
        useBoth = FALSE
        message("Use Permualtions value not interpretable as logical. Did you remember to capitalize? Using TRUE.")
      }
    }
  }else{
    message("Use Permulations value not specified, using TRUE.")
  }
  
  #Permulation File Override 
  if(usePermulations){
    if(!is.na(cmdArgImport('f'))){
      permulationPValOverride = cmdArgImport('f')
      message("Using Manually specified permulation p-value file. NOTE this setting ignores the meta-combination argument.")
    }else{
      message("No Permulation pValue override, using standard CombinedPrunedFastAllPermulationsPValue.rds.")
    } 
  }
  
  #Gene Enrichment Use 
  if(!is.na(cmdArgImport('g'))){
    useGeneEnrichment = cmdArgImport('g')
    useGeneEnrichment = as.logical(useGeneEnrichment)
    if(is.na(useGeneEnrichment)){
      useGeneEnrichment = TRUE
      message("Include geneset enrichment value not interpretable as logical. Did you remember to capitalize? Using TRUE.")
    }
  }else{
    message("Include geneset enrichment value not specified, using TRUE.")
  }
}

#                   ------- Code Body -------- 


# ------ Import the Data ------ 
correlationFileLocation = paste(outputFolderName, filePrefix, "CorrelationFile.rds", sep= "")
correlData = readRDS(correlationFileLocation)                            #Import the correlation data (non-permulated)

# - Permulations - 
if(usePermulations){
  if(!is.null(permulationPValOverride)){
    permulationFileLocation = permulationPValOverride
  }else{
    permulationFileLocation = paste(outputFolderName, filePrefix, permulationDefaultFilename, sep= "")
  }
  
  correlData$permPValue = readRDS(permulationFileLocation)                       #Add a collumn to the data with the permulation p Values
}

# - Q values - 
qValueNoperm = qvalue(p = correlData$p.adj)
correlData$qValueNoperm = qValueNoperm$qvalues

if(usePermulations){
  qValuePerm = qvalue(p = correlData$permPValue)
  correlData$qValuePerm = qValuePerm$qvalues
}

# - Negative and Positive Rho 
correlDataPositive = correlData[which(correlData$Rho > 0),]
correlDataNegative = correlData[which(correlData$Rho < 0),]
correlDataZero = correlData[which(correlData$Rho == 0),]

# - Enrichment results - 
if(useGeneEnrichment){
  outputFolderFilenames = list.files(paste("./",outputFolderName, sep=""))
  enrichmentFileNames = outputFolderFilenames[grep("Enrichment-", outputFolderFilenames)]
  if(length(enrichmentFileNames) == 0){
    enrichmentRange = NA
    useGeneEnrichment = FALSE
    message("No enrichment files found. Deactivating enrichment analysis processing.")
  }else{
    enrichmentRange = length(enrichmentFileNames)
    enrichmentResultSets = NULL
    for(i in 1:enrichmentRange){
      EnrichmentSetNumber = paste("enrichment", i, sep="")
      EnrichmentData = readRDS(paste(outputFolderName, enrichmentFileNames[i], sep=""))
      enrichmentResultSets[i] = EnrichmentData
      names(enrichmentResultSets)[i] = names(EnrichmentData)
      assign(EnrichmentSetNumber, EnrichmentData)
    }
  }
}


# ------ Make plots ------ 

# - p value histograms - 
makePHistogram = function(data, column, titleVal){
  phist = ggplot(data, aes(x=data[[column]]))+
    geom_histogram(binwidth = 0.02, alpha=0.9, col="white", boundary=0)+
    scale_x_continuous(breaks = seq(0,1,0.1), lim = c(0,1))+
    labs(title = titleVal)+
    xlab("p Values")+
    theme(
          #axis.title.x = element_blank(),
          plot.title = element_text(size=18, hjust = 0.5),
          plot.title.position = "plot"
        )
}

postiveRhoNonpermHistogram = makePHistogram(correlDataPositive, "p.adj", "Positive Rho Non-permulated")
negativeRhoNonpermHistogram = makePHistogram(correlDataNegative, "p.adj", "Positive Rho Non-permulated")

if(usePermulations){
  postiveRhoPermHistogram = makePHistogram(correlDataPositive, "permPValue", "Positive Rho Permulated")
  negativeRhoPermHistogram = makePHistogram(correlDataNegative, "permPValue", "Positive Rho Permulated")
  
}
#I think I want to display positive and negative side-by-side, and permulated non-permulated above/below eachother.


# - List plots - 
makePvListPlot = function(data, column, length, titleVal){
  data = data[order(data[[column]]),]
  dataHead = data[1:length,]
  dataFront = dataHead[,1:2]
  dataBack = dataHead[,3:ncol(dataHead)]
  dataBack = format_table(dataBack, pretty_names = T, digits = "scientific3")
  dataMain = data.frame(rownames(dataHead))
  dataMain = append(dataMain, dataFront)
  names(dataMain)[1] = "Gene"
  dataMain = append(dataMain, dataBack)
  dataMain = as.data.frame(dataMain)
  listPlot = ggplot()+
    theme_void()+
    annotate(geom = "table",
             x=1,
             y=1,
             label = dataMain)+
    labs(title = titleVal)+
    theme(plot.title = element_text(size=18, hjust = 0.5, vjust=0))
}

topGenesPadj = makePvListPlot(correlData, "p.adj", 25, "Top genes non-permulated")
topPositiveGenesPadj = makePvListPlot(correlDataPositive, "p.adj", 10,"Top Positive genes non-permulated")
topNegativeGenesPadj = makePvListPlot(correlDataNegative, "p.adj", 10, "Top Negative genes non-permulated")

topGenesNPQ = makePvListPlot(correlData, "qValueNoperm", 25, "Top genes Q-Value non-permulated")
topPositiveGenesNPQ = makePvListPlot(correlDataPositive,  "qValueNoperm",10, "Top Positive genes Q-Value non-permulated")
topNegativeGenesNPQ = makePvListPlot(correlDataNegative,  "qValueNoperm", 10,"Top Negative genes Q-Value non-permulated")

if(usePermulations){
  topGenesPerm = makePvListPlot(correlData, "permPValue", 25, "Top genes Permulated")
  topPositiveGenesPerm = makePvListPlot(correlDataPositive,  "permPValue", 10,"Top Positive genes Permulated")
  topNegativeGenesPerm = makePvListPlot(correlDataNegative,  "permPValue", 10,"Top Negative genes Permulated")
  
  topGenesPermQ = makePvListPlot(correlData, "qValueNoperm", 25, "Top genes Q-Value Permulated")
  topPositiveGenesPermQ = makePvListPlot(correlDataPositive,  "qValueNoperm",10, "Top Positive genes Q-Value Permulated")
  topNegativeGenesPermQ = makePvListPlot(correlDataNegative,  "qValueNoperm",10, "Top Negative genes Q-Value Permulated")
}

# - Gene Enrichment Plots - 
if(useGeneEnrichment){
  enrichmentPlotSet = list()
  makeGeListPlot = function(data, column, length, titleVal, decreasing){
    setData = data[[1]]
    if(decreasing){
    setData = setData[order(abs(setData[[column]]), decreasing = T),]
    }else{
      setData = setData[order(abs(setData[[column]])),]
    }
    dataHead = setData[1:length,]
    dataHead$gene.vals= strsplit(dataHead$gene.vals, ",")
    for(i in 1:length){
      genesetlist = dataHead$gene.vals[[i]][1:6]
      genesetlistsingle = paste(genesetlist[1], genesetlist[2], genesetlist[3], genesetlist[4], genesetlist[5],genesetlist[6])
      dataHead$gene.vals[i] = genesetlistsingle
    }
    dataFront = dataHead[,c(1,4)]
    dataBack = dataHead[,c(2,3,5)]
    dataBack = format_table(dataBack, pretty_names = T, digits = "scientific3")
    dataMain = data.frame(substring(rownames(dataHead), 1, 40))
    dataMain = append(dataMain, dataFront)
    names(dataMain)[1] = "Geneset"
    dataMain = append(dataMain, dataBack)
    dataMain = as.data.frame(dataMain)
    listPlot = ggplot()+
      theme_void()+
      annotate(geom = "table",
               x=1,
               y=1,
               label = dataMain)+
      labs(title = paste(names(data), titleVal))+
      theme(plot.title = element_text(size=18, hjust = 0.5, vjust=1))
  }
  for(i in 1:enrichmentRange){
    genesetPlotName = paste("genesetPlot", i, sep="")
    if(usePermulations){
      genesetPlot = makeGeListPlot(enrichmentResultSets[i], "stat", 40, "Top pathways by permulation", T)
    }else{
      genesetPlot = makeGeListPlot(enrichmentResultSets[i], "stat", 40, "Top pathways by non-permulation", T)
    }
    enrichmentPlotSet[i] = list(genesetPlot)
    assign(genesetPlotName, genesetPlot)
    
  }
}
# ------- Make Plots ------- 

# - Plot the top genes - 
if(useBoth){
  headlineGenes = plot_grid(topGenesPerm, topGenesPermQ, topGenesPadj, topGenesNPQ, ncol = 2, nrow = 2)
  headlineRows = 2
}else if(usePermulations){
  headlineGenes = plot_grid(topGenesPerm, topGenesPermQ, ncol = 2, nrow = 1)
  headlineRows = 1
}else{
  headlineGenes = plot_grid(topGenesPadj, topGenesNPQ, ncol = 2, nrow = 1)
  headlineRows = 1
}

# - plot the histograms - 
histograms = NULL
if(useBoth){
  histograms = plot_grid(postiveRhoPermHistogram, negativeRhoPermHistogram, postiveRhoNonpermHistogram, negativeRhoNonpermHistogram, ncol = 2, nrow = 2)
  histogramRows = 2
}else if(usePermulations){
  histograms = plot_grid(postiveRhoPermHistogram, negativeRhoPermHistogram, ncol = 2, nrow = 1)
  histogramRows = 1
}else{
  histograms = plot_grid(postiveRhoNonpermHistogram, negativeRhoNonpermHistogram, ncol = 2, nrow = 1)
  histogramRows = 1
}

# - plot the signed genes - 
if(useBoth){
  signedGenes = plot_grid(topPositiveGenesPerm, topNegativeGenesPerm, topPositiveGenesPadj, topNegativeGenesPadj, ncol = 2, nrow = 2)
  signedRows = 2
}else if(usePermulations){
  signedGenes = plot_grid(topPositiveGenesPerm, topNegativeGenesPerm, ncol = 2, nrow = 1)
  signedRows = 1
}else{
  signedGenes = plot_grid(topPositiveGenesPadj, topNegativeGenesPadj, ncol = 2, nrow = 1)
  signedRows = 1
}

# - plot enrichments - 
enrichmentRows = 0 
if(useGeneEnrichment){
  enrichmentPlots = enrichmentPlotSet[1][1]
  for(i in 2:enrichmentRange){
    enrichmentPlots = plot_grid(genesetPlot1, genesetPlot2, genesetPlot3, genesetPlot4, genesetPlot5, ncol = 1, nrow = 3)
  }
  enrichmentRows = length(enrichmentRange)
}

# ------ Output to pdf ------ 
pdfRows = headlineRows+histogramRows+signedRows+enrichmentRows
pdfLengthPerRow = 6.5
pdfLength = pdfRows*pdfLengthPerRow

outputPDFLocation = paste(outputFolderName, filePrefix, "VisualizeOutput.pdf", sep= "") # this could be improved to be more dynamic
pdf(file = outputPDFLocation, width = 15, height = pdfLength)

plot_grid(headlineGenes, histograms, signedGenes, ncol = 1, nrow = 3)

if(useGeneEnrichment){
  enrichmentPlots
}
dev.off()

obj = paste("genesetPlot", i, sep="")
eval(parse(text = obj))
#improvements for this script: 
  #Update output title to be more dynamic
  #fix justification of figure titles
  #Change row length to match length of content (shorten the signed row) 
pdf(file = outputPDFLocation, width = 12, height = 12)
genesetPlot5
dev.off()
