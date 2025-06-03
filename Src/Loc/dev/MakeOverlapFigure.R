library(RERconverge)
library(ggvenn)
source("Src/Reu/cmdArgImport.R")
clusterRun = F
clusterRun = T


# -- Making new venn diagrams -- 
significanceCutoff = 0.05
prefix = "CategoricalInsvertivoreTree"
pairwiseSets = c("Herbivore-Insectivore", "Herbivore-Vertivore")

args = c("r=CategoricalInsvertivoreTree")

# -- Standard Startup code -- 
if(clusterRun)args = commandArgs(trailingOnly = TRUE)
{  # Bracket used for collapsing purposes
  #File Prefix
  if(!is.na(cmdArgImport('r'))){
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
    message("Force update not specified, not forcing update")
  }
}



# -- Make central RERData object 
pairwiseCorrelationFileName = paste(outputFolderName, filePrefix, "PairwiseCorrelationFile.rds", sep= "") #make a name for the pairwise comparisons based on prefix
correlationResults = readRDS(pairwiseCorrelationFileName)

combinedResults = NA
combinedDrivers = NA

for(i in 1:length(pairwiseSets)){
  currentSet = pairwiseSets[i]
  correlationSubsetName = gsub("-", " - ", currentSet)
  correlationPrefix = paste(substr(strsplit(currentSet, split = "-")[[1]],1,1), collapse = '')
  which(names(correlationResults) == correlationSubsetName) 
  
  currentResults = correlationResults[[i]]
  names(currentResults) = paste0(correlationPrefix, "-", names(currentResults))
  
  combinedResults = cbind(combinedResults, currentResults)
  
  driverFilename = paste0(outputFolderName, currentSet, "/", filePrefix, currentSet, "DirectionalityTable.rds")
  if(file.exists(driverFilename)){
    driverTable = readRDS(driverFilename)
    driverTable = driverTable[,-grep("main", names(driverTable))]
    names(driverTable)[which(names(driverTable) == "directionality")] = paste0(correlationPrefix, "-", "directionality")
    names(driverTable)[which(names(driverTable) == "directionNumeric")] = paste0(correlationPrefix, "-", "directionalityNumeric")    
    
    combinedDrivers = cbind(combinedDrivers, driverTable)
  }
  
}
combinedResults = combinedResults[,-1]
combinedDrivers = combinedDrivers[,-1]

combinedResults = cbind(combinedResults, combinedDrivers)




HICorrelations = RERResults$`Herbivore - Insectivore`
HICorrelations$gene = rownames(HICorrelations)
HVCorrelations = RERResults$`Herbivore - Vertivore`
HOCorrelations = RERResults$`Herbivore - Omnivore`

HIDrivingData = readRDS("Output/CategoricalInsvertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreDirectionalityTable.rds")
HIDrivingData$gene = rownames(HIDrivingData)
HVDrivingData = readRDS("Output/CategoricalInsvertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreDriverTable.rds")

CombinedRERData = merge(HICorrelations, HIDrivingData) 
?merge


HISignificantGenes = rownames(HICorrelations)[which(HICorrelations$p.adj < significanceCutoff)]
HVSignificantGenes = rownames(HVCorrelations)[which(HVCorrelations$p.adj < significanceCutoff)]
HOSignificantGenes = rownames(HOCorrelations)[which(HOCorrelations$p.adj < significanceCutoff)]
















# ---- Old code ------ 
RERResults = readRDS("Output/CategoricalInsvertivoreTree/CategoricalInsVertivoreTreePairwiseCorrelationFile.rds")

HICorrelations = RERResults$`Herbivore - Insectivore`
HVCorrelations = RERResults$`Herbivore - Vertivore`
HOCorrelations = RERResults$`Herbivore - Omnivore`


HISignificantGenes = rownames(HICorrelations)[which(HICorrelations$p.adj < 0.05)]
HVSignificantGenes = rownames(HVCorrelations)[which(HVCorrelations$p.adj < 0.05)]
HOSignificantGenes = rownames(HOCorrelations)[which(HOCorrelations$p.adj < 0.05)]

sharedGenes = HISignificantGenes[which(HISignificantGenes %in% HVSignificantGenes)]

triSharedGenes = HOSignificantGenes[which(HOSignificantGenes %in% sharedGenes)]


HIDrivingData = readRDS("Output/CategoricalInsvertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreDirectionalityTable.rds")
HIDrivingData[which(rownames(HIDrivingData) %in% sharedGenes),]


HVDrivingData = read.csv("Output/CategoricalInsvertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreDirectionalityTable.csv")
rownames(HVDrivingData) = HVDrivingData$X
HVDrivingData[which(rownames(HVDrivingData) %in% sharedGenes),]

CombinedDrivingData = HIDrivingData
colnames(CombinedDrivingData)[10] = "HIDirectionality"
colnames(CombinedDrivingData)[11] = "HIDirectionalityNumeric"
CombinedDrivingData = cbind(CombinedDrivingData, HVDrivingData[,c(11,12)])
colnames(CombinedDrivingData)[12] = "HVDirectionality"
colnames(CombinedDrivingData)[13] = "HVDirectionalityNumeric"

all.equal(rownames(HVDrivingData), rownames(HIDrivingData))


SharedCombinedDrivingData = CombinedDrivingData[which(rownames(CombinedDrivingData) %in% sharedGenes),]

nrow(SharedCombinedDrivingData[which(SharedCombinedDrivingData$HIDirectionality == "Herbivore" & SharedCombinedDrivingData$HVDirectionality =="Herbivore"),])
nrow(SharedCombinedDrivingData[which(SharedCombinedDrivingData$HIDirectionality == "Herbivore" | SharedCombinedDrivingData$HVDirectionality =="Herbivore"),])


HISharedDriving = HIDrivingData[which(rownames(HIDrivingData) %in% sharedGenes),]
HVSharedDriving = HVDrivingData[which(rownames(HVDrivingData) %in% sharedGenes),]

rownames(HISharedDriving[which(HISharedDriving$directionality == "Herbivore"),])
rownames(HVSharedDriving[which(HVSharedDriving$directionality == "Herbivore"),])



length(which(rownames(HISharedDriving[which(HISharedDriving$directionality == "Herbivore"),]) %in% rownames(HVSharedDriving[which(HVSharedDriving$directionality == "Herbivore"),])))




# ----------- Group of code 2 

HiGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreCorrelationFile.rds")
HvGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreCorrelationFile.rds")
HoGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Omnivore/CategoricalInsVertivoreTreeHerbivore-OmnivoreCorrelationFile.rds")
IoGeneData = readRDS("Output/CategoricalInsVertivoreTree/Insectivore-Omnivore/CategoricalInsVertivoreTreeInsectivore-OmnivoreCorrelationFile.rds")
IvGeneData = readRDS("Output/CategoricalInsVertivoreTree/Insectivore-Vertivore/CategoricalInsVertivoreTreeInsectivore-VertivoreCorrelationFile.rds")
OvGeneData = readRDS("Output/CategoricalInsVertivoreTree/Omnivore-Vertivore/CategoricalInsVertivoreTreeOmnivore-VertivoreCorrelationFile.rds")
DhiGeneData = readRDS("Output/CategoricalDownsampledInsvertTree/Herbivore-Insectivore/CategoricalDownsampledInsvertTreeHerbivore-InsectivoreCorrelationFile.rds")

arrangeData = function(data, prefix){
  data$index = 1:nrow(data)
  data= data[order(data$p.adj),]
  data$rank = 1:nrow(data)
  data= data[order(data$index),]
  data$index=NULL
  colnames(data) = paste0(prefix, colnames(data))
  return(data)
}

HiGeneData = arrangeData(HiGeneData, "Hi_")
HvGeneData = arrangeData(HvGeneData, "Hv_")
HoGeneData = arrangeData(HoGeneData, "Ho_")
IoGeneData = arrangeData(IoGeneData, "Io_")
IvGeneData = arrangeData(IvGeneData, "Iv_")
OvGeneData = arrangeData(OvGeneData, "Ov_")
DhiGeneData = arrangeData(DhiGeneData, "Dhi_")


combinedData = cbind(HiGeneData, HvGeneData, HoGeneData, IoGeneData, IvGeneData, OvGeneData, DhiGeneData)
combinedData$index = 1:nrow(combinedData)


combinedDataPval = combinedData[,c(2,6,10,14,18,22,26)]
combinedDataPadjVal = combinedData[,c(3,7,11,15,19,23,27)]



length(which(combinedDataPval$Dhi_p.adj < 0.05))

apply(combinedDataPval, MARGIN = 2, mean)

apply(combinedDataPval, 2, function(column) length(which(column < 0.02)))

colnames(combinedDataPval) = c("Herbivore-Insectivore", "Herbivore-Vertivore", "Herbivore-Omnivore", "Insectivore-Omnivore", "Insectivore-Vertivore", "Omnivore-Vertivore", "Downsampled Insectivore-Herbivore")

colnames(combinedDataPval) = c("H-I", "H-V", "H-O", "I-O", "I-V", "O-V", "Downsampled H-I")
combinedDataPval$`Downsampled H-I` = NULL

colnames(combinedDataPadjVal) = c("H-I", "H-V", "H-O", "I-O", "I-V", "O-V", "Downsampled H-I")
combinedDataPadjVal$`Downsampled H-I` = NULL
sigGenesData =data.frame(category = colnames(combinedDataPadjVal), value = apply(combinedDataPadjVal, 2, function(column) length(which(column < 0.02)))) 


library(ggplot2)
library(ggpattern)

ggplot(sigGenesData, aes(x = category, y = value, fill = category, pattern)) +
  geom_bar(stat = "identity", color = "black", show.legend = FALSE) +  # Base bar color
  geom_bar_pattern(
    stat = "identity",
    pattern = "stripe",  # Options: "stripe", "crosshatch", "dots", etc.
    pattern_density = 0.25,
    pattern_fill = c("darkblue", "red", "black","black", "red", "red"),  # Pattern color
    aes(pattern = Category),  # Apply pattern per category
    show.legend = FALSE
  ) +
  theme_minimal() +
  labs(title = "Significant Genes per Pairwise Analysis",
       x = "Category", y = "Number of significant genes") +
  scale_fill_manual(values = c("darkgreen", "darkgreen", "darkgreen", "darkblue", "darkblue", "black"))+
  scale_pattern_fill_manual(values = c("darkblue", "black", "red", "black", "red", "red"))




?barplot
geom_ba
CarnivoreGeneData = readRDS("Output/CategoricalPrunedCarnivoreTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreCorrelationFile.rds")


InsectivoreGeneData$index = 1:nrow(InsectivoreGeneData)
VertivoreGeneData$index = 1:nrow(VertivoreGeneData)
CarnivoreGeneData$index = 1:nrow(CarnivoreGeneData)

InsectivoreGeneData = InsectivoreGeneData[order(InsectivoreGeneData$p.adj),]
VertivoreGeneData = VertivoreGeneData[order(VertivoreGeneData$p.adj),]
CarnivoreGeneData = CarnivoreGeneData[order(CarnivoreGeneData$p.adj),]

InsectivoreGeneData$rank = 1:nrow(InsectivoreGeneData)
VertivoreGeneData$rank = 1:nrow(VertivoreGeneData)
CarnivoreGeneData$rank = 1:nrow(CarnivoreGeneData)

InsectivoreGeneData = InsectivoreGeneData[order(InsectivoreGeneData$index),]
VertivoreGeneData = VertivoreGeneData[order(VertivoreGeneData$index),]
CarnivoreGeneData = CarnivoreGeneData[order(CarnivoreGeneData$index),]

colnames(InsectivoreGeneData) = paste0("I_", colnames(InsectivoreGeneData))
colnames(VertivoreGeneData) = paste0("V_", colnames(VertivoreGeneData))
colnames(CarnivoreGeneData) = paste0("C_", colnames(CarnivoreGeneData))

combinedData = cbind(InsectivoreGeneData, VertivoreGeneData, CarnivoreGeneData)


# ---- part 3


InsectivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-GO_Biological_Process_2023.rds")
write.csv(InsectivoreGoData, "Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-GO_Biological_Process_2023.csv")

?correlateWithBinaryPhenotype
install.packages()

#----------------------------------------

InsectivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-GO_Biological_Process_2023.rds")
VertivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreEnrichment-GO_Biological_Process_2023.rds")
CarnivoreGoData = readRDS("Output/CategoricalPrunedCarnivoryTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreEnrichment-GO_Biological_Process_2023.rds")


#----------------------------------------
library(ggvenn)
InsectivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-GO_Biological_Process_2023.rds")[[1]]
VertivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreEnrichment-GO_Biological_Process_2023.rds")[[1]]
CarnivoreGoData = readRDS("Output/CategoricalPrunedCarnivoreTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreEnrichment-GO_Biological_Process_2023.rds")[[1]]

InsectivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreEnrichment-KeggReactome.rds")[[1]]
VertivoreGoData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreEnrichment-KeggReactome.rds")[[1]]
CarnivoreGoData = readRDS("Output/CategoricalPrunedCarnivoreTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreEnrichment-KeggReactome.rds")[[1]]



InsectivoreGoData = InsectivoreGoData[order(InsectivoreGoData$p.adj),]
VertivoreGoData = VertivoreGoData[order(VertivoreGoData$p.adj),]
CarnivoreGoData = CarnivoreGoData[order(CarnivoreGoData$p.adj),]


signficiantInsectivore = InsectivoreGoData[which(InsectivoreGoData$p.adj <0.05),]
signficiantVertivore = VertivoreGoData[which(VertivoreGoData$p.adj <0.05),]
signficiantCarnivore = CarnivoreGoData[which(CarnivoreGoData$p.adj <0.05),]


valuedInsectivore = InsectivoreGoData[which(InsectivoreGoData$pval <1),]
valuedVertivore = VertivoreGoData[which(VertivoreGoData$pval <1),]
valuedCarnivore = CarnivoreGoData[which(CarnivoreGoData$pval <1),]

vennData = list(
  Insectivore = rownames(valuedInsectivore),
  Vertivore = rownames(valuedVertivore),
  Carnivore = rownames(valuedCarnivore)
)

signficianterInsectivore = InsectivoreGoData[which(InsectivoreGoData$p.adj <0.05),]
signficianterVertivore = VertivoreGoData[which(VertivoreGoData$p.adj <0.05),]
signficianterCarnivore = CarnivoreGoData[which(CarnivoreGoData$p.adj <0.05),]

topInsectivore = InsectivoreGoData[1:100,]
topVertivore = VertivoreGoData[1:100,]
topCarnivore = CarnivoreGoData[1:100,]

vennData = list(
  Insectivore = rownames(signficiantInsectivore),
  Vertivore = rownames(signficiantVertivore),
  Carnivore = rownames(signficiantCarnivore)
)
ggvenn(vennData, fill_color = c("blue", "red", "orange"))


# --------

InsectivoreGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Insectivore/CategoricalInsVertivoreTreeHerbivore-InsectivoreCorrelationFile.rds")
VertivoreGeneData = readRDS("Output/CategoricalInsVertivoreTree/Herbivore-Vertivore/CategoricalInsVertivoreTreeHerbivore-VertivoreCorrelationFile.rds")
CarnivoreGeneData = readRDS("Output/CategoricalPrunedCarnivoreTree/Carnivore-Herbivore/CategoricalPrunedCarnivoreTreeCarnivore-HerbivoreCorrelationFile.rds")

InsectivoreGeneData$index = 1:nrow(InsectivoreGeneData)
VertivoreGeneData$index = 1:nrow(VertivoreGeneData)
CarnivoreGeneData$index = 1:nrow(CarnivoreGeneData)

InsectivoreGeneData = InsectivoreGeneData[order(InsectivoreGeneData$p.adj),]
VertivoreGeneData = VertivoreGeneData[order(VertivoreGeneData$p.adj),]
CarnivoreGeneData = CarnivoreGeneData[order(CarnivoreGeneData$p.adj),]

InsectivoreGeneData$rank = 1:nrow(InsectivoreGeneData)
VertivoreGeneData$rank = 1:nrow(VertivoreGeneData)
CarnivoreGeneData$rank = 1:nrow(CarnivoreGeneData)

InsectivoreGeneData = InsectivoreGeneData[order(InsectivoreGeneData$index),]
VertivoreGeneData = VertivoreGeneData[order(VertivoreGeneData$index),]
CarnivoreGeneData = CarnivoreGeneData[order(CarnivoreGeneData$index),]

colnames(InsectivoreGeneData) = paste0("I_", colnames(InsectivoreGeneData))
colnames(VertivoreGeneData) = paste0("V_", colnames(VertivoreGeneData))
colnames(CarnivoreGeneData) = paste0("C_", colnames(CarnivoreGeneData))

combinedData = cbind(InsectivoreGeneData, VertivoreGeneData, CarnivoreGeneData)
combinedData$V_index = NULL
combinedData$C_index = NULL

combinedData = combinedData[order(combinedData$I_p.adj),]


equationLinePlot = function(data, xIn, yIn){
  linearModel = lm(yIn ~ xIn, data = data) 
  
  equation = paste0("y=", round(coef(linearModel)[2], 2), "*x", round(coef(linearModel)[1], 2))
  rSquared = paste("R² =", round(summary(linearModel)$r.squared, 2))
  
  ggplot(data, aes(x = xIn, y = yIn)) + 
    geom_point()+
    geom_smooth(method = "lm")  
}

equationLinePlot(combinedData, "I_Rho", "C_Rho")

linearModel = lm(I_Rho ~ C_Rho, data = combinedData) 

library(gridExtra)



linearModel = lm(-C_Rho ~ I_Rho, data = combinedData) 
equation = paste0("y=", round(coef(linearModel)[2], 2), "*x", round(coef(linearModel)[1], 2))
rSquared = paste("R² =", round(summary(linearModel)$r.squared, 2))
plot1 = ggplot(combinedData, aes(x = I_Rho, y = -C_Rho)) + 
  geom_point()+
  geom_smooth(method = "lm")+
  annotate("text", x = 3, y = 9, label = paste(equation, rSquared, sep = "\n"), color = "blue", size = 10)+
  theme_minimal()

linearModel = lm(-C_Rho ~ V_Rho, data = combinedData) 
equation = paste0("y=", round(coef(linearModel)[2], 2), "*x", round(coef(linearModel)[1], 2))
rSquared = paste("R² =", round(summary(linearModel)$r.squared, 2))
plot2 = ggplot(combinedData, aes(x = V_Rho, y = -C_Rho)) + 
  geom_point()+
  geom_smooth(method = "lm")+
  annotate("text", x = 3, y = 9, label = paste(equation, rSquared, sep = "\n"), color = "blue", size = 10)+
  theme_minimal()


grid.arrange(plot1, plot2, ncol =2)



?lm

library(ggplot2)
ggplot(combinedData, aes(x = I_Rho, y = -C_Rho)) + 
  geom_point()+
  geom_smooth(method = "lm")


plot(combinedData$I_Rho, combinedData$V_Rho)
plot(-combinedData$C_Rho, combinedData$I_Rho)
plot(-combinedData$C_Rho, combinedData$V_Rho)




plot(combinedData$I_rank, combinedData$V_rank, xlim = c(0,3500), ylim = c(0,3500))
plot(combinedData$`I_rank`, combinedData$`C_rank`, xlim = c(0,4000), ylim = c(0,4000))
plot(combinedData$`V_rank`, combinedData$`C_rank`, xlim = c(0,4000), ylim = c(0,4000))
plot(combinedData$V_index, combinedData$V_p.adj)
hist(combinedData$V_p.adj)
length(combinedData$V_p.adj[which(combinedData$V_p.adj < 1)])
length(combinedData$I_p.adj[which(combinedData$I_p.adj < 1)])

?cbind
sigInsectGenes = InsectivoreGeneData[which(InsectivoreGeneData$p.adj <0.05),]
sigVertGenes = VertivoreGeneData[which(VertivoreGeneData$p.adj <0.05),]
sigCarnGenes = CarnivoreGeneData[which(CarnivoreGeneData$p.adj <0.05),]

geneVennData = list(
  Insectivore = rownames(sigInsectGenes),
  Vertivore = rownames(sigVertGenes),
  Carnivore = rownames(sigCarnGenes)
)
ggvenn(geneVennData, fill_color = c("blue", "red", "orange"))
require(venneuler)
v <- venneuler(c(Insectivore=194, Vertivore=145, Carnviore=598, "Insectivore&Vertivore"=0, "Insectivore&Carnviore"=294, "Vertivore&Carnviore"=125, "Carnviore&Vertivore&Insectivore"=75))
plot(v)

valuedInsectGenes = combinedData[which(combinedData$I_p.adj <1),]
valuedVertGenes = combinedData[which(combinedData$V_p.adj <1),]
valuedCarnGenes = combinedData[which(combinedData$C_p.adj <1),]

geneVennData = list(
  Insectivore = rownames(valuedInsectGenes),
  Vertivore = rownames(valuedVertGenes),
  Carnivore = rownames(valuedCarnGenes)
)

ggvenn(geneVennData, fill_color = c("blue", "red", "pink"))