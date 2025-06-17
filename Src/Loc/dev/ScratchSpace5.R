a = b #prevent full runs


# ------------------------------------------------------------------
# --- Work on making a tree figure for lalitha  ----- 
# ------------------------------------------------------------------
lalithaData = read.csv("C:/Users/mit221/Downloads/SpeciesNamesAndPhenosForDE72.csv")


laltihaSpecies = lalithaData[,1]

#laltihaSpecies = unlist(laltihaSpecies)
#laltihaSpecies = laltihaSpecies[-1]
#laltihaSpecies = gsub("[0-9]+$", "", laltihaSpecies)
#laltihaSpecies = unique(laltihaSpecies)

mergeData = read.csv("Data/mergedData.csv")

for(i in 1:length(laltihaSpecies)){
laltihaSpecies[i] = gsub(" ", "_", laltihaSpecies[i]) #replace spaces with underscores 
laltihaSpecies[i] = sub('^([^_]+_[^_]+).*', '\\1', laltihaSpecies[i]) #remove anything  after a second underscore
laltihaSpecies[i] = tolower(laltihaSpecies[i])
}

lalithaData[,1] = laltihaSpecies
write.csv(lalithaData, "C:/Users/mit221/Downloads/SpeciesNamesAndPhenosForDE72.csv")

sum(!laltihaSpecies %in% mergeData$Scientific_Binomial)


ggTreeOut = ggtree(commonCategoricalTree) +scale_color_manual(values=palette()) 


# ------------------------------------------------------------------
# --- Making a make basic radial tree script  ----- 
# ------------------------------------------------------------------


commonMainTrees = mainTrees
commonMainTrees$masterTree = ZoonomTreeNameToCommon(commonMainTrees$masterTree, manualAnnotLocation = spreadSheetLocation, tipCol = nameColumn)
commonMasterTree = commonMainTrees$masterTree


categoricalCommonTreeFilename = paste(outputFolderName, filePrefix, "CategoricalCommonTree.rds", sep="") #make a filename based on the prefix
commonCategoricalTree = readRDS(categoricalCommonTreeFilename)
scientificCategoricalTreeFilename = paste(outputFolderName, filePrefix, "CategoricalScientificTree.rds", sep="") #make a filename based on the prefix
scientificCategoricalTree = readRDS(scientificCategoricalTreeFilename)


commonMasterTrimmed = drop.tip(commonMasterTree, commonMasterTree$tip.label[!commonMasterTree$tip.label %in% commonCategoricalTree$tip.label])
scientificMasterTrimmed = ZoonomTreeNameToCommon(commonMasterTrimmed, manualAnnotLocation = spreadSheetLocation, tipCol = "CommonName", scientific = T, scientificCol = "Scientific_Binomial")

commonCategoricalTreeEdgeLengths = commonCategoricalTree$edge.length
commonCategoricalTreeEdgeLengths = as.character(commonCategoricalTreeEdgeLengths)
edge=data.frame(commonCategoricalTree$edge, edge_num=1:nrow(commonCategoricalTree$edge))
colnames(edge)=c("parent", "node", "edge_num")
edge$Categorylength = commonCategoricalTree$edge.length
edge$CategorylengthChar = as.character(edge$Categorylength)
if(!is.null(CategoryReplacements)){
  for(i in 1:length(unique(edge$CategorylengthChar))){
    edge$CategorylengthChar[edge$CategorylengthChar == i] = CategoryReplacements[i]
  }
}

commonCategoricalTree$edge.length = commonMasterTrimmed$edge.length
scientificCategoricalTree$edge.length = scientificMasterTrimmed$edge.length




phylopicNames = NULL
uuidList = NULL
missingPictures = NA
uuidListFilename =  paste(outputFolderName, filePrefix, "UuidList.rds", sep="") #make a filename based on the prefix
if(!file.exists(uuidListFilename) | forceUpdate){                             #if it does not exist, or update is forced 
  inTips = scientificCategoricalTree$tip.label
  for(i in 1:length(inTips)){
    print(i)
    phylopicNames[i] = tryCatch({autocomplete_name(inTips[i])[1,2]}, error = function(msg) {return(inTips[i])})
    
    uuidList[i] = tryCatch(
      {get_uuid(phylopicNames[i])}, 
      error = function(msg){
        if(length(grep(" ", phylopicNames[i])) > 0){genusName = strsplit(phylopicNames[i], " ")[[1]][1]}else{
          if(length(grep("_", phylopicNames[i])) > 0){genusName = strsplit(phylopicNames[i], "_")[[1]][1]}
        }
        tryCatch(
          {get_uuid(genusName)},
          error = function(msg){
            message(paste("No pic found for ", phylopicNames[i]))
            missingReport = phylopicNames[i]
            names(missingReport) = i 
            append(missingPictures, missingReport)
            return("NULL")
          }
        )
      }
    )
  }
  saveRDS(uuidList, uuidListFilename)
}else{                                                                          #Otherwise
  uuidList = readRDS(uuidListFilename)                                              #Use the existing ones
}


# Replace any missing UUIDs with tardigrades as debug images 
uuidList[uuidList == "NULL"] = get_uuid("tardigrades")
uuidList[21] = get_uuid("tardigrades")

#uuidList[uuidList == "NULL"] = NULL

tip_data = data.frame(
  scientificlabel = scientificCategoricalTree$tip.label,
  uuid = uuidList,
  node = 1:length(scientificCategoricalTree$tip.label),
  stringsAsFactors = FALSE
)
tip_data$uuid[tip_data$uuid == "NULL"] = NULL


ggTreeOut = ggtree(commonCategoricalTree, layout = "circular") +scale_color_manual(values=palette()) 
ggTreeOut = ggTreeOut %<+% edge + aes(color=CategorylengthChar)
ggTreeOut = ggTreeOut %<+% tip_data 
ggTreeOut$data$label = paste(ggTreeOut$data$label, "-", ggTreeOut$data$node, sep="")
ggTreeOut = ggTreeOut + geom_tiplab()
ggTreeOut = ggTreeOut + geom_tiplab(geom = "phylopic", aes(image = uuid))
#ggTreeOut + geom_phylopic(aes(uuid = uuid), color = "black", alpha = 1, size = 0.08)
ggTreeOut


# ------------------------------------------------------------------
# --- Assess driving diet of unqieu carnivory results  ----- 
# ------------------------------------------------------------------

combinedData = readRDS(paste0(combinedDataFilename, ".rds"))

carnivoryUNique = combinedData[which(combinedData$`CH-significant` & !combinedData$`HI-significant` & !combinedData$`HV-significant`),]


table(carnivoryUNique$`CH-Driver`)
table(sign(carnivoryUNique$`CH-Rho`))
# Right. I don't current have driver information because I haven't run the driver analysis for the 
# I need to see about making that. 

# quick fix for the driver analysis to retarget because this is not technically a part of the main analysis 

# ------------------------------------------------------------------
# --- Get set of genes most different between HI and HV  ----- 
# ------------------------------------------------------------------

combinedResults # from MakeOverlapFigure
combinedResultsModified = combinedResults


getComparisionDifference = function(dataframe, colOne, colTwo){
  colOneIndex = names(dataframe)[which(names(dataframe) == colOne)]
  colTwoIndex = names(dataframe)[which(names(dataframe) == colTwo)]
  
  distanceFromEqual = abs(dataframe[colOneIndex] - dataframe[colTwoIndex]) / sqrt(2)
  distanceFromEqual
}



getComparisionDifference(combinedResultsModified, "HI-Rho", "HV-Rho")
combinedResults$`HI-HV-Delta` = getComparisionDifference(combinedResultsModified, "HI-Rho", "HV-Rho")


test = combinedResultsModified[order(combinedResultsModified$`HI-HV-Delta`, decreasing = T),]


df$distance_from_y_eq_x <- abs(df$x - df$y) / sqrt(2)




# ------------------------------------------------------------------
# --- examining the GO sets in specific overlap sections ----- 
# ------------------------------------------------------------------
GoSignificanceResults # from MakeOverlapFigure

View(GoSignificanceResults)

# get pathways that show up in compionents but not carnivory 
IVnCpathways = rownames(GoSignificanceResults[which(GoSignificanceResults$`HI-significant` & GoSignificanceResults$`HV-significant` & !GoSignificanceResults$`CH-significant`),])

which(rownames(GoCombinedResults) %in% IVnCpathways)
IVnCResults = GoCombinedResults[which(rownames(GoCombinedResults) %in% IVnCpathways), ]
View(IVnCResults)
cat(IVnCpathways)

# get the pathways which are unique to vertivory 
VonlyPathways = rownames(GoSignificanceResults[which(!GoSignificanceResults$`HI-significant` & GoSignificanceResults$`HV-significant` & !GoSignificanceResults$`CH-significant`),])
IonlyPathways = rownames(GoSignificanceResults[which(GoSignificanceResults$`HI-significant` & !GoSignificanceResults$`HV-significant` & !GoSignificanceResults$`CH-significant`),])


# ------------------------------------------------------------------
# --- getting permualtion RERResult data fro rho plot creation ----- 
# ------------------------------------------------------------------

permulationIntermediate = readRDS(paste0(outputFolderName, "CategoricalInsVertivoreTreeCategoricalPermulationsIntermediates101.rds"))


permulationHIStatsValues = permulationIntermediate$Peffsize$`1 - 3`
permulationHVStatsValues = permulationIntermediate$Peffsize$`1 - 4`

permulationHIStatsCol = permulationHIStatsValues[,1]
permulationHVStatsCol = permulationHVStatsValues[,1]


rhoValuesPerm = data.frame(permulationHIStatsCol, permulationHVStatsCol)

for(i in 1:length(rhoValuesPerm)){
  xName = names(rhoValuesPerm)[i]
  if(i+1 <= length(rhoValuesPerm)){
    for(j in (i+1):length(rhoValuesPerm)){
      yName = names(rhoValuesPerm)[j]
      yLabel =  paste0(replacePrefixWithName(addDashes(gsub("-Rho", "", yName))), " Stat")
      xLabel =  paste0(replacePrefixWithName(addDashes(gsub("-Rho", "", xName))), " Stat")
      
      rhoCorrellPlot = ggplot(rhoValuesPerm, aes(x = .data[[xName]], y = .data[[yName]])) + 
        geom_point() + geom_pointdensity() + scale_color_viridis(name = "Number of nearby genes", limits = densityScale) + 
        stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),formula = y ~ x,parse = TRUE) +
        theme_classic()+
        xlab(xLabel) + ylab(yLabel)
    
      
      netIndex = netIndex +1
      rhoPlotSet[[netIndex]] = rhoCorrellPlot
      names(rhoPlotSet)[netIndex] = paste(xName, yName, sep="-")
    }
  }
}

rhoDfList = list()
permulationRhoDataframes = for(i in 1:2500){
  rhoDf = data.frame(permulationHIStatsValues[,i], permulationHVStatsValues[,i])
  rhoDfList[[i]] = rhoDf
}

rSquaredSet = NULL
for(k in 1:length(rhoDfList)){
  rhoDf = rhoDfList[[k]]
  names(rhoDf) = c("HI", "HV")
  model <- lm(HV ~ HI, data = rhoDf)
  rSquared = summary(model)$r.squared
  rSquaredSet = append(rSquaredSet, rSquared)
  #cat(rSquared, "\n")
}
pdf()
hist(rSquaredSet)
dev.off()

mean(rSquaredSet)
