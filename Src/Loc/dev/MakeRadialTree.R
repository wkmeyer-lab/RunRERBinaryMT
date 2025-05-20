library(ggtree)
library(ggimage)
library(rphylopic)
library(RERconverge)
source("Src/Reu/ZoonomTreeNameToCommon.R")
source("Src/Reu/cmdArgImport.R")
{
  nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
  rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
  offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
  offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
  child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
  parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")
}

args =c("r=CategoricalInsVertivoreTree")

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


palette(c( "darkgreen", "darkblue", "black", "red"))
CategoryReplacements = c("Herbivore", "Invertivore", "Omnivore", "Vertivore")
mainTreesLocation = "data/zoonomiaAllMammalsTrees.rds"
spreadSheetLocation = "Data/mergedData.csv"
nameColumn = "ZoonomiaTip"

if(!exists("mainTrees")){mainTrees = readRDS(mainTreesLocation)}
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
for(i in 1:length(unique(edge$CategorylengthChar))){
  edge$CategorylengthChar[edge$CategorylengthChar == i] = CategoryReplacements[i]
}


commonCategoricalTree$edge.length = commonMasterTrimmed$edge.length
scientificCategoricalTree$edge.length = scientificMasterTrimmed$edge.length




phylopicNames = NULL
uuidList = NULL
missingPictures = NA

uuidListFilename =  paste(outputFolderName, filePrefix, "UuidList.rds", sep="") #make a filename based on the prefix
if(!uuidListFilename | forceUpdate){                             #if it does not exist, or update is forced 
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

uuidList[uuidList == "NULL"] = NULL

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
#ggTreeOut = ggTreeOut + geom_tiplab(geom = "phylopic", aes(image = uuid))
#ggTreeOut + geom_phylopic(aes(uuid = uuid), color = "black", alpha = 1, size = 0.08)

pdf()
ggTreeOut
dev.off()

{
  collapsedClades = data.frame()
  collapsedClades[1,] = NA
  
  
  
  collapsedClades$Cats = MRCA(commonCategoricalTree, c("Jaguar", "Lion", "Cheetah"))
}

ggTreeClades = ggtree(commonCategoricalTree, layout = "circular") +scale_color_manual(values=palette())
ggTreeClades = ggTreeClades %<+% edge + aes(color=CategorylengthChar)


for(i in 1:ncol(collapsedClades)){
  #ggTreeClades = ggTreeClades + geom_cladelabel(collapsedClades[1,i], names(collapsedClades)[i])
  ggTreeClades = ggTreeClades + geom_cladelabel(collapsedClades[1,i], NA, color = "gray", barsize = 1)
}
ggTreeClades


?geom_cladelabel
# Collapse some clades
{
  collapsedClades = data.frame()
  collapsedClades[1,] = NA
  
  collapsedClades$Bovidae = MRCA(hiller4MasterCommonTrimmed, c("Zebu", "Asian Water Buffalo"))
  
  collapsedClades$Caprids = MRCA(hiller4MasterCommonTrimmed, c("Bighorn sheep", "Tibetan antelope"))
  
  collapsedClades$Deer = MRCA(hiller4MasterCommonTrimmed, c("Central European Red Deer", "White-tailed Deer"))
  
  collapsedClades$Cetacea = MRCA(hiller4MasterCommonTrimmed, c("Bottle-nose Dolphin", "Sperm whale"))
  
  collapsedClades$Camelidae = MRCA(hiller4MasterCommonTrimmed, c("Ferus Camel", "Alpaca"))
  
  collapsedClades$Pinnipeds = MRCA(hiller4MasterCommonTrimmed, c("Hawaiian Monk Seal", "Pacific walrus"))
  
  collapsedClades$Mustelidae = MRCA(hiller4MasterCommonTrimmed, c("Ferret", "Red Panda "))
  
  collapsedClades$Bears = MRCA(hiller4MasterCommonTrimmed, c("Polar Bear", "Panda"))
  
  collapsedClades$Dogs = MRCA(hiller4MasterCommonTrimmed, c("African Hunting Dog ", "Dog"))
  
  collapsedClades$Cats = MRCA(hiller4MasterCommonTrimmed, c("Tiger", "Cat"))
  
  collapsedClades$Pangolins = MRCA(hiller4MasterCommonTrimmed, c("Chinese pangolin", "Sunda pangolin"))
  
  collapsedClades$Horses = MRCA(hiller4MasterCommonTrimmed, c("Horse", "Wild Donkey"))
  
  collapsedClades$Megabats = MRCA(hiller4MasterCommonTrimmed, c("Black flying-fox", "Rousette Fruit Bats"))
  
  collapsedClades$leaftnosedBats = MRCA(hiller4MasterCommonTrimmed, c("Horseshoe Bats", "Old World Leaf-nosed Bats"))
  
  collapsedClades$bigBrownBat = MRCA(hiller4MasterCommonTrimmed, c("Microbat", "Long-winged bats"))
  
  collapsedClades$Eulipotyphla = MRCA(hiller4MasterCommonTrimmed, c("Hedgehog", "Shrew"))
  
  collapsedClades$Maccaca = MRCA(hiller4MasterCommonTrimmed, c("Rhesus Macaca", "Green monkey"))
  
  collapsedClades$Langur = MRCA(hiller4MasterCommonTrimmed, c("Proboscis Monkey", "Black-and-white Colobus Monkey"))
  
  collapsedClades$Apes = MRCA(hiller4MasterCommonTrimmed, c("Chimp", "Gibbon"))
  
  collapsedClades$Monkeys = MRCA(hiller4MasterCommonTrimmed, c("Squirrel monkey", "Night Monkey"))
  
  collapsedClades$Lemur = MRCA(hiller4MasterCommonTrimmed, c("Sifakas", "Bushbaby"))
  
  collapsedClades$Rat = MRCA(hiller4MasterCommonTrimmed, c("Ryukyu mouse", "Rat"))
  
  collapsedClades$Hamster = MRCA(hiller4MasterCommonTrimmed, c("Chinese hamster", "Deer mouse"))
  
  collapsedClades$Beaver = MRCA(hiller4MasterCommonTrimmed, c("Ord Kangaroo Rat", "Beaver"))
  
  collapsedClades$Chinchilla = MRCA(hiller4MasterCommonTrimmed, c("Brush-tailed rat", "Domestic guinea pig"))
  
  collapsedClades$Molerats = MRCA(hiller4MasterCommonTrimmed, c("Naked mole-rat", "Damara mole rat"))
  
  collapsedClades$Squirrel = MRCA(hiller4MasterCommonTrimmed, c("Ground Squirrel ", "Marmot"))
  
  collapsedClades$Rabbit = MRCA(hiller4MasterCommonTrimmed, c("Rabbit", "Pika"))
  
  collapsedClades$Elephant = MRCA(hiller4MasterCommonTrimmed, c("Manatee", "Hydrax"))
  
  collapsedClades$Aardvark = MRCA(hiller4MasterCommonTrimmed, c("Cape elephant shrew", "Aardvark"))
  
  collapsedClades$Marsupials = MRCA(hiller4MasterCommonTrimmed, c("Tasmanian devil", "Opossum"))
}
#

ggTreeClades = ggtree(hiller4MasterCommonTrimmed, layout = "circular", linewidth = 1) + scale_color_manual(values=c("black", "firebrick", "lightgreen", "steelblue"))

ggTreeClades = ggTreeClades %<+% edge + aes(color=CategorylengthChar)

for(i in 1:ncol(collapsedClades)){
  #  ggTreeClades = ggTreeClades + geom_cladelabel(collapsedClades[1,i], names(collapsedClades)[i])
  ggTreeClades = ggTreeClades + geom_cladelabel(collapsedClades[1,i], NA)
}
ggTreeClades

png("Output/NewHiller4Phen/NewHiller4PhenCicleTree.png", width = 2000, height = 2000)
ggTreeClades
dev.off()

?png
?geom_cladelabel

hiller4MasterCommonTrimmed$tip.label

collapsedClades = marsupials

ggTreeOut2 = ggtree(hiller4MasterCommonTrimmed) + geom_tree()

ggTreeOut 
