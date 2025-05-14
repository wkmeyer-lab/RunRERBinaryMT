#------ Run driver assessment code if not already run  --- 
mainPrefix = "CategoricalInsVertivoreTree"
subdirectory = "Herbivore-Vertivore"
BinaryTreeOne = "CategoricalBinaryHerbivoreTree"
BinaryPhenotype = "Herbivore"
BinaryTreeTwo = "CategoricalBinaryInsectivoreTree"
BinaryPhenotypeTwo = "Insectivore"
GOGroup = "KeggReactome"


driverTableFilename = paste0("Output/", mainPrefix, "/", subdirectory, "/", mainPrefix, subdirectory, "DriverTable")
if(!file.exists(paste0(driverTableFilename, ".rds"))){
  source("Src/Reu/AssessRERDirection.R")
  driverTable = AssessRERDirection(mainPrefix,subdirectory , BinaryTreeOne, BinaryPhenotype,  BinaryTreeTwo, BinaryPhenotypeTwo)
  write.csv(driverTable, paste0(driverTableFilename, ".csv"))
  saveRDS(driverTable, paste0(driverTableFilename, ".rds"))
}else{
  driverTable = readRDS(paste0(driverTableFilename, ".rds"))
}

goDriverTableFilename = paste0("Output/", mainPrefix, "/", subdirectory, "/", mainPrefix, subdirectory, "GoDriverTable")
if(!file.exists(paste0(goDriverTableFilename, ".csv"))){
  GODriver = AssessGoCategoryDirection(mainPrefix, subdirectory, GOGroup, 0.1, 1, F)
  View(GODriver)
  write.csv(GODriver, paste0(goDriverTableFilename, ".csv"))
}

#--- split GO Driver into positive and negative

goDriverTableFilename = paste0("Output/", mainPrefix, "/", subdirectory, "/", mainPrefix, subdirectory, "GoDriverTable")
GODriver = read.csv(paste0(goDriverTableFilename, ".csv"))


GoDriverPositive = GODriver[which(GODriver$stat > 0),]
GoDriverNegative = GODriver[which(GODriver$stat < 0),]

GODriverPositiveColored = GODriver
GODriverPositiveColored$p.adj[GODriverPositiveColored$stat < 0] = 0.11
GODriverPositiveColored$pval[GODriverPositiveColored$stat < 0] = 0.11


cytoscapeDirectory = paste0("Output/", mainPrefix, "/", subdirectory, "/", "Cytoscape")
if(!dir.exists(cytoscapeDirectory)){                       #create that directory if it does not exist
  dir.create(cytoscapeDirectory)
}


goPositiveGoDriverTableFilename = paste0("Output/", mainPrefix, "/", subdirectory, "/", "Cytoscape/", mainPrefix, subdirectory, "GoDriverPositiveTable")
goNegativeGoDriverTableFilename = paste0("Output/", mainPrefix, "/", subdirectory, "/", "Cytoscape/", mainPrefix, subdirectory, "GoDriverNegativeTable")
goColoredGoDriverTableFilename = paste0("Output/", mainPrefix, "/", subdirectory, "/", "Cytoscape/", mainPrefix, subdirectory, "GoDriverPositiveColoredTable")

write.table(GODriverPositiveColored, paste0(goColoredGoDriverTableFilename, ".txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
write.csv(GoDriverPositive, paste0(goPositiveGoDriverTableFilename, ".csv"), row.names = F)
write.csv(GoDriverNegative, paste0(goNegativeGoDriverTableFilename, ".csv"), row.names = F)
write.csv(GODriverPositiveColored, paste0(goColoredGoDriverTableFilename, ".csv"), row.names = F)



# ---- Convert a Driver table to cytoscape format 

GODriver
gmtFilename = paste0("Data/", GOGroup, ".gmt")
file.copy(gmtFilename, paste0(cytoscapeDirectory, "/", GOGroup, ".gmt")) #make a copy in the cytoscape direcotry for easy cytoscape work 
GOCytoscape = GODriver[,c(1,7,3,4,2,6)]
colnames(GOCytoscape) = c("GO.ID", "Description", "pVal", "p.adj", "Phenotype", "Gene.vals")
GOCytoscape$Phenotype = sign(GOCytoscape$Phenotype)
#GOCytoscape$Phenotype[GOCytoscape$Phenotype == 1] = "1"
write.table(GOCytoscape, paste0(cytoscapeDirectory, "/CytoscapeInput.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  #Use that cytoscape input file as the data in Cytoscape, on the generic setting


library(GSEABase)
gmtData = getGmt(gmtFilename)

gmtNames = names(gmtData)
gmtDriver = GODriver$Directionality[match(gmtNames, GODriver$X)]
gmtUpdate = data.frame(gmtNames, gmtDriver)
write.csv(gmtUpdate, paste0(cytoscapeDirectory, "/GmtDirectionColumn.csv"), row.names = F)
  #This creates a column with the driving phenotype information. 
  #Open the gmt file in excel, and replace the description column with the produced driver column. 
  #This is set up this way because R is bad at handling variable row lengths, so it is easier to edit the gmt file in excel.

# --- Guide on how to convert these outputs into a cytoscape figure: ----

#Run the NewCytoscapeMaker script
#Copy the Driver column over to the local gmt file 
#Make sure to remove the first row holding the column names before copying
#Make a new Enrichment map (At either 0.1 cutoff or 1 cutoff, see below)
#Input the gmt file as a shared file
#Input cytoscape input at dataset 1
#Set the style on the enrichment map page to chartdata = None
#Make style changes
  #Change border paint mapping
    #Set mapping to fdr_qvalue
    #Set mapping type to continuous
    #Set min to 0 and max to 0.1
    #Set 0.1 side to white
      #From the black-to-white column, not the orange column (the orange column top box is not pure white)
    #Set 0.05 to mid-orange (Two down from top)
    #Set 0 to dark orange
  #Change border width to 7.0
  #Change Fill color settings 
    #Set to be based on Discrete mapping
    #Set to be based on description
    #Use colors based on phenotypes 
  #Change label to be based on Name
  #Set shape to be based on driver 
    #Set to continuous mapping
    #Set to be based on colouring (stand in for phenotype) 
    #Set min to -1
    #Set max to 1
    #Add node (this defaults to zero, which is correct) 
    #Set the shape on one side to different from the other side
  #Use auto-annotate to generate clusters
    #Using a p value cutoff of 1 (all sets) in the original enrichmentMap run is useful for positioning pathways reasonably. 
    #However, the annotations produced are much less helpful.
  #Manually rename autogenerated clusters based on cluster contents 
          