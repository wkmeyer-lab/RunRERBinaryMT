library(RERconverge)
source("Src/Loc/dev/MergeDataCreation.R")


# The first thing we need is a set of genetrees with a mastertree based on the average of the gene trees. This set of trees is called a Trees Object, or a MainTrees. 
# In this case, we will be using the demo shiped with RERConerge. 

rerPath = find.package("RERconverge")
demoTreeFile = paste(rerPath, "/extdata/", "SubsetMammalGeneTrees.txt", sep="")

demoTrees = readTrees(demoTreeFile)

demoTrees$masterTree$tip.label


# The second thing we needs is a phenotype, and a way to connect the phenotype data with the tips on the tree. This is called a MergeData. 
# We will be creating one now. 
# MergeData have a few critical components:
  #The species scientific name
  # The species name as it appears on the tips in the tree
  # Phenotype data for species. 
#In this walkthrough, you have been provided diet phenotype data for the species in the tree. 

#MergeData is made using the "CombineDatasets" script. 
#CombineDatasets requires two things: a column with scientific names of the species, and a column of data you want to add. 
#Your final mergeData requires at least the following: 
  #A column of scientific names
  #A column of common names
  #A column of names as they appear in the Maintrees
  #A column containing your phenotype data. 


#The first table to add into CombineDatasets is the scientific names of species. 
ScientificName = c("Ailuropoda melanoleuca", "Allactaga bullata", "Bos taurus", 
                           "Callithrix jacchus", "Camelus bactrianus", "Canis lupus familiaris", 
                           "Capra hircus", "Cavia porcellus", "Chinchilla lanigera", "Chlorocebus sabaeus", 
                           "Chrysochloris asiatica", "Condylura cristata", "Cricetulus griseus", 
                           "Dasypus novemcinctus", "Delphinus delphis", "Diceros bicornis", 
                           "Echinops telfairi", "Elephantulus edwardii", "Eptesicus fuscus", 
                           "Equus caballus", "Erinaceus europaeus", "Felis catus", "Gorilla gorilla", 
                           "Heterocephalus glaber", "Homo sapiens", "Loxodonta africana", 
                           "Macaca fascicularis", "Macaca mulatta", "Macropus eugenii", 
                           "Mesocricetus auratus", "Monodelphis domestica", "Mus musculus", 
                           "Mustela putorius", "Myotis davidii", "Myotis lucifugus", "Nomascus leucogenys", 
                           "Ochotona princeps", "Octodon degus", "Odobenus rosmarus", "Orcinus orca", 
                           "Ornithorhynchus anatinus", "Orycteropus afer afer", "Oryctolagus cuniculus", 
                           "Otolemur garnettii", "Ovis aries", "Pan troglodytes", "Pantholops hodgsonii", 
                           "Papio anubis", "Pongo abelii", "Pteropus vampyrus", "Rattus norvegicus", 
                           "Saimiri boliviensis", "Sarcophilus harrisii", "Sorex araneus", 
                           "Spermophilus tridecemlineatus", "Sus scrofa", "Trichechus manatus latirostris", 
                           "Tupaia chinensis", "Vicugna pacos", "Halichoerus grypus", "Microtus arvalis", 
                           "Pteropus giganteus")
SpeciesScientificNamesDataframe = data.frame(ScientificName)


#We then use the createSortedTable function to create a table with a column with a consistently formatted scientific name that can be used to connect across datasets (explain space, sub spcies, ect)

starterData = createSortedTable(SpeciesScientificNamesDataframe, "ScientificName")

#Next, we add the common names 
CommonName = c("Panda", "Gobi jerboa", "Cow", "Marmoset", "Bactrian Camel", 
                       "Dog", "Domestic goat", "Domestic guinea pig", "Chinchilla", 
                       "Green monkey", "Cape golden mole", "Star-nosed mole", "Chinese hamster", 
                       "Armadillo", "Dolphin", "Black Rhinoceros", "Tenrec", "Cape elephant shrew", 
                       "Big brown bat", "Horse", "Hedgehog", "Cat", "Gorilla", "Naked mole-rat", 
                       "Human", "Elephant", "Crab-eating macaque", "Rhesus Macaca", 
                       "Wallaby", "Golden hamster", "Opossum", "Common mouse", "Ferret", 
                       "David's myotis bat", "Microbat", "Gibbon", "Pika", "Brush-tailed rat", 
                       "Pacific walrus", "Killer whale", "Platypus", "Aardvark", "Rabbit", 
                       "Bushbaby", "Domestic sheep", "Chimp", "Tibetan antelope", "Olive Baboon", 
                       "Orangutan", "Megabat", "Rat", "Squirrel monkey", "Tasmanian devil", 
                       "Shrew", "Squirrel", "Pig", "Manatee", "Chinese tree shrew", 
                       "Alpaca", "Grey seal", "Common vole", "Indian flying fox")

SpeciesCommonNamesDataframe = data.frame(CommonName, ScientificName) #Note that all added data must have scientific name included in the added data


#### DEV NOTE: Figure out why a global combinedData has to exist for this to work. #######
#We then combine these into a single dataset using the CombineDatsets Function. By listing CommonName as a Name column, it is put before the divider between name columns and data columns. 
combinedData = starterData
combinedData = CombineDatasets(combinedDataInput = starterData, newDatasetInput = SpeciesCommonNamesDataframe, newDataScientificNameColumn = "ScientificName", addNewSpeciesValue = T, nameColumns = "CommonName")

#We then need to add phenotype data for our species to the mergedData 

DemoDietPhenotype = c("Omnivore", "Omnivore", "Herbivore", "Omnivore", "Omnivore", 
                      "Carnivore", "Herbivore", "Herbivore", "Herbivore", "Omnivore", 
                      "Carnivore", "Carnivore", "Herbivore", "Carnivore", NA, "Herbivore", 
                      NA, "Carnivore", "Carnivore", "Herbivore", "Omnivore", "Carnivore", 
                      "Herbivore", "Herbivore", "Omnivore", "Herbivore", "Omnivore", 
                      "Omnivore", "Herbivore", "Omnivore", "Omnivore", "Omnivore", 
                      "Carnivore", "Carnivore", "Carnivore", "Herbivore", "Herbivore", 
                      "Herbivore", "Carnivore", "Carnivore", "Carnivore", "Carnivore", 
                      "Herbivore", "Herbivore", "Herbivore", "Omnivore", "Herbivore", 
                      "Omnivore", "Omnivore", "Herbivore", "Omnivore", "Omnivore", 
                      "Omnivore", "Omnivore", "Omnivore", "Omnivore", "Herbivore", 
                      NA, "Herbivore", NA, NA, NA)

DietDataframe = data.frame(DemoDietPhenotype, ScientificName)
newDataset = DietDataframe
combinedData = CombineDatasets(combinedDataInput = combinedData, newDatasetInput = DietDataframe, newDataScientificNameColumn = "ScientificName", addNewSpeciesValue = T)

#Finally, we need to add the column of the tip names as they appear in MainTrees. 

demoTreeTipName = demoTrees$masterTree$tip.label
tipNameDataframe = data.frame(demoTreeTipName, ScientificName)

combinedData = CombineDatasets(combinedDataInput = combinedData, newDatasetInput = tipNameDataframe, newDataScientificNameColumn = "ScientificName", addNewSpeciesValue = T, nameColumns = "demoTreeTipName")


MergedData = combinedData 

#Here, you would save your mergedData as a result. 
write.csv(MergedData, "Results/DemoMergedData.csv")



#With your dataset and merged data complete, we can now begin using the functions of runRER. 

#The first script to use is MakeCategoricalPhenotypeTree. If you are not using a categorical phenotype, you would use the appropriate phenotype tree creator. 


#These scripts are operated by providing them a list of arguments at the start of the script. This is for compatibility with running on a cluster, which will be covered later. 
#When running a script, the first step is to read the list of args, and what they do. 

#There are a few args which are common to many scripts: 
# r   This inidcates the prefix being used. Each prefix has its own output folder, and allows you to run the script on mulitple analyses without overwriting the other analyses. 
  #In this case, out prefix will be "demo"
# v   This prefix will force the script to de-generate all of its input files, instead of using ones generated in the past. 
  #Useful if you changed something, and need to make updates. If you don't, leave v to False (F), because regenerating files wastes a lot of time. 
#m    This is the filename of your maintrees file. 
  #In this case, this will be the what is stored in the <DemoTreeFile> variaible.

#There are many arguments in Make Categorical Phenotype Tree. other than the three above, here are the ones we will use: 
#d    This is the location of your MergedData Spreadsheet. 
#a    This is the name of the column with the phenotype data. In this case, "DemoDietPhenotype". 
#c    This is a list of the categories you want to include in the analysis. Any categories not listed here will be ignored. This allows the script to ignore species either without phenoytpe data or with phenotypes not relevant to the analysis.
  #In this case, "Carnivore", "Herbivore", "Omnivore". 
#n    This is the column in your spreadsheet with the tip names in MainTrees. 
  #In this case, "demoTreeTipName". 

#The other arguments are involved in remaning phenotype values (u and o), using only a subset of the species in the mainTrees with the relevant phenotypes (s), details of ancestral state reconsturction (t and g) and tree pruning (z, x, y, and p).
  #We can ignore all of these for now -- if you do not include the in your arg string, they will use their default values, which are fine.  


## TEST YOURSELF: What should the args you use for this script be? ##
#####DEV Note: In full version, better explain why args are a string (need to use c), why use ' quotes, why args with mulitple values need another c and use " quotes. ########
##Spoiler this##
#args = c('r=Demo', 'v=F', 'm=demoTreeFile', 'd=Results/DemoMergedData.csv', 'a=DemoDietPhenotype', 'c=c("Carnivore", "Herbivore", "Omnivore")', 'n=demoTreeTipName')
##end spoiler##

#If you run that script, you should find a copy of your tree in "Ouput/Demo/DemoCategoricalTree.pdf
#If you look at it, and the branches are colored, good job, it worked! 
#However, the colors might seem a little mis-aligned. Run the following command, and then run the code that generates the pdf (and ONLY that part of the code), and see if which phenotype has which color has changed. 
palette(c("red", "darkgreen", "black"))




