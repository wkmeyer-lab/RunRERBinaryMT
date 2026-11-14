#This is a control file, designed to help make inputting arguments easier. 


# ---- REQUIRED ARGUMENTS ----
#These are the arguments which are STRICTLY REQUIRED to run the pipeline. 
#They are currently set to the demo used in the walk through. You will need to change these each analysis. 

requiredArgs = c(
    #Input your File Prefix (the name of the analysis) two lines below
  'r=     
    Demo
  ',  #Input your Analysis Type two lines below
  'a=    
    Categorical
  ',  #Input your Phenotype Data Spreadsheet Filepath two lines below
  'td=    
    Data/examplePhenotypeData.csv
  ',  #Input your Phenotype Data Column's Name two lines below
  'tp=    
    myPhenotype
  ',  #Input your Species' Names In Alignments ( or On The Tips Of Tree) two lines below
  'tt=    
    toyTips
  ',  #Input your Geneset .gmt File Filepaths two lines below
  'g=    
    c("gene","Data/exampleGeneset1.gmt","Data/exampleGeneset2.gmt")
  ',  #Input your MainTrees File Path two lines below
  'm=    
    Data/toyTrees.rds
  '
)

# -- CONTROL FILE FOR CLUSTER --
#If you want the values from the control file to be used by a cluster-run script:
#Simply run scripts on a cluster without providing any trailing arguments. 
#Alternitavely, set the following boolean to TRUE, 
forceControlFileCluster = FALSE
forceControlFileAlways = FALSE


#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#If you simply want to run the most basic version of the pipeline, feel free to stop here. 
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------








#Template to add new arg to one of the arg vectors:
#,  #Input your [] two lines below
#  '[]=    
#    []
#  '



# ---- FREQUENTLY USED ARGUMENTS ----
#These are arguments which a user is likely to want to change. 
#Any arguments with a value of NA act as though the arugment was not provided; the default behavior is used. 
#In order for these arguments to be proccessed by this control file, the boolean on the next line must be set to TRUE
useControlFileCommonArgs = TRUE

#commonAnalysisArgs
{ # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
commonAnalysisArgs = c(
     #Input your Force Update value two lines below
  'f=    
    FALSE
  ',  #Input your Minimum Species For Gene To Be Included In RER Analysis Value two lines below
  'l=    
    NA
  ', #Input your continuous Metric two lines below
  'k=    
    diff
  ',  #Input your Subdirectories To Use In The Analysis two lines below
  's=    
    NA
  ',  #Input your UsedAlternates Value two lines below (Unless the number of alternates is less than 10, recommended to NOT use this on a local machine.)
  'i=    
    all
  ',  #Input your Categories To Have A Second Analysis Run On Their Combination two lines below
  'e=    
    NA
  ',  #Input your Name Of The Combined Phenotype two lines below
  't=    
    NA
  '
)
}

#commonTreeArgs
{ # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
commonTreeArgs = c(
     #Input your mainTreesLocation two lines below
  'm=    
    Data/toyTrees.rds
  ',  #Input your Color Of Phenotypes In Tree Plot (ordered alphabetically) two lines below
  'tb=    
    c("darkgreen","darkblue","black","red","yellow")
  ',  #Input your Phenotype Data Spreadsheet's Common Name Column's Name two lines below
  'tn=    
    NA
  ',  #Input your Set Of Categories To Include In The Analysis (species with any other categories not included in the tree) two lines below
  'tc=    
    NA
  ',  #Input your Phenotype Data Spreadsheet's Screen Column's Name (species must have a value of 1 in this column, or they are excluded) two lines below
  'ts=    
    NA
  ',  #Input your Generate Alternate Species Sets Boolean two lines below
  'tl=    
    T
  ',  #Input your Number Of Generated Alternate Species Sets two lines below
  'ta=    
    5
  ',  #Input your Ancestral State Reconstruction Model two lines below
  'q=    
    NA
  ',  #Input your Minimum Branch Length For Automatic Pruning (removes branches until all are longer than this value) two lines below
  'pm=    
    NA
  ',  #Input your Tip Names Of Species To Be Pruned From The Tree two lines below
  'pd=    
    NA
  '
)
}



# ---- RARELY USED ARGUMENTS ----
#These are arguments which are rarely used when launching via control file. 
#Any arguments with a value of NA act as though the arugment was not provided; the default behavior is used. 
#In order for these arguments to be proccessed by this control file, the boolean on the next line must be set to TRUE
useControlFileRareArgs = FALSE

#rareAnalysisArgs
{ # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
  rareAnalysisArgs = c(
    #Input your Print Arguments Instead Of Run Boolean two lines below (This arg is here, because it is rarely used via control file.)
    'h=    
    NA
  ',  #Input your writeData Boolean two lines below
  'w=    
    NA
  ',  #Input your Alternates Ratio For Robustness Value two lines below
  'o=    
    NA
  ',  #Input your Significance Cutoff Value two lines below
  'c=    
    NA
  '
  )
}

#rareTreeArgs
{ # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
rareTreeArgs = c(
     #Input your Forced Interal Branch Category two lines below
  'ta=    
    NA
  ',  #Input your Use Liam Inference Boolean two lines below
  'tz=    
    NA
  ',  #Input your Tips To Block From Automatic Pruning two lines below
  'ps=    
    NA
  ',  #Input your Pruning Preference Column's Name (species with a value of 1 in this column will be preferentially kept during automatic pruning) two lines below
  'pk=    
    NA
  ',  #Input your List Of Categories To Be Replaced With The Second Value In The List two lines below
  'tr=    
    list(c("replaceExample1", "with1"),c("replaceExample2, with2"))
  ',  #Input your List Of Categories To Have Combinations Combined two lines below
  'tm=    
    list(c("example1", "intoExample1"), c("example2", "intoExample2"))
  '
)
}


# ---- PERMULATIONS ARGUMENTS ----
#These are arguments which are used to control permulations. 
#Permulations are time-impractical without parrallelization on a cluster, which requires job scheduling. 
#As such, it is STRONGLY ADVISED that you do not use a control file for permulation scripts. They are included here for testing and completeness.
#Any arguments with a value of NA act as though the arugment was not provided; the default behavior is used. 
#In order for these arguments to be proccessed by this control file, the boolean on the next line must be set to TRUE
useControlFilePermulationArgs = FALSE

#permulationArgs
{ # This bracket allows this code to be collapsed in Rstudio; it serves no code purpose. 
permulationArgs = c(
  #Input your Use Permulations Value two lines below
  'p=    
    NA
  ',  #Input your Number Of Permulations Per Script Instance Value two lines below
  'n=    
    3
  ',  #Input your Permulations Run Instance Value two lines below
  'j=    
    LocalDebugTest
  ',  #Input your Root State Probability two lines below
  'u=    
    NA
  ',  #Input your Permulation Combination Start Value two lines below
  'st=    
    NA
  ',  #Input your Continuous Permualtions Outgroup Tip Name two lines below
  'po=    
    NA
  ',  #Input your Categorical Permulation Relaxation Value two lines below
  'x=    
    NA
  ',  #Input your Permulations Metacombination Boolean two lines below (NOTE THAT THIS ARGUMENT WILL NOT BEHAVE CORRECTLY VIA CONTROL FILE IF SET TO TRUE)
  'b=    
    F
  ',  #Input your Permualtions Calculate P-Value Boolean two lines below
  'd=    
    NA
  ',  #Input your Permualtions Only Calculate P-Value Boolean two lines below
  'v=    
    NA
  '
)
}




# ---------------------------------------------
# --- Process Above Values Into Args Vector ---
# ---------------------------------------------

#This function converts the control file format above to the required input format
argPolishing = function(vector){
  naIndexes = which(grepl("=\\s*NA\\s*$", vector, ignore.case = TRUE))
  if(length(naIndexes >0)){vector = vector[-naIndexes]}
  
  
  vector = gsub("[\t\r\n]", "", vector)
  vector = gsub("=\\s+", "=", vector)
  vector = trimws(vector)
  return(vector)
}

controlFileArgs = argPolishing(requiredArgs)

if(useControlFileCommonArgs){
  commonAnalysisArgs = argPolishing(commonAnalysisArgs)
  commonTreeArgs = argPolishing(commonTreeArgs)
  
  controlFileArgs = append(controlFileArgs, commonAnalysisArgs)
  controlFileArgs = append(controlFileArgs, commonTreeArgs)
}

if(useControlFileRareArgs){
  rareAnalysisArgs = argPolishing(rareAnalysisArgs)
  rareTreeArgs = argPolishing(rareTreeArgs)
  
  controlFileArgs = append(controlFileArgs, rareAnalysisArgs)
  controlFileArgs = append(controlFileArgs, rareTreeArgs)
}

if(useControlFilePermulationArgs){
  permulationArgs = argPolishing(permulationArgs)

  controlFileArgs = append(controlFileArgs, permulationArgs)
}

#Remove any duplicated args 
repeatArgs = which(duplicated(trimws(sub("=.*", "", controlFileArgs)), fromLast = FALSE)) #This marks the value further down the file as the duplicate
if(length(repeatArgs)>0){controlFileArgs = controlFileArgs[-repeatArgs]}



#If one of the force variables was toggled on, create the stable forcing marker. 
if(forceControlFileCluster){
  controlFileForceCluster = "exists"
}
if(forceControlFileAlways){
  controlFileForceAlways = "exists"
}


#Clean up all of the non-output vairables made by this file 
rm(forceControlFileCluster)
rm(forceControlFileAlways)
rm(requiredArgs)
rm(useControlFileCommonArgs)
rm(commonAnalysisArgs)
rm(commonTreeArgs)
rm(useControlFileRareArgs)
rm(rareAnalysisArgs)
rm(rareTreeArgs)
rm(useControlFilePermulationArgs)
rm(permulationArgs)
rm(argPolishing)
rm(repeatArgs)




controlFileArgs