clusterRun = F
clusterRun = T
if(clusterRun){.libPaths("/share/ceph/wym219group/shared/libraries/R4")} #add path to custom libraries to searched locations
library(jsonlite)
library(dplyr)
source("Src/Reu/cmdArgImport.R")

# -- Command arguments list
# r = filePrefix                                                               This is a prefix used to organize and separate files by analysis run. Always required. 
# g = geneName                                                                 This sets the gene for hyphy to be run on 
# f = foregroundValue                                                          This is the foreground value used for file name matching
# t = hyphyTool                                                               This is the tool used for file matching

args = c("r=CategoricalInsvertivoreTree", "g=EHHADH", "f=2", "t=absrel")

if(clusterRun)args = commandArgs(trailingOnly = TRUE)
filePrefix = NULL
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
#Make hyphySpecific output folder
outputHyphyFolderNameNoSlash = paste0(outputFolderName, "Hyphy")
if(!dir.exists(outputHyphyFolderNameNoSlash)){                       #create that directory if it does not exist
  dir.create(outputHyphyFolderNameNoSlash)
}
outputHyphyFolderName = paste0(outputHyphyFolderNameNoSlash, "/")

# --- Argument Imports ---
# Defaults
geneName = NULL
foregroundValue = NULL
hyphyTool = NULL

#Target gene
if(!is.na(cmdArgImport('g'))){
  geneName = cmdArgImport('g')
}else{
  stop("THIS IS AN ISSUE MESSAGE; SPECIFY TARGET GENE")
}

#ForegroundValue
if(!is.na(cmdArgImport('f'))){
  foregroundValue = cmdArgImport('f')
}else{
  stop("THIS IS AN ISSUE MESSAGE; SPECIFY ForegroundValue")
}

#hyphy tool
if(!is.na(cmdArgImport('t'))){
  hyphyTool = cmdArgImport('t')
}else{
  stop("THIS IS AN ISSUE MESSAGE; SPECIFY HYPHY TOOL")
}


hyphyFileNoExtension = paste0(outputHyphyFolderName, filePrefix, "-Hyphy-", hyphyTool, "-", geneName, "-Foreground_", foregroundValue)
hyphyFile = paste0(hyphyFileNoExtension, ".json")

# --- Cody body, from chat GPT --- 

json_data <- fromJSON(hyphyFile)

# Branch for different HyPhy tools
if(hyphyTool == "absrel"){
	# Extract branch attributes
  branch_data <- json_data$`branch attributes`$`0`
  
  # Convert to data frame
	df <- do.call(rbind, lapply(names(branch_data), function(node) {
  	data.frame(Node = node,
             Uncorrected_P_value = ifelse(is.null(branch_data[[node]]$`Uncorrected P-value`), NA, branch_data[[node]]$`Uncorrected P-value`),
             Corrected_P_value = ifelse(is.null(branch_data[[node]]$`Corrected P-value`), NA, branch_data[[node]]$`Corrected P-value`))
	}))
	colnames(df) = paste0(hyphyTool, "_", colnames(df))

	df = df[order(df$absrel_Corrected_P_value),]

	outFileName = paste0(hyphyFileNoExtension, ".csv")

	# Save as CSV
	write.csv(df, outFileName, row.names = FALSE)

	# Print message
	cat(paste("CSV file with p-values has been saved as '", outFileName, "'\n"))
}
if(hyphyTool == "relax"){
	#p-value is inside 'test results'.
	p_value <- json_data[['test results']][['p-value']]

	#extract Omega (ω) values
	#these are under 'fits' -> 'Alternative' -> 'rate distributions'.
  omega_reference <- json_data$fits[['RELAX alternative']][['Rate Distributions']][['Reference']]
    
    omega_test <- json_data$fits[['RELAX alternative']][['Rate Distributions']][['Test']]

	#extract the Tree File (Newick String)
	tree_string <- json_data[['branch attributes']][['0']]

	#combine all extracted data into a single list
	output_data <- list(
	  p_value = p_value,
    omega_reference = omega_reference, 
    omega_test = omega_test,           
    newick_tree_string = tree_string
	)	

	#define output file path
	rds_file_path <- paste0(hyphyFileNoExtension, ".Rds")

	#save the list object to the .Rds file
	saveRDS(output_data, file = rds_file_path)

	#print confirmation message
	print(paste0("All data successfully saved to:", rds_file_path))
}else{
	print(paste0("The ", hyphyTool, " tool has not yet been implemented into this pipeline"))
}
