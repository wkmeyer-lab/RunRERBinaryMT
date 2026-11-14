# --- Standard start-up code ---

# - Load in args from either the command line or control file -
if(clusterRun){args = commandArgs(trailingOnly = TRUE)}

if(!is.vector(args)){
  message("No ARGs supplied, using control file ARGs.")
  args = controlFileArgs
}else if(exists("controlFileForceAlways")){ #Load args from the control file if none exist, or Always Force is on
  message("Control File always-force is active, using control file ARGs instead of supplied ARGs!")
  args = controlFileArgs
}else if(clusterRun && exists("controlFileForceCluster")){ #Load args from the control file if running on a cluster, and Cluster Force is on
  message("Control File cluster-force is active, using control file ARGs instead of supplied ARGs!")
  args = controlFileArgs
} 


# - Import core arguments used by all scripts -
{  # Bracket used for collapsing purposes
  #Arg help 
  if(isTRUE(as.logical(cmdArgImport('h')))){
    message("Not Running Script, instead returning ARGs for script.")
    printArgDocumentation()
    stop()
  }
  
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
  if(!is.na(cmdArgImport('f'))){                                 #Import if update being forced with argument 
    forceUpdate = cmdArgImport('f')
    forceUpdate = as.logical(forceUpdate)
  }else{
    message("Force update not specified, not forcing update")
  }
  
}
