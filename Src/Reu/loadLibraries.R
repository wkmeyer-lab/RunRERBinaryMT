controlFileLocation = "Src/Loc/dev/controlFile.R"
loadLibraries = function(returnClusterRun = T, install = F, useVisualizePackages = F, additionalLibraryLocations = NULL, loadAllScripts = F){
  
  #This line of code checks a global variable which can be set to let it know that libraries have already been loaded, to prevent re-runs on mulitple scripts.
  
  if(!isTRUE(.GlobalEnv$globalHasLoadedLibrariesAlready)){
    runningLocally=identical(Sys.getenv("RSTUDIO"), "1")
    clusterRunLibLoad = F
    if(!runningLocally){clusterRunLibLoad = T}
    
    newLibraryLocations = c("/share/ceph/wym219group/shared/libraries/R4") # add the Meyer cluster as a default additional location, so that dir.exists works properly
    
    if(!is.null(additionalLibraryLocations)){ #add path to custom libraries to searched locations
      newLibraryLocations = append(newLibraryLocations, additionalLibraryLocations)
    }
    
    .libPaths(c( #set the library paths to the new locations, in addition to any existing library locations. 
      newLibraryLocations[dir.exists(newLibraryLocations)],
      .libPaths()
    ))
    
    
    if(returnClusterRun){
      .GlobalEnv$clusterRun = clusterRunLibLoad
    }
    
    # --- Packages ---
    
    # List of required libraries
    #Main packages
    mainPackages = c(
      "RERconverge",
      "tools",
      "data.table",
      "purrr",
      "dplyr"
    )
    requiredPackages = mainPackages
    
    
    #localOnly packages 
    localyOnlyPackages = c(
      "xlsx"
    )
    if(!clusterRunLibLoad){
      requiredPackages = append(requiredPackages, localyOnlyPackages)
    }
    
    #Visualization Packages
    visualizationPackages = c(
      "ggtree",
      "ggimage", 
      "rphylopic",
      "utils",
      "tidytree",
      "ape",
      "gridExtra",
      "ggpmisc",
      "ggvenn",
      "ggrepel",
      "stringr",
      "tibble"
    )
    if(useVisualizePackages){
      requiredPackages = append(requiredPackages, visualizationPackages)
    }
    
    
    # Check which packages are not installed
    missingPackages = requiredPackages[
      !requiredPackages %in% installed.packages()[, "Package"]
    ]
    
    # Install missing packages
    
    
    
    if(length(missingPackages) > 0){
      if(install){
        install.packages(missingPackages)
      }else{
        stop(
          "Missing required packages: ",
          paste(missingPackages, collapse = ", ")
        )
      }
    }
    
    
    # Load all libraries
    invisible(lapply(requiredPackages, library, character.only = TRUE))
    
    # ---- Loading R scripts ---- 
    #Only the most core scripts are loaded here, to make other dependencies clear in the relevant script
    #However, any other scripts which are sourced by the included scripts are included here, as a centralized resource, but not laoded by default. 
    source("Src/Reu/cmdArgImport.R")
    source(controlFileLocation)
    if(loadAllScripts){
      #Main additional scripts
      source("Src/Reu/ZonomNameConvertVectorCommon.R")
      source("Src/Reu/trimGmtFilenameToSetname.R")
      
      #conditional Scripts:
      {
        #Conditional on running pruning:
        source("Src/Reu/autoPruner.R") #Make[continuous/Categorical/Binary]PhenotypeTree.R
        #Conditional on alternate species generation:
        source("Src/Reu/CategoricalDropTip.R") #MakeCategoricalTree.R
        #Conditional on using relaxed Cateogrical Permulations:
        source("Src/Reu/RelaxedRejectionPermFuncs2.R") #CategoricalPermulationsGeneration.R
      }
    }
  }else{
    message("Global variable indicating pre-loaded libraries, skipping library loading.")
  }

}

