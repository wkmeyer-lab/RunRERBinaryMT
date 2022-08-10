#Remove the pre-existing gitHUb remote
system('git remote remove origin')

#Load all source files in /reu and /reu/import
file.sources = list.files(c("src/reu", "src/reu/import"), 
                          pattern="*.R$", full.names=TRUE, 
                          ignore.case=TRUE)
#Below code would auto-import the results from file.sources. I've turned it off for now, to avoid cultting the environment. 

#sapply(file.sources,source,.GlobalEnv) 
