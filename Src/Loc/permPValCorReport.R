permPValCorReport = function (realcor, permvals) {
  timeStart = Sys.time()
  permcor = permvals$corRho                                                     #Rho values from the permulations
  timePermExtractEnd = Sys.time()
  message("Time to extract permulation values: ", timePermExtractEnd - timeStart, attr(timePermExtractEnd - timeStart, "units"))
  
  realstat = realcor$Rho                                                        #Rho values from the real data 
  timeStatExtractEnd = Sys.time()
  message("Time to extract real values: ", timeStatExtractEnd - timePermExtractEnd, attr(timeStatExtractEnd - timePermExtractEnd, "units"))
  names(realstat) = rownames(realcor)                                           #Name the values after the row names 
  
  permcor = permcor[match(names(realstat), rownames(permcor)),]                 #trim the permulations' Rho values to only ones in the real data 
  timeMatchEnd = Sys.time()
  message("Time to match names: ", timeMatchEnd - timeStatExtractEnd, attr(timeMatchEnd - timeStatExtractEnd, "units"))
  
  permpvals = vector(length = length(realstat))                                 #Make a vector to store the pValues of a length equal to the number of realstat genes
  names(permpvals) = names(realstat)                                            #Copy the names from the real Rho values onto that vector
  count = 1                                                                     #Set a count equal to one 
  
  timeBeforeLoopStart = Sys.time()
  message("Total time before loop start: ", timeBeforeLoopStart - timeStart, attr(timeBeforeLoopStart - timeStart, "units"))
  
  while (count <= length(realstat)) {                                           #While the count hasn't done all of the entries                                         
    timeLoopBegin = Sys.time()
    message("Gene number: ", count)
    if (is.na(realstat[count])) {                                               #if the real correlation is NA
      permpvals[count] = NA                                                     #Make the permulation value NA
      message("is NA")
    }
    else {                                                                      #If not 
      num = sum(abs(permcor[count, ]) > abs(realstat[count]),                   #Make a numerator which is the sum of the permulated values greater than the actual values; after removing NA entires. 
                na.rm = T)
      timeNumeratorCalc = Sys.time()
      message("Numerator Calculate time: ", timeNumeratorCalc - timeLoopBegin, attr(timeNumeratorCalc - timeLoopBegin, "units"))
      
      denom = sum(!is.na(permcor[count, ]))                                     #Make a denominator which is the sum of the non-NA permulation values
      timeDenomSum = Sys.time()
      message("Denominator sum time: ", timeDenomSum - timeNumeratorCalc, attr(timeDenomSum - timeNumeratorCalc, "units"))
      permpvals[count] = num/denom                                              #p-value = numerator/denominator
    }
    timeLoopEnd = Sys.time()
    message("Total loop time: ", timeLoopEnd - timeLoopBegin, attr(timeLoopEnd - timeLoopBegin, "units"))
    count = count + 1                                                           #increase count by one
  }
  message("calculation compelte.")
  timeEnd = Sys.time()
  message("Total calculation time: ", timeEnd - timeStart, attr(timeEnd - timeStart, "units"))
  permpvals                                                                     #Return the vector of pValues 
}