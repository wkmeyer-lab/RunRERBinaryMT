#This function will combine the list of the intermediates produced by CategoricalPermulationGetCor; to be fed into CategoricalCalculatePermulationPValues

combineCategoricalPermulationIntermediates = function(intermediate1, intermediate2){
  intermediatesCombined = intermediate1
  intermediatesCombined[[1]] = cbind(intermediatesCombined[[1]], intermediate2[[1]])
  for(i in 1:length(intermediatesCombined[[2]])){
    intermediatesCombined[[2]][[i]] = cbind(intermediatesCombined[[2]][[i]], intermediate2[[2]][[i]])
  }
  intermediatesCombined[[3]] = cbind(intermediatesCombined[[3]], intermediate2[[3]])
  for(i in 1:length(intermediatesCombined[[4]])){
    intermediatesCombined[[4]][[i]] = cbind(intermediatesCombined[[4]][[i]], intermediate2[[4]][[i]])
  }
  return(intermediatesCombined)
}