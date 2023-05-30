il_objfunlogistic <- function (parms,partial_gro_list,observed){
  
  #these predicted values will be in vectors
  #although they are initialized as a single value
  predicted = list()
  predicted$stem = 0
  predicted$root = 0
  predicted$rhizome = 0
  predicted$leaf = 0
  for (i in 1:length(partial_gro_list)){
    partial_gro_function = partial_gro_list[[i]]
    current_res <- partial_gro_function(parms)
    predicted$stem    = predicted$stem    + current_res$Stem[observed$resultindex_for_optim]
    predicted$root    = predicted$root    + current_res$Root[observed$resultindex_for_optim]
    predicted$rhizome = predicted$rhizome + current_res$Rhizome[observed$resultindex_for_optim]
    predicted$leaf    = predicted$leaf    + current_res$Leaf[observed$resultindex_for_optim]
  }
  #3-year averages
  predicted$stem    = predicted$stem/length(partial_gro_list)
  predicted$root    = predicted$root/length(partial_gro_list)
  predicted$rhizome = predicted$rhizome/length(partial_gro_list)
  predicted$leaf    = predicted$leaf/length(partial_gro_list)
  
# First data point of root should be zero. Observation does not make difference between alive and dead roots
  observed$root[1] = 0

  stemE    = observed$stem - predicted$stem 
  rootE    = observed$root - predicted$root 
  rhizomeE = observed$rhizome - predicted$rhizome 
  leafE    = observed$leaf - predicted$leaf 

#this weights apply to the DOYs' of errors where the last DOY is assigned a low weight
  E = sum(stemE^2* observed$wt)  + sum(rootE^2* observed$wt) + sum(rhizomeE^2* observed$wt)  + sum(leafE^2* observed$wt)

  if(is.na(E)) {
    print(parms)
    E=999
  }
 #if kLeaf_emr + kStem_emr >1, discard it 
  if(parms[2] + parms[3] >1 ) E=999
 return(E)
}
