il_objfunlogistic <- function (parms,partial_gro_list,observed){
  
  #these predicted values will be in vectors
  #although they are initialized as a single value
  predicted = list()
  predicted$stem = 0
  predicted$root = 0
  predicted$rhizome = 0
  predicted$leaf = 0
  obs_doys = observed$doy #get the obs' DOYs
  
  #hard-coded here!
  gro_2006_2008 = partial_gro_list#[c(5:7)]
  
  time_ind = c(1,1273,2713,4201,5665)

  for (i in 1:length(gro_2006_2008))
  {
    partial_gro_function = gro_2006_2008[[i]]
    current_res <- partial_gro_function(parms)
  #  time = current_res$doy + current_res$hour/24
  #  time_ind = which(time%in%obs_doys)
  #  if(length(time_ind)>length(obs_doys)) {
  #    stop(c('doy not matching!',time_ind))
  #  } 
    predicted$stem    = predicted$stem    + current_res$Stem[time_ind]
    predicted$root    = predicted$root    + current_res$Root[time_ind]
    predicted$rhizome = predicted$rhizome + current_res$Rhizome[time_ind]
    predicted$leaf    = predicted$leaf    + current_res$Leaf[time_ind]
  }
  #3-year averages
  predicted$stem    = predicted$stem/length(gro_2006_2008)
  predicted$root    = predicted$root/length(gro_2006_2008)
  predicted$rhizome = predicted$rhizome/length(gro_2006_2008)
  predicted$leaf    = predicted$leaf/length(gro_2006_2008)
  
# First data point of root should be zero. Observation does not make difference between alive and dead roots
  observed$root[1] = 0

#be careful of the variable names of the observation
#if they don't match, the program won't stop!
#you may see warning message like "stack imbalance". Very annoying!
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
