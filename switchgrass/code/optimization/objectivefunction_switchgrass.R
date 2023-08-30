objectivefunction_switchgrass <- function (parms,soybean_optsolver,optsolver_site2,observation_switchgrass_IL,observation_switchgrass_ND,wt)
{
         
  date_2006 = c("2006-04-15","2006-06-15","2006-08-15","2006-10-15")
  date_2006 = as.Date(date_2006)
  doy_2006  = as.numeric(format(date_2006,"%j"))
  date_2007 = c("2007-04-15","2007-06-15","2007-08-15","2007-10-15")
  date_2007 = as.Date(date_2007)
  doy_2007  = as.numeric(format(date_2007,"%j"))
  date_2008 = c("2008-04-15","2008-06-15","2008-08-15","2008-10-15")
  date_2008 = as.Date(date_2008)
  doy_2008  = as.numeric(format(date_2008,"%j"))
  doy_all = list(doy_2006,doy_2007,doy_2008)

  solver <- soybean_optsolver[[1]]
  result <- solver(parms)
  result2006 <- data.frame(doy=result$doy,stem = result$Stem, leaf = result$Leaf, leaflitter = result$LeafLitter, 
                           rhizome = result$Rhizome, root = result$Root)
  solver <- soybean_optsolver[[2]]
  result <- solver(parms)
  result2007 <- data.frame(doy=result$doy,stem = result$Stem, leaf = result$Leaf, leaflitter = result$LeafLitter, 
                           rhizome = result$Rhizome, root = result$Root)
  solver <- soybean_optsolver[[3]]
  result <- solver(parms)
  result2008 <- data.frame(doy=result$doy,stem = result$Stem, leaf = result$Leaf, leaflitter = result$LeafLitter, 
                                rhizome = result$Rhizome, root = result$Root)

#get the values at the observed DOY. There are four days excluding Dec.
#average over three years,but it's still hourly values
    avgstem =  (result2006$stem[result2006$doy %in% doy_all[[1]]] +
                result2007$stem[result2007$doy %in% doy_all[[2]]] +
                result2008$stem[result2008$doy %in% doy_all[[3]]] )*(1/3)
    
    avgleaf =  (result2006$leaf[result2006$doy %in% doy_all[[1]]] +
                result2007$leaf[result2007$doy %in% doy_all[[2]]] +
                result2008$leaf[result2008$doy %in% doy_all[[3]]] )*(1/3)
       
    avgleaflitter =  (result2006$leaflitter[result2006$doy %in% doy_all[[1]]] +
                      result2007$leaflitter[result2007$doy %in% doy_all[[2]]] +
                      result2008$leaflitter[result2008$doy %in% doy_all[[3]]] )*(1/3)
    
    avgrhizome =  (result2006$rhizome[result2006$doy %in% doy_all[[1]]] +
                   result2007$rhizome[result2007$doy %in% doy_all[[2]]] +
                   result2008$rhizome[result2008$doy %in% doy_all[[3]]])*(1/3)
    
    avgroot =     (result2006$root[result2006$doy %in% doy_all[[1]]] +
                   result2007$root[result2007$doy %in% doy_all[[2]]] +
                   result2008$root[result2008$doy %in% doy_all[[3]]])*(1/3)

    Aboveground_pred = avgstem + avgleaf #+ avgleaflitter
    Leaf_pred    = avgleaf
    Stem_pred    = avgstem
    Rhizome_pred = avgrhizome
    Root_pred    = avgroot
   #hourly to daily to match OBS's
    Aboveground_pred = colMeans(matrix(as.vector(Aboveground_pred), nrow=24))
    Leaf_pred        = colMeans(matrix(as.vector(Leaf_pred), nrow=24))
    Stem_pred        = colMeans(matrix(as.vector(Stem_pred), nrow=24))
    Rhizome_pred     = colMeans(matrix(as.vector(Rhizome_pred), nrow=24))
    Root_pred        = colMeans(matrix(as.vector(Root_pred), nrow=24))
   #error 
   LeafE        = (observation_switchgrass_IL$Leaf        - Leaf_pred)
   StemE        = (observation_switchgrass_IL$Stem        - Stem_pred)
   RootE        = (observation_switchgrass_IL$Root        - Root_pred)
   RhizomeE     = (observation_switchgrass_IL$Rhizome     - Rhizome_pred)
   
   E1 = sum(LeafE^2)*wt$leaf + sum(StemE^2)*wt$stem + sum(RootE^2)*wt$root + sum(RhizomeE^2)*wt$rhizome

   E1 = E1/length(LeafE)
    
   #try avoid overestimation of the LAST aboveground by adding a penalty
   penalty =  max(0,tail(Aboveground_pred - observation_switchgrass_IL$Aboveground,1))
   E1 = E1 + 20 * penalty^2

   if(is.na(E1)) {
     print(parms) 
     E1 = 9999  
   }

#get the second site's error
   obs_doys = observation_switchgrass_ND$doy
   obs_stem = observation_switchgrass_ND$stem
   obs_leaf = observation_switchgrass_ND$leaf
   solver <- optsolver_site2[[1]] #2001
   result2001 <- solver(parms)
   solver <- optsolver_site2[[2]] #2002
   result2002 <- solver(parms)

   avgstem =  (result2001$Stem[result2001$doy %in% obs_doys] +
               result2002$Stem[result2002$doy %in% obs_doys])/2
    
   avgleaf =  (result2001$Leaf[result2001$doy %in% obs_doys] +
               result2002$Leaf[result2002$doy %in% obs_doys])/2
   #hourly to daily
   avgstem_pred = colMeans(matrix(as.vector(avgstem), nrow=24))
   avgleaf_pred = colMeans(matrix(as.vector(avgleaf), nrow=24))

   stemE = avgstem - obs_stem
   leafE = avgleaf - obs_leaf
   E2 = sum(stemE^2) + sum(leafE^2)
   E2 = E2/length(stemE)
   #print(c(E1,E2))
   E = E1+E2
#   diff = result$kLeaf - result$kStem
#   if(length(which(diff>0)) > 100) E=9999 

   #be careful! here is hard-coded indices, corresponding to kLeaf_emr and kStem_emr
   if(parms[8]+parms[9]>1) E=9999  

   if(max(result2006$leaf) < 2) E=9999

   return(E)
}
