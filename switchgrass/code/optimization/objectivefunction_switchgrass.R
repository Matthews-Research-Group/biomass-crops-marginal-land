objectivefunction_switchgrass <- function (parms,soybean_optsolver,observation_switchgrass_IL,wt)
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
    Rhizome_pred = avgrhizome
    Root_pred = avgroot
   print(Aboveground_pred)
   #hourly to daily
    Aboveground_pred = colMeans(matrix(as.vector(Aboveground_pred), nrow=24))
    Rhizome_pred     = colMeans(matrix(as.vector(Rhizome_pred), nrow=24))
    Root_pred        = colMeans(matrix(as.vector(Root_pred), nrow=24))
   print(Aboveground_pred)
   print(observation_switchgrass_IL$Aboveground)
   #error 
   AbovegroundE = (observation_switchgrass_IL$Aboveground - Aboveground_pred)
   RootE        = (observation_switchgrass_IL$Root        - Root_pred)
   RhizomeE     = (observation_switchgrass_IL$Rhizome     - Rhizome_pred)
   
   E = sum(AbovegroundE^2)*wt$aboveground + sum(RootE^2)*wt$root + sum(RhizomeE^2)*wt$rhizome
    
   #try avoid overestimation of the LAST aboveground by adding a penalty
   penalty =  max(0,tail(Aboveground_pred - observation_switchgrass_IL$Aboveground,1))
   E = E + 20 * penalty^2

   if(is.na(E)) {
     print(parms) 
     E = 9999  
   }

#   diff = result$kLeaf - result$kStem
#   if(length(which(diff>0)) > 100) E=9999 
 
   if(max(result2006$leaf) < 3) E=9999

   return(E)
}
