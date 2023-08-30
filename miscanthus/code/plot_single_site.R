library(ggplot2)
Compare_Observed_and_Predicted <- function(observeddata,predicteddata,xlabtitle,ylabtitle){
  
  observeddata$pointlegend <- as.factor(paste("observed",observeddata$varname,sep =" "))
  predicteddata$linelegend <- as.factor(paste("predicted", predicteddata$varname, sep = " "))
  
  compareplot <- ggplot(data=observeddata)
  compareplot <- compareplot + geom_point(aes(x=x,y=y,color=pointlegend))
  compareplot <- compareplot + geom_line(data = predicteddata, aes(x=x,y=y,color=linelegend))
  
  
  compareplot <- compareplot + xlab(xlabtitle)
  compareplot <- compareplot + ylab(ylabtitle)
  
  N <- length(unique(observeddata$varname))
  if(N==1){
    compareplot <- compareplot + scale_colour_manual(values= c("black","black"))
    shapelist <- c(16,NA)
    linelist <- c(0,1)
  } else {
    cbPalette <- c("#000000", "#E69F00", "#CC79A7", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")
    cbPalette <- cbPalette[1:N]
    compareplot <- compareplot + scale_colour_manual(values = rep(cbPalette,2))
    shapelist <- c(rep(16,N),rep(NA,N))
    linelist <- c(rep(0,N),rep(1,N))
  }
  compareplot <- compareplot + guides(colour = guide_legend(title=NULL,override.aes = list(shape=shapelist,linetype=linelist)),shape=guide_legend(title=NULL)) 
  return(compareplot)  
}
