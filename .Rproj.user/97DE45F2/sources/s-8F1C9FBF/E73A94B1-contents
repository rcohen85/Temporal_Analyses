library(tidyverse)
library(lubridate)
library(geepack)
library(splines2)
library(SimDesign)
library(pracma)
library(mgcv)
source("getPvalues.r")
#library(statmod)
#library(splines)
#library(forecast)
#library(car)


modelDFDir = 'I:/TimeSeries_ScaledByEffortError'
dfList = list.files(path=modelDFDir,pattern=paste('*Master.csv',sep=""),
                                          full.names=TRUE,recursive=FALSE,
                                          include.dirs=FALSE,no..=TRUE)
        
        
for (i in 1:numel(dfList)){
  
  # load file
  thisSite = data.frame(read.csv(dfList[i]))
  
  # # Fit GAM
  # if (j==7){
  #   JDayGAM = gam(thisSite~s(Jday,bs="cc",k=6)+as.factor(yearGroup)+as.factor(hatSite),
  #                 family=tw)
  # } else {JDayGAM = gam(thisSite~s(Jday,bs="cc",k=6)+as.factor(yearGroup),
  #               family=tw)}
  # sinkName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GAMSummary.txt",sep="")
  # sink(sinkName)
  # print(summary(JDayGAM))
  # sink()
  # 
  # # Plot GAM partial residuals
  # saveName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GAM.png",sep="")
  # if (j==7){
  #   png(saveName,width=800,height=500)
  # } else {
  #   png(saveName,width=800,height=400)
  # }
  # 
  # plot.gam(JDayGAM,all.terms=TRUE,pages=1,main=paste(CTname,'at',sites[j]))
  # while (dev.cur()>1) {dev.off()}
  
  # Fit GEEGLM; previously used bs(Jday,knots=mean(Jday))
  if (j==7){
    modJday = geeglm(thisSite~mSpline(Jday,
                                      knots=quantile(Jday, probs=c(0.333,0.666)),
                                      Boundary.knots=c(1,365),
                                      periodic=T)
                     +as.factor(yearGroup)+as.factor(hatSite),
                     family=poisson(link="log"),
                     id=reducedClustID,
                     corstr="ar1")
  } else {
    modJday = geeglm(thisSite~mSpline(Jday,
                                      knots=quantile(Jday, probs=c(0.333,0.666)),
                                      Boundary.knots=c(1,365),
                                      periodic=T)
                     +as.factor(yearGroup),
                     family=poisson(link="log"),
                     id=reducedClustID,
                     corstr="ar1")
  }
  
  sinkName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLMSummary.txt",sep="")
  sink(sinkName)
  #print(summary(modJday))
  getPvalues(modJday)
  sink()
  
  
  # Bootstrap GEEGLM parameter estimates for later construction of confidence intervals
  BootstrapParameters<-rmvnorm(10000, coef(modJday), summary(modJday)$cov.unscaled)
  JDayBootstrapCoefs<- BootstrapParameters[,2:3]
  YearBootstrapCoefs<- BootstrapParameters[,c(1,4:6)]
  if (j==7){ # if HAT, get coefficients for Site
    hatSiteBootstrapCoefs<- BootstrapParameters[,c(1,7)]
    quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
    cisH<-apply(hatSiteBootstrapCoefs, 2, quant.func)
  }
  
  # Predict presence at each X value using model coefficients
  if (j==7){ # if HAT, include Site as factor
    Jx1<-model.matrix(modJday)[,2:3]%*%coef(modJday)[c(2:3)]
    Jx2<-model.matrix(modJday)[,c(1,4:6)]%*%coef(modJday)[c(1,4:6)]
    Jx3<-model.matrix(modJday)[,c(1,7)]%*%coef(modJday)[c(1,7)]
  } else {
    Jx1 = model.matrix(modJday)[,2:3]%*%coef(modJday)[c(2:3)]
    Jx2 = model.matrix(modJday)[,c(1,4:6)]%*%coef(modJday)[c(1,4:6)]
  }
  
  
  # Plot GEEGLM partial residuals
  # Julian Day
  JDayForPlotting<- seq(min(Jday), max(Jday), length=256)
  Basis<- mSpline(JDayForPlotting,  # spline spanning range of X values
                  knots=quantile(Jday, probs=c(0.333,0.666)),
                  Boundary.knots=c(1,365),
                  periodic=T) # basis functions for smooth function
  RealFit<- Basis%*%coef(modJday)[c(2:3)] # multiply basis functions by model coefficients to get values of spline at each X
  RealFitCenterJ<- RealFit-mean(Jx1)-coef(modJday)[1] # adjust offset
  JDayBootstrapFits<- Basis%*%t(JDayBootstrapCoefs) # get spread of spline values at each X based on distributions of each coefficient
  quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
  cisJ<-apply(JDayBootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
  Jcil<-cisJ[1,]-mean(Jx1)-coef(modJday)[1] # lowerCI bound
  Jciu<-cisJ[2,]-mean(Jx1)-coef(modJday)[1] # upper CI bound
  
  
  plotDF = data.frame(JDayForPlotting,RealFitCenterJ)
  colnames(plotDF) = c("Jday","Fit")
  dJday = stats::density(Jday,na.rm = TRUE,n=256,from=1,to=365) # Calculate kernel density of Jday observations
  dens = data.frame(c(dJday$x,rev(dJday$x)),c(dJday$y,rep(0,length(dJday$y))))
  colnames(dens) = c("Day","Density")
  dens$Density = dens$Density/max(dens$Density) # normalize kernel density
  if (min(Jcil)<0){ # set kernel density at bottom of y axis
    dens$Density = dens$Density - abs(min(Jcil)) 
  } else {
    dens$Density = dens$Density + min(Jcil)
  }
  
  #saveName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDayPlot.png",sep="")
  JD = ggplot(plotDF, aes(Jday, Fit),
  ) + geom_polygon(data=dens,
                   aes(Day,Density),
                   fill=4,
                   alpha=0.2
  ) + geom_smooth(fill = "grey",
                  colour = "black",
                  aes(ymin=Jcil, ymax=Jciu),
                  stat ="identity"
  ) + labs(x = "Julian Day",
           y = "s(Julian Day)",
           title = paste(CTname, 'at',sites[j]),
  ) + theme(axis.line = element_line(),
            panel.background = element_blank()
  )
  
  ggsave(saveName,device="png")
  while (dev.cur()>1) {dev.off()}
  
  
  # Year (as boxplot)
  # Center intercept (1st level of year factor) at 0 and show other levels relative to it
  AdjustedYearCoefs = data.frame(c(YearBootstrapCoefs[,1]-mean(YearBootstrapCoefs[,1]),
                                   YearBootstrapCoefs[,2],
                                   YearBootstrapCoefs[,3],
                                   YearBootstrapCoefs[,4]),
                                 as.factor(rep(2016:2019,each=10000)))
  colnames(AdjustedYearCoefs) = c("Coefficient","Year")
  # AdjustedYearCoefs = apply(AdjustedYearCoefs,2,exp) return to units of response var?
  
  #saveName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_YearPlot.png",sep="")
  Yr = ggplot(AdjustedYearCoefs,aes(Year,Coefficient)
  ) + geom_boxplot(
  )+ theme(axis.line = element_line(),
           panel.background = element_blank()
  ) + labs(title = paste(CTname, 'at',sites[j]))
  ggsave(saveName,device="png")
  while (dev.cur()>1) {dev.off()}
  
  # Site (as boxplot), if HAT
  if (j==7){
    AdjustedHATCoefs = data.frame(c(hatSiteBootstrapCoefs[,1]-mean(hatSiteBootstrapCoefs[,1]),
                                    hatSiteBootstrapCoefs[,2]),
                                  as.factor(rep(c("HAT_A","HAT_B"),each=10000)))
    colnames(AdjustedHATCoefs) = c("Coefficient","Site")
    
    saveName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_HATSitePlot.png",sep="")
    ggplot(AdjustedHATCoefs,aes(Site,Coefficient)
    ) + geom_boxplot(
    )+ theme(axis.line = element_line(),
             panel.background = element_blank()
    )+labs(title = paste(CTname, 'at',sites[j]))
    ggsave(saveName,device="png")
    while (dev.cur()>1) {dev.off()}
    
  }
}
