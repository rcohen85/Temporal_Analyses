library(tidyverse)
library(geepack)
#library(statmod)
#library(splines)
library(splines2)
library(SimDesign)
#library(car)
library(pracma)
library(mgcv)
#library(forecast)
source("getPvalues.r")

inDir = 'I:/TimeSeries_ScaledByEffortError'
seasDir = 'I:/Seasonal_Plots'
int = "Hourly"

fileList = list.files(path=inDir,pattern=paste('*_',int,'.csv',sep=""),
                      full.names=TRUE,recursive=FALSE,
                      include.dirs=FALSE,no..=TRUE)
goodIdx = c(1,2,4,8,11,13:19)


for (i in goodIdx){  # for each species' file
  
  thisCT = data.frame(read.csv(fileList[i]))  # load file
  attribs = attributes(thisCT)  # find names of sites (all columns but first)
  sites = attribs$names[-1]
  
  CTname = str_remove(fileList[i],paste(inDir,'/',sep="")) # get the species/CT name
  CTname = str_remove(CTname,paste('_',int,'.csv',sep=""))
  
  #dateVec = as.POSIXlt(thisCT$Date,tz="GMT",format="%d-%b-%Y")
  if (int=="Daily") {
    dateVec = as.Date(thisCT$Date,format="%d-%b-%Y",tz="GMT")
  } else if (int=="Hourly") {
    dateVec = as.Date(thisCT$Hour,format="%d-%b-%Y %H:%M:%S",tz="GMT")
  }
  
  
  # if it doesn't already exist, create directory to save figures
  if (!dir.exists(paste(seasDir,'/',CTname,sep=""))){
    dir.create(paste(seasDir,'/',CTname,sep=""))
  }
  
  for (j in 1:numel(sites)){ # for each site
    
    thisSite = as.numeric(thisCT[,j+1]) # get presence data for this site
    pres = which(thisSite!=0)
    
    if (!numel(pres)==0){ # if there is any presence
      
      # plot presence to determine if it's worth modeling
      plot(dateVec,thisSite,main=paste(CTname,'at',sites[j],sep=" "),
           ylab="# 5-min Bins w Presence")
      
      userVote = c()
      while (numel(userVote)==0){
        
        userVote = readline(prompt="Enter 1 to proceed with modeling, or 0 to skip to next deployment: ")
        
        if (!numel(userVote)==0 && userVote!=1 && userVote!=0){
          message('WARNING: Entry not allowed')
          userVote = c()
          
        }
        else if (numel(userVote)==0){
          message('WARNING: Entry not allowed')
        }
        
      }
      
      if (as.numeric(userVote)==1){ # if user says "yes" to modeling
       
        # Prepare data
        corr = acf(thisSite,lag.max=1500,na.action=na.pass,plot=FALSE) 
        lagID = which(abs(corr$acf)<0.2) # determine lag at which autocorrelation is <0.2
        # itsVal = IntegralTimeScaleCalc(thisSite)
        numClust = length(thisSite)/(lagID[1]-1)
        if (numClust<length(thisSite)){
          clustID = rep(1:ceiling(numClust),each=lagID[1]) # create grouping vector for GEEGLM
          clustID = clustID[1:numel(thisSite)]
        } else {
          clustID = 1:length(thisSite)
        }
        
        
        noDat = which(is.na(thisSite)) # remove days/hours with no effort
        if (!numel(noDat) == 0) {
          thisSite = thisSite[-noDat]
          reducedDateVec = dateVec[-noDat]
          reducedClustID = clustID[-noDat]
          Jday = as.numeric(format(reducedDateVec, "%j"))
          # monthGroup = month(reducedDateVec) # find which month each day of data falls in
          yearGroup = year(reducedDateVec) # find which year each day of data falls in
        } else if (numel(noDat) == 0) {
          reducedDateVec = dateVec
          reducedClustID = clustID
          Jday = as.numeric(format(reducedDateVec, "%j"))
          # monthGroup = month(reducedDateVec) # find which month each day of data falls in
          yearGroup = year(reducedDateVec) # find which year each day of data falls in
        }
        
        if (j == 7) {
          # create variable coding for change of site at HAT
          hatSite = rep(2, length(reducedDateVec))
          hatAdates = which(reducedDateVec <= as.Date("06-Feb-2017", format =
                                                        "%d-%b-%Y"))
          hatSite[hatAdates] = 1
        }
        
        # account for leap day in 2016, shift Julian days by value of 1 
        # (so we don't try modeling Julian day 366 based on 1 data point)
        leapIdx = which(yearGroup==2016)
        Jday[leapIdx] = Jday[leapIdx]-1
        
        # round presence data back to integers so it's Poisson distributed again
        thisSite = round(thisSite)
        
        
        # Fit GAM
        if (j==7){
          JDayGAM = gam(thisSite~s(Jday,bs="cc",k=6)+as.factor(yearGroup)+as.factor(hatSite),
                        family=tw)
        } else {JDayGAM = gam(thisSite~s(Jday,bs="cc",k=6)+as.factor(yearGroup),
                      family=tw)}
        sinkName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GAMSummary.txt",sep="")
        sink(sinkName)
        print(summary(JDayGAM))
        sink()
        
        # Plot GAM partial residuals
        saveName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GAM.png",sep="")
        if (j==7){
          png(saveName,width=800,height=500)
        } else {
          png(saveName,width=800,height=400)
        }
        
        plot.gam(JDayGAM,all.terms=TRUE,pages=1,main=paste(CTname,'at',sites[j]))
        while (dev.cur()>1) {dev.off()}
        
        
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
        JDayForPlotting<- seq(min(Jday), max(Jday), length=512)
        Basis<- mSpline(JDayForPlotting,  # spline spanning range of X values
                        knots=quantile(Jday, probs=c(0.333,0.666)),
                        Boundary.knots=c(1,365),
                        periodic=T) # basis functions for smooth function
        RealFit<- Basis%*%coef(modJday)[c(2:3)] # multiply basis functions by model coefficients to get values of spline at each X
        RealFitCenterJ<- RealFit-mean(Jx1)-coef(modJday)[1] # adjust offset
        JDayBootstrapFits<- Basis%*%t(JDayBootstrapCoefs) # get spread of spline values at each X based on distributions of each coefficient
        quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
        cisJ<-apply(JDayBootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
        Jcil<-cisJ[1,]-mean(Jx1)-coef(modJday)[1] # upper CI bound
        Jciu<-cisJ[2,]-mean(Jx1)-coef(modJday)[1] # lower CI bound
        #MinimumYlimJ<- min(cisJ-mean(Jx1)-coef(modJday)[1])
        #MaximumYlimJ<- max(cisJ-mean(Jx1)-coef(modJday)[1])
        saveName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDayPlot.png",sep="")
        dJday = stats::density(Jday,na.rm = TRUE) # Calculate kernel density of Jday observations
        plotDF = data.frame(JDayForPlotting,RealFitCenterJ,dJday$x,dJday$y)
        colnames(plotDF) = c("Jday","Fit","Day","Density")
        plotDF$Density = plotDF$Density/max(plotDF$Density) # normalize kernel density
        plotDF$Density = plotDF$Density
        # png(saveName,width=600,height=350)
        #ggplot(dens,aes(x=Day,y=Density)) + geom_area(data=plotDF,x=Day,y=Density,fill=4,alpha=0.2)+
        ggplot(plotDF, aes(Jday, Fit),
               ) + geom_smooth(fill = "grey",
                               colour = "black",
                               aes(ymin=Jcil, ymax=Jciu),
                               stat ="identity"
                ) + labs(x = "Julian Day",
                         y = "s(Julian Day)",
                         title = paste(CTname, 'at',sites[j]),
                 ) + theme(axis.line = element_line(),
                           panel.background = element_blank()
                 ) + geom_polygon(data=plotDF,
                                  aes(Day,Density),
                                  fill=4,
                                  alpha=0.2)
          
          
          #geom_area(data = plotDF,
                              # aes(Day, Density),
                              # fill = 4,
                               #alpha = 0.2)
          #+ geom_density(data=plotDF,aes(Jday,Fit,fill=4,alpha=0.2),stat="density")

          
        # qplot(
        #   JDayForPlotting,
        #   RealFitCenterJ,
        #   xlab = "Julian Day",
        #   ylab = "s(Julian Day)",
        #   main = paste(CTname, 'at', sites[j]),
        #   ylim = c(MinimumYlimJ, MaximumYlimJ),
        #   geom = "line"
        # ) + theme(
        #   axis.line = element_line(),
        #   panel.background = element_blank(),
        #   panel.grid.major = element_blank(),
        #   panel.grid.minor = element_blank()
        # ) + geom_smooth(
        #   fill = "grey",
        #   colour =
        #     "black",
        #   aes(ymin =
        #         Jcil, ymax = Jciu),
        #   stat =
        #     "identity") + geom_rug(aes(x = Jday, y = -10000))
        polygon(Jday,thisSite, xlim=range(Jday),col="lightgray",border = FALSE, main = "",xlab =  "",ylab = "")
        ggsave(saveName, device = "png")
        while (dev.cur()>1) {dev.off()}
        
       
        # Year (as boxplot)
        # Center intercept (1st level of year factor) at 0 and show other levels relative to it
        AdjustedYearCoefs = data.frame(YearBootstrapCoefs[,1]-mean(YearBootstrapCoefs[,1]),
                                       YearBootstrapCoefs[,2],
                                       YearBootstrapCoefs[,3],
                                       YearBootstrapCoefs[,4])
        colnames(AdjustedYearCoefs) = c("2016","2017","2018","2019")
        # AdjustedYearCoefs = apply(AdjustedYearCoefs,2,exp) return to units of response var?
        
        
        saveName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_YearPlot.png",sep="")
        png(saveName,width=500,height=400)
        boxplot(AdjustedYearCoefs,main=paste(CTname,'at',sites[j]),outline=FALSE,
                ylab=c("Year Coefficients"))
        #grid.arrange(plot1,plot2,ncol=1,nrow=2)
        while (dev.cur()>1) {dev.off()}
        
        # Site (as boxplot), if HAT
        if (j==7){
          AdjustedHATCoefs = data.frame(hatSiteBootstrapCoefs[,1]-mean(hatSiteBootstrapCoefs[,1]),
                                        hatSiteBootstrapCoefs[,2])
          colnames(AdjustedHATCoefs) = c("HAT_A","HAT_B")
          
          saveName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_HATSitePlot.png",sep="")
          png(saveName,width=500,height=400)
          boxplot(AdjustedHATCoefs,main=paste(CTname,'at',sites[j]),outline=FALSE,
                  ylab=c("HAT Site Coefficients"))
          while (dev.cur()>1) {dev.off()}
          
        }
      }
      
    }
    
  }
}
