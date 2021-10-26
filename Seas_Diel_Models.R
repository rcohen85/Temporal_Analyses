library(geepack)
library(statmod)
library(splines)
library(SimDesign)
library(gridExtra)

inDir = 'I:/TimeSeries'
seasDir = 'I:/Seasonal_Plots'
int = "Daily"

fileList = list.files(path=inDir,pattern=paste('*_',int,'.csv',sep=""),
                      full.names=TRUE,recursive=FALSE,
                      include.dirs=FALSE,no..=TRUE)


for (i in 1:numel(fileList)){  # for each species' file
  
  thisCT = data.frame(read.csv(fileList[i]))  # load file
  attribs = attributes(thisCT)  # find names of sites (all columns but first)
  sites = attribs$names[2:numel(attribs$names)]
  
  CTname = str_remove(fileList[i],paste(inDir,'/',sep="")) # get the species/CT name
  CTname = str_remove(CTname,paste('_',int,'.csv',sep=""))
  
  dateVec = as.POSIXlt(thisCT$Date,tz="GMT",format="%d-%b-%Y")
  
  # if it doesn't already exist, create directory to save figures
  if (!dir.exists(paste(seasDir,'/',CTname,sep=""))){
    dir.create(paste(seasDir,'/',CTname,sep=""))
  }
  
  for (j in 1:numel(sites)){ # for each site
    
    thisSite = as.numeric(thisCT[,j+1]) # get presence data for this site
    corr = acf(thisSite,lag.max=100,na.action=na.pass,plot=FALSE) 
    lagID = which(corr$acf<0.1) # determine lag at which autocorrelation is <0.15
    #itsVal = IntegralTimeScaleCalc(thisSite)
    numClust = length(thisSite)/(lagID[1]-1)
    clustID = rep(1:ceiling(numClust),each=lagID[1]) # create grouping vector for GEEGLM
    clustID = clustID[1:numel(thisSite)]
    
    noDat = which(is.na(thisSite)) # remove days with no presence data
    thisSite = thisSite[-noDat]
    thisSite = round(thisSite)
    reducedDateVec = dateVec[-noDat]
    reducedClustID = clustID[-noDat]
    Jday = reducedDateVec$yday
    monthGroup = month(reducedDateVec) # find which month each day of data falls in
    yearGroup = year(reducedDateVec) # find which year each day of data falls in
    
    pres = which(thisSite!=0)
    
    if (!numel(pres)==0){ # if there is any presence
      
      # Fit model
    modSm = geeglm(thisSite~bs(Jday,knots=mean(Jday))+bs(monthGroup,knots=mean(monthGroup))+as.factor(yearGroup),family=poisson(link="log"),id=reducedClustID,corstr="ar1")
    summary(modSm)
    # mod = geeglm(thisSite~Jday+monthGroup+as.factor(yearGroup),family=poisson,id=reducedClustID,corstr="ar1")
    # summary(mod)
    
    
    # Bootstrap parameter estimate confidence intervals
    JDayForPlotting<- seq(min(Jday), max(Jday), length=50)
    MonthForPlotting<- seq(min(monthGroup), max(monthGroup), length=50)
    BootstrapParameters<-rmvnorm(10000, coef(modSm), summary(modSm)$cov.unscaled)
    test<- glm(thisSite ~ bs(Jday,knots=mean(Jday))+bs(monthGroup,knots=mean(monthGroup))+as.factor(yearGroup),family=poisson)
    x1<-model.matrix(test)[,2:5]%*%coef(modSm)[c(2:5)]
    x2<-model.matrix(test)[,6:9]%*%coef(modSm)[c(6:9)]
    x3<-model.matrix(test)[,10:12]%*%coef(modSm)[c(1,10:12)]
    
    # Plot partial residuals
    # Julian Day
    BootstrapCoefs<- BootstrapParameters[,2:5]
    Basis<- bs(JDayForPlotting, knots=mean(Jday), Boundary.knots=range(Jday))
    RealFit<- Basis%*%coef(modSm)[c(2:5)]
    RealFitCenter1<- RealFit-mean(x1)-coef(modSm)[1]
    BootstrapFits<- Basis%*%t(BootstrapCoefs)
    quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
    cis<-apply(BootstrapFits, 1, quant.func)
    MinimumYlim1<- min(cis-mean(x1)-coef(modSm)[1])
    MaximumYlim1<- max(cis-mean(x1)-coef(modSm)[1])
    cil1<-cis[1,]-mean(x1)-coef(modSm)[1]
    ciu1<-cis[2,]-mean(x1)-coef(modSm)[1]
    plot1<-qplot(JDayForPlotting,RealFitCenter1,xlab="Julian Day", ylab="s(Julian Day)", 
                 main="a",ylim=c(MinimumYlim1, MaximumYlim1),
                 geom="line")+theme(axis.line = element_line(),
                                   panel.background=element_blank(),
                                   panel.grid.major=element_blank(),
                                   panel.grid.minor=element_blank())+geom_smooth(fill="grey",
                                                                               colour="black",
                                                                               aes(ymin=cil1,ymax=ciu1),
                                                                               stat="identity")+geom_rug(aes(x = Jday, y=-10000))
    # Month
    BootstrapCoefs<- BootstrapParameters[,6:9]
    Basis<- bs(MonthForPlotting, knots=mean(monthGroup), Boundary.knots=range(monthGroup))
    RealFit<- Basis%*%coef(modSm)[c(6:9)]
    RealFitCenter2<- RealFit-mean(x2)-coef(modSm)[1]
    BootstrapFits<- Basis%*%t(BootstrapCoefs)
    quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
    cis<-apply(BootstrapFits, 1, quant.func)
    MinimumYlim2<- min(cis-mean(x2)-coef(modSm)[1])
    MaximumYlim2<- max(cis-mean(x2)-coef(modSm)[1])
    cil2<-cis[1,]-mean(x2)-coef(modSm)[1]
    ciu2<-cis[2,]-mean(x2)-coef(modSm)[1]
    plot2<-qplot(MonthForPlotting,RealFitCenter2,
                 xlab="Month", ylab="s(Month)", 
                 main="b",ylim=c(MinimumYlim2, MaximumYlim2),
                 geom="line")+theme(axis.line = element_line(),
                                    panel.background=element_blank(),
                                    panel.grid.major=element_blank(),
                                    panel.grid.minor=element_blank())+geom_smooth(fill="grey",
                                                                                  colour="black",
                                                                                  aes(ymin=cil2,ymax=ciu2),
                                                                                  stat="identity")+geom_rug(aes(x = monthGroup, y=-10000))
   
  
    # Year
    BootstrapCoefs<- BootstrapParameters[,c(1,10:12)]
    #Basis<- bs(MonthForPlotting, knots=mean(monthGroup), Boundary.knots=range(monthGroup))
    facBasis = c(1:4)
    RealFit<- facBasis%*%coef(modSm)[c(1,10:12)]
    RealFitCenter2<- RealFit-mean(x1)-coef(modSm)[1]
    BootstrapFits<- facBasis%*%t(BootstrapCoefs)
    quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
    cis<-apply(BootstrapFits, 1, quant.func)
    MinimumYlim3<- min(cis-mean(x3)-coef(modSm)[1])
    MaximumYlim3<- max(cis-mean(x3)-coef(modSm)[1])
    cil3<-cis[1,]-mean(x3)-coef(modSm)[1]
    ciu3<-cis[2,]-mean(x3)-coef(modSm)[1]
    
    grid.arrange(plot1,plot2,ncol=1,nrow=2)
    
    
    }
  }
}
