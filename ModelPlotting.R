#### PLOT ALL FIXED SPLINE MODELS FOR A GIVEN SPECIES --------------------------------------
library(tidyverse)
library(broom)
library(lubridate)
library(splines2)
library(pracma)
library(SimDesign)
library(boot)
library(gridExtra)
library(ResourceSelection)
library(performance)


modelDFDir = 'J:/Chpt_2/TimeSeries_ScaledByEffortError'
outDir = 'J:/Chpt_2/ModelOutput'
int = "5minBin"

dfList = list.files(path=modelDFDir,pattern=paste('*',int,'_Master.csv',sep=""),
                    full.names=TRUE,recursive=FALSE,
                    include.dirs=FALSE,no..=TRUE)
lunList = list.files(path=modelDFDir,pattern=paste('*',int,'_MasterLun.csv',sep=""),
                     full.names=TRUE,recursive=FALSE,
                     include.dirs=FALSE,no..=TRUE)

species = list.dirs(outDir,recursive=FALSE)

for ( i in 1:numel(species)){
  
  JDmodFiles = list.files(path=species[i],pattern="*Model_JDYrNormTime.Rdata",
                          full.names=TRUE,recursive=FALSE,include.dirs=FALSE,no..=TRUE)
  LunmodFiles = list.files(path=species[i],pattern="*Model_LunPhaseAltPres.Rdata",
                           full.names=TRUE,recursive=FALSE,include.dirs=FALSE,no..=TRUE)
  CTname = str_remove(species[i],paste(outDir,'/',sep=""))
  
  sites = list()
  for (j in 1:numel(JDmodFiles)){
    
    site = str_remove(JDmodFiles[j],paste(outDir,"/",CTname,"/",sep=""))
    sites = c(sites,str_remove(site,"_5minBin_Model_JDYrNormTime.Rdata"))
    if (sites[j]=="NULL"){
      stop("Didn't get site name")}

    load(JDmodFiles[j]) # load JD + NormTime + Yr model

    if (tempMod$geese$error==0){ # not all models converged, don't plot non-converged models

      # find associated master dataframe
      thisSpec = which(str_detect(dfList,CTname))
      atSite = which(str_detect(dfList,unlist(sites[j])))
      thisModInd = intersect(thisSpec,atSite)
      thisSite = data.frame(read.csv(dfList[thisModInd])) # JD model data frame

      # Indices of coefficients for each covar
      JDInd = numeric()
      NTInd = numeric()
      YrInd = numeric()
      if (tempMod$formula =='Pres ~ JDs + NTs + YrF'){
        JDInd = c(2:4)
        NTInd = c(5:7)
        YrInd = c(1,8:9)
      } else if (tempMod$formula =='Pres ~ JDs + NTs'){
        JDInd = c(2:4)
        NTInd = c(5:7)
      } else if (tempMod$formula =='Pres ~ JDs + YrF'){
        JDInd = c(2:4)
        YrInd = c(1,5:6)
      } else if (tempMod$formula =='Pres ~ NTs + YrF'){
        NTInd = c(2:4)
        YrInd = c(1,5:6)
      } else if (tempMod$formula =='Pres ~ JDs'){
        JDInd = c(2:4)
      } else if (tempMod$formula =='Pres ~ NTs'){
        NTInd = c(2:4)
      } else if (tempMod$formula =='Pres ~ YrF'){
        YrInd = c(1:3)
      }

      if (sites[j]=="HAT"){
        startInd = which(thisSite$TimeStamp>=as.POSIXct('2017-05-01 00:00:00',format="%Y-%m-%d %H:%M:%S",tz="GMT"))
        thisSite = thisSite[startInd,]
        if (!isempty(YrInd)){
          YrInd = c(1,length(coef(tempMod))) # only two years of data at HAT
        }
      }

      ## Create model fit plots ---------------------
      modFit = data.frame(Date=as.Date(thisSite$TimeStamp),
                          Pres=thisSite$Presence,
                          Fits=tempMod$fitted.values,
                          Res=tempMod$residuals)

      PresPlot = ggplot(modFit
      ) + geom_point(aes(x=Date,y=Pres),
                     color="#000000",
                     size=2
      )+scale_x_continuous(breaks=c(as.Date("2016-05-01"),
                                    as.Date("2016-11-01"),
                                    as.Date("2017-05-01"),
                                    as.Date("2017-11-01"),
                                    as.Date("2018-05-01"),
                                    as.Date("2018-11-01"),
                                    as.Date("2019-04-30"))
      ) + labs(x="",y="Presence")
      FitPlot = ggplot(modFit
      )+geom_point(aes(x=Date,y=Fits),
                   color="#1976D2",
                   size=2
      )+scale_x_continuous(breaks=c(as.Date("2016-05-01"),
                                    as.Date("2016-11-01"),
                                    as.Date("2017-05-01"),
                                    as.Date("2017-11-01"),
                                    as.Date("2018-05-01"),
                                    as.Date("2018-11-01"),
                                    as.Date("2019-04-30"))
      ) + labs(x="",y="Fitted Values")
      ModRes = ggplot(modFit
      ) + geom_point(aes(x=Date,y=Res),
                     size=2
      ) +scale_x_continuous(breaks=c(as.Date("2016-05-01"),
                                     as.Date("2016-11-01"),
                                     as.Date("2017-05-01"),
                                     as.Date("2017-11-01"),
                                     as.Date("2018-05-01"),
                                     as.Date("2018-11-01"),
                                     as.Date("2019-04-30"))
      ) + labs(x="",y="Residuals")

      HL = hoslem.test(thisSite$Presence,fitted(tempMod))
      HLPlot = data.frame(Expected=HL$expected[,1],Observed=HL$observed[,1])
      HLPlot$X = seq(min(HL$observed[,1]),max(HL$observed[,1]),length.out=length(HLPlot$Expected))
      HLPlot$Y = seq(min(HL$observed[,1]),max(HL$observed[,1]),length.out=length(HLPlot$Expected))

      ObsPred = ggplot(HLPlot
      ) + geom_line(aes(x=X,y=Y)
      ) + geom_point(aes(y=Expected,x=Observed)
      ) + labs(x="Observed",y="Expected"
      ) + annotate("text",
                   label=paste('p-value: ',as.character(HL$p.value),sep=""),
                   size=4,
                   x=1.005*min(HL$observed[,1]),
                   y=0.999*max(HL$expected[,1]))

      png(file=paste(outDir,'/',CTname,'/',sites[j],"_ModelFit.png",sep=""),width = 1400, height = 800, units = "px")
      grid.arrange(PresPlot,FitPlot,ModRes,ObsPred,ncol=1,nrow=4,top=paste(CTname,'at',sites[j]))
      while (dev.cur()>1) {dev.off()}
      
      png(file=paste(outDir,'/',CTname,'/',sites[j],"_BinResid.png",sep=""),width = 400, height = 300, units = "px")
      binned_residuals_NP(tempMod)
      while (dev.cur()>1) {dev.off()}
      
      ## Bootstrap GEEGLM parameter estimates for later construction of confidence intervals ----------------
      BootstrapParameters1<-rmvnorm(10000, coef(tempMod), summary(tempMod)$cov.unscaled)
      JDayBootstrapCoefs<- BootstrapParameters1[,JDInd]
      # DayPhaseBootstrapCoefs<- BootstrapParameters1[,DPInd]
      NTBootstrapCoefs<- BootstrapParameters1[,NTInd]
      YearBootstrapCoefs<- BootstrapParameters1[,YrInd]

      # Predict presence at each X value using model coefficients
      JxJD = model.matrix(tempMod)[,JDInd]%*%coef(tempMod)[JDInd]
      JxNT = model.matrix(tempMod)[,NTInd]%*%coef(tempMod)[NTInd]
      # JxDP = model.matrix(tempMod)[,DPInd]%*%coef(tempMod)[DPInd]
      JxYr = model.matrix(tempMod)[,YrInd]%*%coef(tempMod)[YrInd]

      ### Generate GEEGLM partial residual plots ---------------------------
      # Julian Day ---------------------------
      if (!isempty(JDInd)){
        JDayForPlotting<- seq(min(thisSite$JulianDay), max(thisSite$JulianDay), length=5000)
        JBasis<- mSpline(JDayForPlotting,  # spline spanning range of X values
                         knots=quantile(thisSite$JulianDay, probs=c(0.275,0.5,0.725)),
                         Boundary.knots=c(1,365),
                         periodic=T) # basis functions for smooth function
        RealFitJ<- JBasis%*%(coef(tempMod)[JDInd]) # multiply basis functions by model coefficients to get values of spline at each X
        RealFitCenterJ<- RealFitJ+coef(tempMod)[1]-mean(JxJD) # adjust offset
        RealFitCenterJ<- inv.logit(RealFitCenterJ)
        JDayBootstrapFits<- (JBasis%*%t(JDayBootstrapCoefs))+coef(tempMod)[1] # get spread of spline values at each X based on distributions of each coefficient
        quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
        cisJ<-apply(JDayBootstrapFits, 1, quant.func)-mean(JxJD) # confidence interval of smooth function estimate
        Jcil<-inv.logit(cisJ[1,]) # lowerCI bound
        Jciu<-inv.logit(cisJ[2,]) # upper CI bound

        JplotDF = data.frame(JDayForPlotting,RealFitCenterJ)
        colnames(JplotDF) = c("Jday","Fit")

        JD = ggplot(JplotDF, aes(Jday, Fit),
        ) + geom_smooth(fill = "grey",
                        colour = "black",
                        aes(ymin=Jcil, ymax=1.05*Jciu),
                        stat ="identity"
        ) + labs(x = "Julian Day",
                 y = "Probability",
        ) + scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),
                               label=c("J","F","M","A","M","J","J","A","S","O","N","D")
        ) + theme(axis.line = element_line(size=0.2),
                  panel.background = element_blank()
        )

        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDayPlot_Probs.png",sep="")
        ggsave(saveName,device="png", width=2, scale=3, height=0.5, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDayPlot_Probs.pdf",sep="")
        ggsave(saveName,device="pdf", width=2, scale=3, height=0.5, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}

        JDdens = ggplot(thisSite,aes(x=JulianDay)
        )+geom_histogram(aes(y=..ncount..),
                         fill='#66B2FF',
                         binwidth = 0.5,
                         alpha=0.5
        ) + scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),
                               label=c("J","F","M","A","M","J","J","A","S","O","N","D")
        ) + labs(y="Count",x=NULL
        ) + theme_minimal()
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDDataDensity.png",sep="")
        ggsave(saveName,device="png", width=2, scale=4, height=0.5, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDDataDensity.pdf",sep="")
        ggsave(saveName,device="pdf", width=2, scale=4, height=0.5, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
      }
      # Norm Time ---------------------------
      if (!isempty(NTInd)){
        NTForPlotting<- seq(min(thisSite$NormTime), max(thisSite$NormTime), length=5000)
        NTBasis<- mSpline(NTForPlotting,  # spline spanning range of X values
                          knots=quantile(thisSite$NormTime, probs=c(0.275,0.5,0.725)),
                          Boundary.knots=c(-1,1),
                          periodic=T) # basis functions for smooth function
        RealFitNT<- NTBasis%*%(coef(tempMod)[NTInd]) # multiply basis functions by model coefficients to get values of spline at each X
        RealFitCenterNT<- RealFitNT+coef(tempMod)[1]-mean(JxNT) # adjust offset
        RealFitCenterNT<- inv.logit(RealFitCenterNT)
        NTBootstrapFits<- (NTBasis%*%t(NTBootstrapCoefs))+coef(tempMod)[1] # get spread of spline values at each X based on distributions of each coefficient
        quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
        cisNT<-apply(NTBootstrapFits, 1, quant.func)-mean(JxNT) # confidence interval of smooth function estimate
        NTcil<-inv.logit(cisNT[1,]) # lowerCI bound
        NTciu<-inv.logit(cisNT[2,]) # upper CI bound

        NTplotDF = data.frame(NTForPlotting,RealFitCenterNT)
        colnames(NTplotDF) = c("NormTime","Fit")

        NT = ggplot(NTplotDF, aes(NormTime, Fit),
        ) + geom_smooth(fill = "grey",
                        colour = "black",
                        aes(ymin=NTcil, ymax=1.05*NTciu),
                        stat ="identity"
        ) + labs(x = NULL,
                 y = "Probability"
        ) + coord_cartesian(xlim=c(-1,1)
        ) + scale_x_continuous(breaks=c(-1,0,1),
                               labels=c("Sunrise","Sunset","Sunrise")
        ) + theme(axis.line = element_line(size=0.2),
                  panel.background = element_blank()
        )

        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NormTimePlot_Probs.png",sep="")
        ggsave(saveName,device="png", width=2, scale=3, height=0.5, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NormTimePlot_Probs.pdf",sep="")
        ggsave(saveName,device="pdf", width=2, scale=3, height=0.5, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}

        NTdens = ggplot(thisSite,aes(x=NormTime)
        )+geom_histogram(aes(y=..ncount..),
                         fill='#66B2FF',
                         binwidth = 0.05,
                         alpha=0.5
        ) + scale_x_continuous(breaks=c(-1,0,1),
                               label=c("Sunrise","Sunset","Sunrise")
        ) + coord_cartesian(xlim=c(-1,1)
        ) + labs(y="Count",x=NULL
        ) + theme_minimal()
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NTDataDensity.png",sep="")
        ggsave(saveName,device="png", width=1, scale=4, height=0.5, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NTDataDensity.pdf",sep="")
        ggsave(saveName,device="pdf", width=1, scale=4, height=0.5, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
      }
      # Year (as boxplot) -----------------------------------
      if (!isempty(YrInd)){
        # Center intercept (1st level of year factor) at 0 and show other levels relative to it
        if (sites[j]=="HAT"){ # only 2 years of data from HATB
          AdjustedYearCoefs = data.frame(c(YearBootstrapCoefs[,1],
                                           YearBootstrapCoefs[,2]+mean(YearBootstrapCoefs[,1])),
                                         as.factor(rep(1:2,each=10000)))
        } else {
          AdjustedYearCoefs = data.frame(c(YearBootstrapCoefs[,1],
                                           YearBootstrapCoefs[,2]+mean(YearBootstrapCoefs[,1]),
                                           YearBootstrapCoefs[,3]+mean(YearBootstrapCoefs[,1])),
                                         as.factor(rep(1:3,each=10000)))}
        colnames(AdjustedYearCoefs) = c("Coefficient","Year")
        AdjustedYearCoefs$Prob = inv.logit(AdjustedYearCoefs$Coefficient)
        
        # calculate quantiles for plotting limits
        quants = AdjustedYearCoefs %>%
          group_by(Year) %>%
          summarize(q25 = quantile(Prob,probs=0.25),
                    q75 = quantile(Prob,probs=0.75))
        iqr = quants$q75-quants$q25
        
        Yr = ggplot(AdjustedYearCoefs,aes(Year,Prob)
        ) + geom_boxplot(varwidth=TRUE,
                         outlier.shape=NA,
                         lwd=0.2
        ) + labs(x='Study Year',y="Probability"
        ) + coord_cartesian(ylim = c(0,(1.5*max(iqr))+max(quants$q75))
        ) + theme(axis.line = element_line(size=0.2),
                  panel.background = element_blank())

        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_YearPlot_Probs.png",sep="")
        ggsave(saveName,device="png", width=2, scale=3, height=0.75, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_YearPlot_Probs.pdf",sep="")
        ggsave(saveName,device="pdf", width=2, scale=3, height=0.75, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}

        Yrdens = ggplot(thisSite,aes(x=StudyYear)
        )+geom_histogram(aes(y=..ncount..),
                         fill='#66B2FF',
                         binwidth = 0.5,
                         alpha=0.5
        ) + scale_x_continuous(breaks=c(1,2,3)
        ) + labs(y="Count",x="Study Year"
        ) + theme_minimal()
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_YrDataDensity.png",sep="")
        ggsave(saveName,device="png", width=1, scale=4, height=0.5, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_YrDataDensity.pdf",sep="")
        ggsave(saveName,device="pdf", width=1, scale=4, height=0.5, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}

      }
    }
  }
  
  sites=list()
  for (k in 1:numel(LunmodFiles)){

  #   site = str_remove(LunmodFiles[k],paste(outDir,"/",CTname,"/",sep=""))
  #   sites = c(sites,str_remove(site,"_5minBin_Model_LunPhaseAltPres.Rdata"))
  #   if (sites[k]=="NULL"){
  #     stop("Didn't get site name")}
  # 
  #   load(LunmodFiles[k]) # load LunPhase + Alt + Pres model
  # 
  #   if (lunMod$geese$error==0){
  # 
  #     # find associated master dataframe
  #     thisSpec = which(str_detect(lunList,CTname))
  #     atSite = which(str_detect(lunList,unlist(sites[k])))
  #     thisModInd = intersect(thisSpec,atSite)
  #     masterLun = data.frame(read.csv(lunList[thisModInd])) # lunar model data frame
  # 
  #     masterLun$MoonAltitude = masterLun$MoonAltitude*(180/pi) # convert altitude from radians to degrees
  # 
  #     if (sites[k]=="HAT"){
  #       startInd = which(masterLun$NightBinTimes>=as.POSIXct('2017-05-01 00:00:00',format="%Y-%m-%d %H:%M:%S",tz="GMT"))
  #       masterLun = masterLun[startInd,]
  #     }
  # 
  #     # Indices of coefficients for each covar
  #     MPInd = numeric()
  #     AltInd = numeric()
  #     PresInd = numeric()
  #     if (lunMod$formula =='NightPres ~ MPhs + MAs + MPrF'){
  #       MPInd = c(2:4)
  #       AltInd = c(5:9)
  #       PresInd = c(1,10:11)
  #     } else if (lunMod$formula =='NightPres ~ MPhs + MAs'){
  #       MPInd = c(2:4)
  #       AltInd = c(5:9)
  #     } else if (lunMod$formula =='NightPres ~ MPhs + MprF'){
  #       MPInd = c(2:4)
  #       PresInd = c(1,5:6)
  #     } else if (lunMod$formula =='NightPres ~ MAs + MprF'){
  #       AltInd = c(2,3)
  #       PresInd = c(1,4:5)
  #     } else if (lunMod$formula =='NightPres ~ MPhs'){
  #       MPInd = c(2:4)
  #     } else if (lunMod$formula =='NightPres ~ MAs'){
  #       AltInd = c(2:6)
  #     } else if (lunMod$formula =='NightPres ~ MPrF'){
  #       PresInd = c(1:3)
  #     }
  # 
  #     ## Create model fit plots --------------------------
  #     modFit = data.frame(Date=as.Date(masterLun$NightBinTimes),
  #                         Pres=masterLun$NightPres,
  #                         Fits=lunMod$fitted.values,
  #                         Res=lunMod$residuals)
  # 
  #     PresPlot = ggplot(modFit
  #     ) + geom_point(aes(x=Date,y=Pres),
  #                    color="#000000",
  #                    size=2
  #     )+scale_x_continuous(breaks=c(as.Date("2016-05-01"),
  #                                   as.Date("2016-11-01"),
  #                                   as.Date("2017-05-01"),
  #                                   as.Date("2017-11-01"),
  #                                   as.Date("2018-05-01"),
  #                                   as.Date("2018-11-01"),
  #                                   as.Date("2019-04-30"))
  #     ) + labs(x="",y="Presence")
  #     FitPlot = ggplot(modFit
  #     )+geom_point(aes(x=Date,y=Fits),
  #                  color="#1976D2",
  #                  size=2
  #     )+scale_x_continuous(breaks=c(as.Date("2016-05-01"),
  #                                   as.Date("2016-11-01"),
  #                                   as.Date("2017-05-01"),
  #                                   as.Date("2017-11-01"),
  #                                   as.Date("2018-05-01"),
  #                                   as.Date("2018-11-01"),
  #                                   as.Date("2019-04-30"))
  #     ) + labs(x="",y="Fitted Values")
  #     ModRes = ggplot(modFit
  #     ) + geom_point(aes(x=Date,y=Res),
  #                    size=2
  #     ) +scale_x_continuous(breaks=c(as.Date("2016-05-01"),
  #                                    as.Date("2016-11-01"),
  #                                    as.Date("2017-05-01"),
  #                                    as.Date("2017-11-01"),
  #                                    as.Date("2018-05-01"),
  #                                    as.Date("2018-11-01"),
  #                                    as.Date("2019-04-30"))
  #     ) + labs(x="",y="Residuals")
  #     HL = hoslem.test(masterLun$NightPres,fitted(lunMod))
  #     HLPlot = data.frame(Expected=HL$expected[,1],Observed=HL$observed[,1])
  #     HLPlot$X = seq(min(HL$observed[,1]),max(HL$observed[,1]),length.out=length(HLPlot$Expected))
  #     HLPlot$Y = seq(min(HL$observed[,1]),max(HL$observed[,1]),length.out=length(HLPlot$Expected))
  # 
  #     ObsPred = ggplot(HLPlot
  #     ) + geom_line(aes(x=X,y=Y)
  #     ) + geom_point(aes(y=Expected,x=Observed)
  #     ) + labs(x="Observed",y="Expected"
  #     ) + annotate("text",
  #                  label=paste('p-value: ',as.character(HL$p.value),sep=""),
  #                  size=4,
  #                  x=1.005*min(HL$observed[,1]),
  #                  y=0.999*max(HL$expected[,1]))
  #     png(file=paste(outDir,'/',CTname,'/',sites[k],"_LunModelFit.png",sep=""),width = 1400, height = 800, units = "px")
  #     grid.arrange(PresPlot,FitPlot,ModRes,ObsPred,BinRes,ncol=1,nrow=5,top=paste(CTname,'at',sites[k]))
  #     while (dev.cur()>1) {dev.off()}
  # 
  #     png(file=paste(outDir,'/',CTname,'/',sites[j],"_LunBinResid.png",sep=""),width = 400, height = 300, units = "px")
  #     binnedplot(fitted(lunMod),resid(lunMod,type="response"))
  #     while (dev.cur()>1) {dev.off()}
  # 
  #     ## Bootstrap GEEGLM parameter estimates for later construction of confidence intervals ----------------
  # 
  #     BootstrapParameters2<-rmvnorm(10000, coef(lunMod), summary(lunMod)$cov.unscaled)
  #     MoonPhaseBootstrapCoefs<- BootstrapParameters2[,MPInd]
  #     MoonAltBootstrapCoefs<- BootstrapParameters2[,AltInd]
  #     MoonPresBootstrapCoefs<- BootstrapParameters2[,PresInd]
  # 
  #     # Predict presence at each X value using model coefficients
  #     JxMP = model.matrix(lunMod)[,MPInd]*coef(lunMod)[MPInd]
  #     JxAlt = model.matrix(lunMod)[,AltInd]*coef(lunMod)[AltInd]
  #     JxPres = model.matrix(lunMod)[,PresInd]*coef(lunMod)[PresInd]
  # 
  #     ## MoonPhase ---------------------------
  #     if (!isempty(MPInd)){
  #       MoonPhaseForPlotting<- seq(0, 1, length=5000)
  #       MPBasis<- mSpline(MoonPhaseForPlotting,
  #                         knots=quantile(masterLun$MoonPhase,probs=c(0.275,0.5,0.725)),
  #                         Boundary.knots=c(0,1),
  #                         periodic=T)
  #       RealFitMP<- MPBasis%*%coef(lunMod)[MPInd] # multiply basis functions by model coefficients to get values of spline at each X
  #       RealFitCenterMP<- RealFitMP+coef(lunMod)[1]-mean(JxMP)# adjust offset
  #       RealFitCenterMP = inv.logit(RealFitCenterMP)
  #       MPBootstrapFits<- (MPBasis%*%t(MoonPhaseBootstrapCoefs))+coef(lunMod)[1] # get spread of spline values at each X based on distributions of each coefficient
  #       quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
  #       cisMP<-apply(MPBootstrapFits, 1, quant.func)-mean(JxMP) # confidence interval of smooth function estimate
  #       MPcil = inv.logit(cisMP[1,])
  #       MPciu = inv.logit(cisMP[2,])
  # 
  #       MPplotDF = data.frame(MoonPhaseForPlotting,RealFitCenterMP)
  #       colnames(MPplotDF) = c("MoonPhase","Fit")
  # 
  #       MP = ggplot(MPplotDF, aes(MoonPhase, Fit),
  #       ) + geom_smooth(fill = "grey",
  #                       colour = "black",
  #                       aes(ymin=MPcil, ymax=MPciu),
  #                       stat ="identity"
  #       ) + labs(x = "Moon Phase",
  #                y = "Probability",
  #       ) + theme(axis.line = element_line(size=0.2),
  #                 panel.background = element_blank()
  #       )
  # 
  #       saveName = paste(outDir,'/',CTname,'/',sites[k],"_",int,"_GEEGLM_MoonPhasePlot_Probs.png",sep="")
  #       ggsave(saveName,device="png", width=2, scale=3, height=0.5, units="in",dpi=600)
  #       while (dev.cur()>1) {dev.off()}
  #       saveName = paste(outDir,'/',CTname,'/',sites[k],"_",int,"_GEEGLM_MoonPhasePlot_Probs.pdf",sep="")
  #       ggsave(saveName,device="pdf", width=2, scale=3, height=0.5, units="in",dpi=600)
  #       while (dev.cur()>1) {dev.off()}
  # 
  #       MPdens =  ggplot(masterLun,aes(x=MoonPhase)
  #       )+geom_histogram(aes(y=..ncount..),
  #                        fill='#66B2FF',
  #                        binwidth = 0.5,
  #                        alpha=0.5
  #       ) + scale_x_continuous(breaks=c(0,0.5,0.99),
  #                              label=c("0","0.5","1")
  #       ) + labs(y="Normalized Count",x=NULL
  #       ) + theme_minimal()
  #       saveName = paste(outDir,'/',CTname,'/',sites[k],"_",int,"_GEEGLM_MPDataDensity.png",sep="")
  #       ggsave(saveName,device="png", width=4, scale=4, height=0.5, units="in",dpi=600)
  #       while (dev.cur()>1) {dev.off()}
  #       saveName = paste(outDir,'/',CTname,'/',sites[k],"_",int,"_GEEGLM_MPDataDensity.pdf",sep="")
  #       ggsave(saveName,device="pdf", width=4, scale=4, height=0.5, units="in",dpi=600)
  #       while (dev.cur()>1) {dev.off()}
  #     }
  #     ## Moon Altitude ---------------------------
  #     if (!isempty(AltInd)){
  #       MoonAltForPlotting<- seq(min(masterLun$MoonAltitude),max(masterLun$MoonAltitude), length=5000)
  #       AltBasis<- mSpline(MoonAltForPlotting,
  #                          knots=quantile(masterLun$MoonAltitude,probs=c(0.333,0.666)),
  #                          Boundary.knots=c(min(masterLun$MoonAltitude),max(masterLun$MoonAltitude)))
  #       RealFitAlt<- AltBasis%*%coef(lunMod)[AltInd] # multiply basis functions by model coefficients to get values of spline at each X
  #       RealFitCenterAlt<- RealFitAlt+coef(lunMod)[1]-mean(JxAlt)# adjust offset
  #       RealFitCenterAlt = inv.logit(RealFitCenterAlt)
  #       AltBootstrapFits<- (AltBasis%*%t(MoonAltBootstrapCoefs))+coef(lunMod)[1] # get spread of spline values at each X based on distributions of each coefficient
  #       quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
  #       cisAlt<-apply(AltBootstrapFits, 1, quant.func)-mean(JxAlt) # confidence interval of smooth function estimate
  #       Altcil = inv.logit(cisAlt[1,])
  #       Altciu = inv.logit(cisAlt[2,])
  # 
  #       AltplotDF = data.frame(MoonAltForPlotting,RealFitCenterAlt)
  #       colnames(AltplotDF) = c("MoonAlt","Fit")
  # 
  #       MoonAlt = ggplot(AltplotDF, aes(MoonAlt, Fit),
  #       ) + geom_smooth(fill = "grey",
  #                       colour = "black",
  #                       aes(ymin=Altcil, ymax=Altciu),
  #                       stat ="identity"
  #       ) + labs(x = "Moon Altitude",
  #                y = "Probability",
  #       ) + theme(axis.line = element_line(size=0.2),
  #                 panel.background = element_blank()
  #       )
  #       saveName = paste(outDir,'/',CTname,'/',sites[k],"_",int,"_GEEGLM_MoonAltPlot_Probs.png",sep="")
  #       ggsave(saveName,device="png", width=2, scale=3, height=0.5, units="in",dpi=600)
  #       while (dev.cur()>1) {dev.off()}
  #       saveName = paste(outDir,'/',CTname,'/',sites[k],"_",int,"_GEEGLM_MoonAltPlot_Probs.pdf",sep="")
  #       ggsave(saveName,device="pdf", width=2, scale=3, height=0.5, units="in",dpi=600)
  #       while (dev.cur()>1) {dev.off()}
  # 
  #       Altdens =  ggplot(masterLun,aes(x=MoonAltitude)
  #       )+geom_histogram(aes(y=..ncount..),
  #                        fill='#66B2FF',
  #                        binwidth = 0.5,
  #                        alpha=0.5
  #       ) + scale_x_continuous(breaks=c(0,45,90),
  #                              label=c("0","45","90")
  #       ) + labs(y="Normalized Count",x=NULL
  #       ) + theme_minimal()
  #       saveName = paste(outDir,'/',CTname,'/',sites[k],"_",int,"_GEEGLM_MoonAltDataDensity.png",sep="")
  #       ggsave(saveName,device="png", width=4, scale=4, height=0.5, units="in",dpi=600)
  #       while (dev.cur()>1) {dev.off()}
  #       saveName = paste(outDir,'/',CTname,'/',sites[k],"_",int,"_GEEGLM_MoonAltDataDensity.pdf",sep="")
  #       ggsave(saveName,device="pdf", width=4, scale=4, height=0.5, units="in",dpi=600)
  #       while (dev.cur()>1) {dev.off()}
  #     }
  #     ## Moon Presence ---------------------------
  #     if (!isempty(PresInd)){
  #       if (sites[j]=="HAT"){
  #         AdjustedMoonPresCoefs = data.frame(c(MoonPresBootstrapCoefs[,1],
  #                                              MoonPresBootstrapCoefs[,2]+mean(MoonPresBootstrapCoefs[,1])),
  #                                            as.factor(rep(1:3,each=10000)))
  #       } else {
  #         AdjustedMoonPresCoefs = data.frame(c(MoonPresBootstrapCoefs[,1],
  #                                              MoonPresBootstrapCoefs[,2]+mean(MoonPresBootstrapCoefs[,1]),
  #                                              MoonPresBootstrapCoefs[,3]+mean(MoonPresBootstrapCoefs[,1])),
  #                                            as.factor(rep(1:3,each=10000)))
  #       }
  #       colnames(AdjustedMoonPresCoefs) = c("Coefficient","Presence")
  #       AdjustedMoonPresCoefs$Prob = inv.logit(AdjustedMoonPresCoefs$Coefficient)
  # 
  #       # calculate quantiles for plotting limits
  #       quants = AdjustedMoonPresCoefs %>%
  #         group_by(Presence) %>%
  #         summarize(q25 = quantile(Prob,probs=0.25),
  #                   q75 = quantile(Prob,probs=0.75))
  #       iqr = quants$q75-quants$q25
  # 
  #       MPres = ggplot(AdjustedMoonPresCoefs,aes(Presence,Prob)
  #       ) + geom_boxplot(outlier.shape=NA,
  #                        lwd=0.2
  #       ) + labs(x='MoonPresence',y="Probability"
  #       ) + coord_cartesian(ylim = c(0,(1.5*max(iqr))+max(quants$q75))
  #       ) + scale_x_discrete(breaks=c(1,2,3),label=c("Before","MoonUp","After")
  #       ) + theme(axis.line = element_line(size=0.2),
  #                 panel.background = element_blank())
  #       saveName = paste(outDir,'/',CTname,'/',sites[k],"_",int,"_GEEGLM_MoonPresPlot_Probs.png",sep="")
  #       ggsave(saveName,device="png", width=2, scale=3, height=0.75, units="in",dpi=600)
  #       while (dev.cur()>1) {dev.off()}
  #       saveName = paste(outDir,'/',CTname,'/',sites[k],"_",int,"_GEEGLM_MoonPresPlot_Probs.pdf",sep="")
  #       ggsave(saveName,device="pdf", width=2, scale=3, height=0.75, units="in",dpi=600)
  #       while (dev.cur()>1) {dev.off()}
  # 
  #       Presdens =  ggplot(masterLun,aes(x=as.numeric(as.factor(MoonPres)))
  #       )+geom_histogram(aes(y=..ncount..),
  #                        fill='#66B2FF',
  #                        binwidth = 0.5,
  #                        alpha=0.5
  #       ) + scale_x_discrete(breaks=c(1,2,3),
  #                            label=c("Pre","MoonUp","Post")
  #       ) + labs(y="Normalized Count",x=NULL
  #       ) + theme_minimal()
  #       saveName = paste(outDir,'/',CTname,'/',sites[k],"_",int,"_GEEGLM_MoonPresDataDensity.png",sep="")
  #       ggsave(saveName,device="png", width=4, scale=4, height=0.5, units="in",dpi=600)
  #       while (dev.cur()>1) {dev.off()}
  #       saveName = paste(outDir,'/',CTname,'/',sites[k],"_",int,"_GEEGLM_MoonPresDataDensity.pdf",sep="")
  #       ggsave(saveName,device="pdf", width=4, scale=4, height=0.5, units="in",dpi=600)
  #       while (dev.cur()>1) {dev.off()}
  #     }
  #   }
  }
}

# #### PLOT ALL BEST FIT SPLINE MODELS FOR A GIVEN SPECIES --------------------------------------

# species = list.dirs(outDir,recursive=FALSE)
# 
# for ( i in 1:numel(species)){
# 
# 
#   JDmodFiles = list.files(path=species[i],pattern="*Model.Rdata",
#                         full.names=TRUE,recursive=FALSE,include.dirs=FALSE,no..=TRUE)
#   CTname = str_remove(species[i],paste(outDir,'/',sep=""))
#   sites = list()
# 
#   # JDList = list()
#   # LunIllumList = list()
#   # YrList = list()
# 
#   for (j in 1:numel(JDmodFiles)){
# 
#     site = str_remove(JDmodFiles[j],paste(outDir,"/",CTname,"/",sep=""))
#     sites = c(sites,str_remove(site,"_5minBin_Model.Rdata"))
# 
#     load(JDmodFiles[j]) # load JD + Yr model
#     
#     # find associated master dataframe
#     thisSpec = which(str_detect(dfList,CTname))
#     atSite = which(str_detect(dfList,unlist(sites[j])))
#     thisModInd = intersect(thisSpec,atSite)
#     thisSite = data.frame(read.csv(dfList[thisModInd]))
# 
# 
#     # Get indices of coefficients for each covar
#     if (numel(JDKnots)==2){
#       JDInd = c(2,3)
#     } else {JDInd = c(2:4)}
#     if (numel(LunKnots)==2){
#       LunInd = JDInd[length(JDInd)]+c(1:5)
#     } else if (numel(LunKnots)==3){
#       LunInd = JDInd[length(JDInd)]+c(1:6)
#     }
#     YrInd = c(1,LunInd[length(LunInd)]+c(1:2))
#     #PhsInd = c(1,YrInd[length(YrInd)]+c(1:3))
# 
# 
#     ## Bootstrap GEEGLM parameter estimates for later construction of confidence intervals ----------------
#     BootstrapParameters<-rmvnorm(10000, coef(tempMod), summary(tempMod)$cov.unscaled)
#     JDayBootstrapCoefs<- BootstrapParameters[,JDInd]
#     LunBootstrapCoefs<- BootstrapParameters[,LunInd]
#     YearBootstrapCoefs<- BootstrapParameters[,YrInd]
#     # PhsBootstrapCoefs<- BootstrapParameters[,PhsInd]
#     # if (j==7){ # if HAT, get coefficients for Site
#     #   hatSiteBootstrapCoefs<- BootstrapParameters[,c(1,7)]
#     #   quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
#     #   cisH<-apply(hatSiteBootstrapCoefs, 2, quant.func)
#     # }
# 
#     # Predict presence at each X value using model coefficients
#     # if (j==7){ # if HAT, include Site as factor
#     #   JxJD<-model.matrix(tempMod)[,2:3]%*%coef(tempMod)[c(2:3)]
#     #   Jx2<-model.matrix(tempMod)[,c(1,4:6)]%*%coef(tempMod)[c(1,4:6)]
#     #   Jx3<-model.matrix(tempMod)[,c(1,7)]%*%coef(tempMod)[c(1,7)]
#     # } else {
#     JxJD = model.matrix(tempMod)[,JDInd]%*%coef(tempMod)[JDInd]
#     if (length(LunInd)>1){
#       JxLun = model.matrix(tempMod)[,LunInd]%*%coef(tempMod)[LunInd]
#     } else if (length(LunInd)==1){
#       JxLun = model.matrix(tempMod)[,LunInd]*coef(tempMod)[LunInd]
#     }
#     JxYr = model.matrix(tempMod)[,YrInd]%*%coef(tempMod)[YrInd]
#     # }
# 
# 
#     ### Generate GEEGLM partial residual plots ---------------------------
#     ## Julian Day ---------------------------
#     JDayForPlotting<- seq(min(thisSite$JulianDay), max(thisSite$JulianDay), length=5000)
#     JBasis<- mSpline(JDayForPlotting,  # spline spanning range of X values
#                                    knots=quantile(thisSite$JulianDay, probs=JDKnots),
#                                    Boundary.knots=c(1,365),
#                                    periodic=T) # basis functions for smooth function
#     RealFitJ<- JBasis%*%coef(tempMod)[JDInd] # multiply basis functions by model coefficients to get values of spline at each X
#     RealFitCenterJ<- RealFitJ-mean(JxJD)#-coef(tempMod)[1] # adjust offset
#     RealFitCenterJ<- inv.logit(RealFitCenterJ)
#     JDayBootstrapFits<- JBasis%*%t(JDayBootstrapCoefs) # get spread of spline values at each X based on distributions of each coefficient
#     quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
#     cisJ<-apply(JDayBootstrapFits, 1, quant.func)-mean(JxJD) # confidence interval of smooth function estimate
#     # Jcil<-cisJ[1,]-mean(JxJD)-coef(tempMod)[1] # lowerCI bound
#     # Jciu<-cisJ[2,]-mean(JxJD)-coef(tempMod)[1] # upper CI bound
#     Jcil<-inv.logit(cisJ[1,]) # lowerCI bound
#     Jciu<-inv.logit(cisJ[2,]) # upper CI bound
# 
# 
#     JplotDF = data.frame(JDayForPlotting,RealFitCenterJ)
#     colnames(JplotDF) = c("Jday","Fit")
#     # dJday = stats::density(thisSite$JulianDay,na.rm = TRUE,n=256,from=1,to=365) # Calculate kernel density of Jday observations
#     # Jdens = data.frame(c(dJday$x,rev(dJday$x)),c(dJday$y,rep(0,length(dJday$y))))
#     # colnames(Jdens) = c("Day","Density")
#     # Jdens$Density = rescale(Jdens$Density, to=c(0.95*min(Jcil),min(Jcil))) # rescale density to sit at bottom of y-axis
#     densDF = data.frame(JD = seq(1,365,1))
#     whichDay = histc(thisSite$JulianDay,seq(1,366,1)) # Count up observations on each Julian day
#     densDF$normDens = (whichDay$cnt-min(whichDay$cnt))[1:365] # set min to 0
#     densDF$normDens = densDF$normDens/864#max(densDF$normDens) # normalize by max # obs possible (288 5-min bins per day, times 3 days)
# 
#     JD = ggplot(JplotDF, aes(Jday, Fit),
#     # ) + geom_polygon(data=Jdens,
#     #                  aes(Day,Density),
#     #                  fill=4,
#     #                  alpha=0.2
#     ) + geom_smooth(fill = "grey",
#                     colour = "black",
#                     aes(ymin=Jcil, ymax=1.05*Jciu),
#                     stat ="identity"
#     ) + labs(x = "Julian Day",
#              y = "Probability",
#              #title = paste(CTname, 'at',site),
#     ) + scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),
#                            label=c("J","F","M","A","M","J","J","A","S","O","N","D")
#     ) + theme(axis.line = element_line(size=0.2),
#               panel.background = element_blank()
#     )
#     # if (PV$'p-value'[1]<0.05){
#     #   JD = JD + annotate("text",x=0.97*365,y=max(Jciu),label="*",size=12)
#     # }
# 
#     # JD + JDdens + plot_layout(widths=c(3,1))
#     # saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDayPlot_Probs.png",sep="")
#     # ggsave(saveName,device="png", width=2, scale=3, height=0.5, units="in",dpi=600)
#     saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDayPlot_Probs.pdf",sep="")
#     ggsave(saveName,device="pdf", width=2, scale=3, height=0.5, units="in",dpi=600)
#     while (dev.cur()>1) {dev.off()}
# 
#     JDdens = ggplot(densDF)+geom_area(aes(JD,normDens),
#                                       fill=4
#     ) + scale_x_continuous(breaks=c(1,182,335),
#                            label=c("Jan","Jun","Dec")
#     ) + labs(y="Data Density",x=NULL
#     ) + theme_minimal()
#     # saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDDataDensity.png",sep="")
#     # ggsave(saveName,device="png", width=4, scale=4, height=0.5, units="in",dpi=600)
#     saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDDataDensity.pdf",sep="")
#     ggsave(saveName,device="pdf", width=1, scale=4, height=0.5, units="in",dpi=600)
# 
#     # assign(paste("JD_",j,sep=""),JD)
#     # JDList = c(JDList,assign(paste("JD_",j,sep=""),JD))
# 
#     ## Lunar Illuminance ---------------------------
#     LunIllumForPlotting<- seq(0, 1, length=5000)
#     if (numel(LunKnots)>1){
#       LBasis<- mSpline(LunIllumForPlotting,
#                        knots=quantile(thisSite$LunarIllum,probs=LunKnots),
#                        Boundary.knots=c(0,1))
#       RealFitL<- LBasis%*%coef(tempMod)[LunInd] # multiply basis functions by model coefficients to get values of spline at each X
#     } else {LBasis = LunIllumForPlotting
#     RealFitL<- LBasis*coef(tempMod)[LunInd]}
#     RealFitCenterL<- RealFitL-mean(JxLun)#-coef(tempMod)[1] # adjust offset
#     RealFitCenterL = inv.logit(RealFitCenterL)
#     LunBootstrapFits<- LBasis%*%t(LunBootstrapCoefs) # get spread of spline values at each X based on distributions of each coefficient
#     quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
#     cisL<-apply(LunBootstrapFits, 1, quant.func)-mean(JxLun) # confidence interval of smooth function estimate
#     # Lcil<-cisL[1,]-mean(JxLun)-coef(tempMod)[1] # lower CI bound
#     # Lciu<-cisL[2,]-mean(JxLun)-coef(tempMod)[1] # upper CI bound
#     Lcil = inv.logit(cisL[1,])
#     Lciu = inv.logit(cisL[2,])
# 
# 
#     LplotDF = data.frame(LunIllumForPlotting,RealFitCenterL)
#     colnames(LplotDF) = c("LunIllum","Fit")
#     # dLun = stats::density(thisSite$LunarIllum,na.rm = TRUE,n=256,from=0,to=1) # Calculate kernel density of LunIllum observations
#     # Ldens = data.frame(c(dLun$x,rev(dLun$x)),c(dLun$y,rep(0,length(dLun$y))))
#     # colnames(Ldens) = c("Illumination","Density")
#     # if (min(Lcil)>0){
#     #   Ldens$Density = rescale(Ldens$Density, to=c(0.95*min(Lcil),min(Lcil))) # rescale density to sit at bottom of y-axis
#     # } else {
#     #   Ldens$Density = rescale(Ldens$Density, to=c(1.05*min(Lcil),min(Lcil)))
#     # }
#     densDF = data.frame(LunFrac = seq(0,0.99,by=0.01))
#     whichDay = histc(thisSite$LunarIllum,seq(0,1.01,by=0.01)) # Count up observations in each lunar illuminance bin
#     densDF$normDens = (whichDay$cnt-min(whichDay$cnt))[1:100] # set min to 0
#     densDF$normDens = densDF$normDens/max(densDF$normDens) # normalize by max # obs possible (288 5-min bins per day, times 3 days)
# 
#     # dens$Density = dens$Density/max(dens$Density) # normalize kernel density
#     # if (min(Lcil)<0){ # set kernel density at bottom of y axis
#     #   dens$Density = dens$Density - abs(min(Lcil))
#     # } else {
#     #   dens$Density = dens$Density + min(Lcil)
#     # }
# 
#     if (numel(LunKnots)>1){ # plot smooth
#       LunIllum = ggplot(LplotDF, aes(LunIllum, Fit),
#       # ) + geom_polygon(data=Ldens,
#       #                  aes(Illumination,Density),
#       #                  fill=4,
#       #                  alpha=0.2
#       ) + geom_smooth(fill = "grey",
#                       colour = "black",
#                       aes(ymin=Lcil, ymax=Lciu),
#                       stat ="identity"
#       ) + labs(x = "Lunar Illuminance",
#                y = "Probability",
#                #title = paste(CTname, 'at',site),
#       ) + theme(axis.line = element_line(size=0.2),
#                 panel.background = element_blank()
#       ) } else { # plot linear
#         LunIllum = ggplot(LplotDF, aes(LunIllum, Fit),
#         ) + geom_polygon(data=Ldens,
#                          aes(Illumination,Density),
#                          fill=4,
#                          alpha=0.2
#         ) + geom_smooth(fill = "grey",
#                         colour = "black",
#                         aes(ymin=Lcil, ymax=Lciu),
#                         stat ="identity"
#         ) + labs(x = "Lunar Illuminance",
#                  y = "Probability",
#                  #title = paste(CTname, 'at',site),
#         ) + theme(axis.line = element_line(size=0.2),
#                   panel.background = element_blank()
#         )
#       }
#     # if (PV$'p-value'[2]<0.05){
#     #   LunIllum = LunIllum + annotate("text",x=0.97,y=0.95*max(Lciu),label="*",size=12)
#     # }
#     # saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_LunIllumPlot_Probs.png",sep="")
#     # ggsave(saveName,device="png", width=2, scale=3, height=0.5, units="in",dpi=600)
#     saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_LunIllumPlot_Probs.pdf",sep="")
#     ggsave(saveName,device="pdf", width=2, scale=3, height=0.5, units="in",dpi=600)
#     while (dev.cur()>1) {dev.off()}
# 
#     # assign(paste("LunIllum_",j,sep=""),LunIllum)
#     # LunIllumList = c(LunIllumList,assign(paste("LunIllum_",j,sep=""),LunIllum))
# 
#     LIDdens = ggplot(densDF)+geom_area(aes(LunFrac,normDens),
#                                       fill=4
#     ) + scale_x_continuous(breaks=c(0,0.5,0.99),
#                            label=c("0","0.5","1")
#     ) + labs(y="Data Density",x=NULL
#     ) + theme_minimal()
#     # saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_LIDataDensity.png",sep="")
#     # ggsave(saveName,device="png", width=4, scale=4, height=0.5, units="in",dpi=600)
#     saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_LIDataDensity.pdf",sep="")
#     ggsave(saveName,device="pdf", width=4, scale=4, height=0.5, units="in",dpi=600)
# 
# 
#     ## Year (as boxplot) -----------------------------------
#     # Center intercept (1st level of year factor) at 0 and show other levels relative to it
#     AdjustedYearCoefs = data.frame(c(YearBootstrapCoefs[,1]-mean(YearBootstrapCoefs[,1]),
#                                      YearBootstrapCoefs[,2],
#                                      YearBootstrapCoefs[,3]),
#                                    # as.factor(rep(2016:2019,each=10000)))
#                                    as.factor(rep(1:3,each=10000)))
#     colnames(AdjustedYearCoefs) = c("Coefficient","Year")
#     # AdjustedYearCoefs$Prob = AdjustedYearCoefs$Coefficient*as.numeric(AdjustedYearCoefs$Year) # go from coefficients to probability of presence
#     AdjustedYearCoefs$Prob = inv.logit(AdjustedYearCoefs$Coefficient)
# 
#     # calculate quantiles for plotting limits
#     quants = AdjustedYearCoefs %>%
#       group_by(Year) %>%
#       summarize(q25 = quantile(Prob,probs=0.25),
#                 q75 = quantile(Prob,probs=0.75))
#     iqr = quants$q75-quants$q25
# 
#     Yr = ggplot(AdjustedYearCoefs,aes(Year,Prob)
#     ) + geom_boxplot(outlier.shape=NA,
#                      lwd=0.2
#     ) + labs(x='Year',y="Probability"
#     ) + coord_cartesian(ylim = c(min(quants$q25-(1.75*max(iqr))),(1.75*max(iqr))+max(quants$q75))
#     ) + theme(axis.line = element_line(size=0.2),
#              panel.background = element_blank()
#     ) #+ labs(title = paste(CTname, 'at',site))
#     # if (PV$'p-value'[3]<0.05){
#     #   Yr = Yr + annotate("text",x=3,y=0.95*max(AdjustedYearCoefs$Prob),label="*",size=12)
#     # }
#     # saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_YearPlot_Probs.png",sep="")
#     # ggsave(saveName,device="png", width=2, scale=3, height=1, units="in",dpi=600)
#     saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_YearPlot_Probs.pdf",sep="")
#     ggsave(saveName,device="pdf", width=2, scale=3, height=0.75, units="in",dpi=600)
#     while (dev.cur()>1) {dev.off()}

# assign(paste("Yr_",j,sep=""),Yr)
# YrList = c(YrList,assign(paste("Yr_",j,sep=""),Yr))

## Day Phase (as boxplot) -----------------------------
# # Center intercept (1st level of day phase factor) at 0 and show other levels relative to it
# AdjustedPhsCoefs = data.frame(c(PhsBootstrapCoefs[,1]-mean(PhsBootstrapCoefs[,1]),
#                                 PhsBootstrapCoefs[,2],
#                                 PhsBootstrapCoefs[,3],
#                                 PhsBootstrapCoefs[,4]),
#                                as.factor(rep(c("Dawn","Day","Dusk","Night"),each=10000)))
# colnames(AdjustedPhsCoefs) = c("Coefficient","Phase")
# # AdjustedPhsCoefs = apply(AdjustedPhsCoefs,2,exp) return to units of response var?
#
# # calculate quantiles for plotting limits
# quants = AdjustedPhsCoefs %>%
#   group_by(Phase) %>%
#   summarize(q25 = quantile(Coefficient,probs=0.25),
#             q75 = quantile(Coefficient,probs=0.75))
# iqr = quants$q75-quants$q25
#
# #saveName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_YearPlot.png",sep="")
# Phs = ggplot(AdjustedPhsCoefs,aes(Phase,Coefficient)
# ) + geom_boxplot(outlier.shape=NA
# )+ coord_cartesian(ylim = c(min(quants$q25-(1.75*max(iqr))),(1.75*max(iqr))+max(quants$q75))
# )+ theme(axis.line = element_line(),
#          panel.background = element_blank()
# )# + labs(title = paste(CTname, 'at',site))
# # ggsave(saveName,device="png")
# # while (dev.cur()>1) {dev.off()}

# # Site (as boxplot), if HAT
# if (j==7){
#   AdjustedHATCoefs = data.frame(c(hatSiteBootstrapCoefs[,1]-mean(hatSiteBootstrapCoefs[,1]),
#                                   hatSiteBootstrapCoefs[,2]),
#                                 as.factor(rep(c("HAT_A","HAT_B"),each=10000)))
#   colnames(AdjustedHATCoefs) = c("Coefficient","Site")
#
#   saveName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_HATSitePlot.png",sep="")
#   ggplot(AdjustedHATCoefs,aes(Site,Coefficient)
#   ) + geom_boxplot(
#   )+ theme(axis.line = element_line(),
#            panel.background = element_blank()
#   )+labs(title = paste(CTname, 'at',sites[j]))
#   ggsave(saveName,device="png")
#   while (dev.cur()>1) {dev.off()}
#
# }

# }


#### ARRANGE PLOTS AND SAVE ---------------------

# JD + LunIllum + Yr + Phs +plot_annotation(
#   title=paste(CTname, 'at',site)
# )

# JD / LunIllum + plot_annotation(
#   title=paste(CTname, 'at',site)
# )
# saveName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLM_SmoothPlots.png",sep="")
#   ggsave(saveName,device="png")
#   while (dev.cur()>1) {dev.off()}
#
#   Yr + plot_annotation(
#     title=paste(CTname, 'at',site)
#   )
#   saveName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLM_YearPlots.png",sep="")
#   ggsave(saveName,device="png")
#   while (dev.cur()>1) {dev.off()}

# JD_1 / JD_2 / JD_3
#
#     wrap_plots(JDList,ncol=1) + plot_layout(guides='collect')
#   wrap_plots(LunIllumList,ncol=1)
#   wrap_plots(YrList,ncol=1)

# }



#### Old code for plotting Day Phase (as factor and as smooth) -----------------
# Day Phase (as smooth) ---------------------------
# DPForPlotting<- seq(min(thisSite$numDayPhase), max(thisSite$numDayPhase), length=5000)
# DPBasis<- mSpline(DPForPlotting,  # spline spanning range of X values
#                  knots=quantile(thisSite$numDayPhase, probs=c(0.333,0.666)),
#                  Boundary.knots=c(0.5,4.5),
#                  periodic=T) # basis functions for smooth function
# RealFitDP<- DPBasis%*%coef(tempMod)[DPInd] # multiply basis functions by model coefficients to get values of spline at each X
# RealFitCenterDP<- RealFitDP-mean(JxDP) # adjust offset
# RealFitCenterDP<- inv.logit(RealFitCenterDP)
# DPBootstrapFits<- DPBasis%*%t(DayPhaseBootstrapCoefs) # get spread of spline values at each X based on distributions of each coefficient
# quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
# cisDP<-apply(DPBootstrapFits, 1, quant.func)-mean(JxDP) # confidence interval of smooth function estimate
# DPcil<-inv.logit(cisDP[1,]) # lowerCI bound
# DPciu<-inv.logit(cisDP[2,]) # upper CI bound
# 
# DPplotDF = data.frame(DPForPlotting,RealFitCenterDP)
# colnames(DPplotDF) = c("DayPhase","Fit")
# 
# # # calculate data density
# # densDF = data.frame(DP = seq(1,4,1))
# # whichDP = histc(thisSite$numDayPhase,seq(1,5,1)) # Count up observations on each Julian day
# # densDF$normDens = (whichDP$cnt-min(whichDP$cnt))[1:4] # set min to 0
# # # normalize by max # obs possible (288 5-min bins per day, times 3 observations of each Julian day in the study period)
# # densDF$normDens = densDF$normDens/max(densDF$normDens)
# 
# DP = ggplot(DPplotDF, aes(DayPhase, Fit),
# ) + geom_smooth(fill = "grey",
#                 colour = "black",
#                 aes(ymin=DPcil, ymax=1.05*DPciu),
#                 stat ="identity"
# ) + labs(x = "Day Phase",
#          y = "Probability",
#          #title = paste(CTname, 'at',site),
# ) + scale_x_continuous(breaks=c(1,2,3,4),
#                        label=c("Dawn","Day","Dusk","Night")
# ) + theme(axis.line = element_line(size=0.2),
#           panel.background = element_blank()
# )
# 
# saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_DPSmoothPlot_Probs.pdf",sep="")
# ggsave(saveName,device="pdf", width=2, scale=3, height=0.5, units="in",dpi=600)
# while (dev.cur()>1) {dev.off()}
# 
# # DPdens = ggplot(densDF)+geom_area(aes(DP,normDens),
# #                                   fill=4
# # ) + scale_x_continuous(breaks=c(1,2,3,4),
# #                        label=c("Dawn","Day","Dusk","Night")
# # ) + labs(y="Data Density",x=NULL
# # ) + theme_minimal()
# DPdens = ggplot(thisSite)+geom_col(aes(numDayPhase,Presence)
# ) + scale_x_continuous(breaks=c(1,2,3,4),
#                        label=c("Dawn","Day","Dusk","Night")
# ) + labs(y="Data Density",x=NULL
# ) + theme_minimal()
# saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_DPSmoothDataDensity.pdf",sep="")
# ggsave(saveName,device="pdf", width=1, scale=4, height=0.5, units="in",dpi=600)

## Day Phase (as boxplot) -----------------------------
# # Center intercept (1st level of day phase factor) at 0 and show other levels relative to it
# AdjustedDPCoefs = data.frame(c(DayPhaseBootstrapCoefs[,1]-mean(DayPhaseBootstrapCoefs[,1]),
#                                DayPhaseBootstrapCoefs[,2],
#                                DayPhaseBootstrapCoefs[,3],
#                                DayPhaseBootstrapCoefs[,4]),
#                              as.factor(rep(c("Dawn","Day","Dusk","Night"),each=10000)))
# colnames(AdjustedDPCoefs) = c("Coefficient","DayPhase")
# AdjustedDPCoefs$Prob = inv.logit(AdjustedDPCoefs$Coefficient)
#
# # calculate quantiles for plotting limits
# quants = AdjustedDPCoefs %>%
#   group_by(DayPhase) %>%
#   summarize(q25 = quantile(Prob,probs=0.25),
#             q75 = quantile(Prob,probs=0.75))
# iqr = quants$q75-quants$q25
#
# sigComps = logical()
# if (PV$'p-value'[2]<0.05){ # ifday phase was significant, do pairwise comparison
#   KW = kruskal.test(AdjustedDPCoefs$Coefficient,AdjustedDPCoefs$DayPhase,na.action=na.pass)
#   pairCompCI = conover.test(AdjustedDPCoefs$Coefficient,AdjustedDPCoefs$DayPhase,method='bonferroni',label=TRUE,
#                             wrap=TRUE,table=TRUE,alpha=0.05)
#   my_comparisons <- list( c("Dawn", "Day"), c("Dawn", "Dusk"),c("Day","Dusk"),
#                           c("Dawn", "Night"), c("Day","Night"),c("Dusk","Night"))
#   sigComps = which(pairCompCI$P.adjusted<0.05)
#
#   ndim=numel(levels(AdjustedDPCoefs$DayPhase))-1
#   CITable = matrix(NA,ndim,ndim)
#   CITable[upper.tri(CITable,diag=TRUE)] = t(pairCompCI$P.adjusted)
#   CITable = t(CITable)
#   CITable = signif(CITable,digits=3)
#   Comparisons = matrix(NA,ndim,ndim)
#   Comparisons[upper.tri(Comparisons,diag=TRUE)] = t(pairCompCI$comparisons)
#   Comparisons = t(Comparisons)
#   CITable = rbind(CITable,Comparisons)
#   saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_DPCI_Comparisons.csv",sep="")
#   write.csv(CITable,file=saveName,row.names=TRUE)
# }
#
# Phs = ggplot(AdjustedDPCoefs,aes(DayPhase,Prob)
# ) + geom_boxplot(varwidth=TRUE,
#                  outlier.shape=NA,
#                  lwd=0.2
# )+ labs(x='Day Phase',y="Probability"
# )+ coord_cartesian(ylim = c(min(quants$q25-(1.75*max(iqr))),(3.5*max(iqr))+max(quants$q75))
# )+ theme(axis.line = element_line(size=0.2),
#          panel.background = element_blank()
# )+ geom_signif(comparisons=my_comparisons[sigComps],
#                y_position=(max(iqr)*c(1.5,1.8,2.1,2.4,2.7,3))+max(quants$q75),
#                test=wilcox.test,
#                textsize=6,
#                tip_length=0.01,
#                vjust=0.7,
#                map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05))
# saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_DPPlot_Probs.pdf",sep="")
# ggsave(saveName,device="pdf", width=2, scale=3, height=0.75, units="in",dpi=600)
# while (dev.cur()>1) {dev.off()}


