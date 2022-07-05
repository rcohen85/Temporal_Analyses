#### PLOT ALL FIXED SPLINE MODELS FOR A GIVEN SPECIES --------------------------------------
library(tidyverse)
library(stringr)
library(pracma)
library(broom)
library(lubridate)
library(splines2)
library(pracma)
library(SimDesign)
library(boot)
library(gridExtra)
library(ResourceSelection)
library(performance)
# library(fixest)
source("binned_residuals_RC.R")
# source("r2_tjur_RC.R")


modelDFDir = 'J:/Chpt_2/TimeSeries_ScaledByEffortError'
outDir = 'J:/Chpt_2/ModelOutput'
int = "5minBin"

dfList = list.files(path=modelDFDir,pattern=paste('*',int,'_MasterTempLun.csv',sep=""),
                    full.names=TRUE,recursive=FALSE,
                    include.dirs=FALSE,no..=TRUE)
# lunList = list.files(path=modelDFDir,pattern=paste('*',int,'_MasterLun.csv',sep=""),
#                      full.names=TRUE,recursive=FALSE,
#                      include.dirs=FALSE,no..=TRUE)

species = list.dirs(outDir,recursive=FALSE)
species = species[-which(!is.na(str_match(species,"plots")))]
modPerf = data.frame(Species=as.character(),
                     Site=as.character(),
                     PropGoodResid=as.numeric())

PresStats = list()
siteNames = c("HZ","OC","NC","BC","WC","NFC","HAT","GS","BP","BS","JAX")

for ( i in c(1:8,10:numel(species))){
  
  JDmodFiles = list.files(path=species[i],pattern="*Model_TempLun3.Rdata",
                          full.names=TRUE,recursive=FALSE,include.dirs=FALSE,no..=TRUE)
  # LunmodFiles = list.files(path=species[i],pattern="*Model_LunPhaseAltPres.Rdata",
  #                          full.names=TRUE,recursive=FALSE,include.dirs=FALSE,no..=TRUE)
  CTname = str_remove(species[i],paste(outDir,'/',sep=""))
  PresStats[[CTname]] = data.frame("ModPerf"=numeric(11),
                                   "JDSignif"=numeric(11),
                                   "JDPeak"=numeric(11),
                                   "MPhSignif"=numeric(11),
                                   "LunarPeak"=numeric(11),
                                   "NTSignif"=numeric(11),
                                   "NTJDSignif"=numeric(11),
                                   "NTMPhSignif"=numeric(11),
                                   "DielProp"=numeric(11),
                                   "DielPeak"=numeric(11))
  rownames(PresStats[[CTname]]) = siteNames
  sites = list()
  
  for (j in 1:numel(JDmodFiles)){
    
    site = str_remove(JDmodFiles[j],paste(outDir,"/",CTname,"/",sep=""))
    sites = c(sites,str_remove(site,"_5minBin_Model_TempLun3.Rdata"))
    if (sites[j]=="NULL"){
      stop("Didn't get site name")}
    siteInd = which(!is.na(str_match(siteNames,unlist(sites[j]))))
    
    load(JDmodFiles[j]) # load model
    
    # Save term significance
    if (any(PV$Variable=="JDs")) {
      if(PV$'p-value'[PV$Variable=="JDs"]=="<0.0001"){PresStats[[CTname]]$JDSignif[siteInd] = 0.0001
      } else {PresStats[[CTname]]$JDSignif[siteInd] = as.numeric(PV$'p-value'[PV$Variable=="JDs"])}
    } else {PresStats[[CTname]]$JDSignif[siteInd] = "NS" }
    if (any(PV$Variable=="NTs")) {
      if(PV$'p-value'[PV$Variable=="NTs"]=="<0.0001"){PresStats[[CTname]]$NTSignif[siteInd] = 0.0001
      } else {PresStats[[CTname]]$NTSignif[siteInd] = as.numeric(PV$'p-value'[PV$Variable=="NTs"])}
    } else {PresStats[[CTname]]$NTSignif[siteInd] = "NS" }
    if (any(PV$Variable=="NTs:JDs")) {
      if(PV$'p-value'[PV$Variable=="NTs:JDs"]=="<0.0001"){PresStats[[CTname]]$NTJDSignif[siteInd] = 0.0001
      } else {PresStats[[CTname]]$NTJDSignif[siteInd] = as.numeric(PV$'p-value'[PV$Variable=="NTs:JDs"])}
    } else {PresStats[[CTname]]$NTJDSignif[siteInd] = "NS" }
    if (any(PV$Variable=="NTs:MPhs")) {
      if(PV$'p-value'[PV$Variable=="NTs:MPhs"]=="<0.0001"){PresStats[[CTname]]$NTMPhSignif[siteInd] = 0.0001
      } else {PresStats[[CTname]]$NTMPhSignif[siteInd] = as.numeric(PV$'p-value'[PV$Variable=="NTs:MPhs"])}
    } else {PresStats[[CTname]]$NTMPhSignif[siteInd] = "NS" }
    
    if (tempMod$geese$error==0){ # not all models converged, don't plot non-converged models
      
      # find associated master dataframe
      thisSpec = which(str_detect(dfList,CTname))
      atSite = which(str_detect(dfList,unlist(sites[j])))
      thisModInd = intersect(thisSpec,atSite)
      thisSite = data.frame(read.csv(dfList[thisModInd])) # load master data frame
      
      # Indices of coefficients for each covar
      JDInd = numeric()
      MPhInd = numeric()
      NTInd = numeric()
      # LFInd = numeric()
      YrInd = numeric()
      JDNTInd = numeric()
      NTMPhInd = numeric()
      # MPhLFInd = numeric()
      JDInd = which(str_detect(names(tempMod$coefficients),"JDs")&!str_detect(names(tempMod$coefficients),"NTs"))
      MPhInd = which(str_detect(names(tempMod$coefficients),"MPhs")&!str_detect(names(tempMod$coefficients),"NTs")&!str_detect(names(tempMod$coefficients),"LF"))
      NTInd = which(str_detect(names(tempMod$coefficients),"NTs")&!str_detect(names(tempMod$coefficients),"JDs")&!str_detect(names(tempMod$coefficients),"MPhs"))
      # LFInd = which(str_detect(names(tempMod$coefficients),"LF")&!str_detect(names(tempMod$coefficients),"MPhs"))
      YrInd = which(str_detect(names(tempMod$coefficients),"YrF"))
      if (!isempty(YrInd)){
        YrInd = c(1,YrInd)
      }
      JDNTInd = which(str_detect(names(tempMod$coefficients),"JDs")&str_detect(names(tempMod$coefficients),"NTs"))
      NTMPhInd = which(str_detect(names(tempMod$coefficients),"MPhs")&str_detect(names(tempMod$coefficients),"NTs"))
      # MPhLFInd = which(str_detect(names(tempMod$coefficients),"MPhs")&str_detect(names(tempMod$coefficients),"LF"))
      
      
      if (sites[j]=="HAT"){
        startInd = which(thisSite$TimeStamp>=as.POSIXct('2017-05-01 00:00:00',format="%Y-%m-%d %H:%M:%S",tz="GMT"))
        thisSite = thisSite[startInd,]
      }
      
      # calculation and plotting of binned residuals stopped working after a system update, unable to debug so far
      # binRes = binned_residuals_RC(tempMod)
      # resid_ok = sum(binRes$group == "yes")/length(binRes$group)
      # thisModPer = cbind(CTname,sites[j],round(resid_ok,digits=5))
      # modPerf = rbind(modPerf,thisModPer)
      # PresStats[[CTname]]$ModPerf[siteInd] = unlist(thisModPer[3])
      # 
      # png(file=paste(outDir,'/',CTname,'/',sites[j],"_BinResidTL.png",sep=""),width = 400, height = 300, units = "px")
      # print(binned_residuals_RC(tempMod))
      # while (dev.cur()>1) {dev.off()}
      
      ## Bootstrap GEEGLM parameter estimates for later construction of confidence intervals ----------------
      BootstrapParameters1<-rmvnorm(10000, coef(tempMod), summary(tempMod)$cov.unscaled)
      NTBootstrapCoefs<- BootstrapParameters1[,NTInd]
      JDayBootstrapCoefs<- BootstrapParameters1[,JDInd]
      MPhBootstrapCoefs<- BootstrapParameters1[,MPhInd]
      # LFBootstrapCoefs<- BootstrapParameters1[,LFInd]
      YearBootstrapCoefs<- BootstrapParameters1[,YrInd]
      JDNTBootstrapCoefs<- BootstrapParameters1[,JDNTInd]
      NTMPhBootstrapCoefs<- BootstrapParameters1[,NTMPhInd]
      # MPhLFBootstrapCoefs<- BootstrapParameters1[,MPhLFInd]
      
      ## Recreate basis functions for smooths
      JDayForPlotting<- seq(min(thisSite$JulianDay), max(thisSite$JulianDay), length=1000)
      JDBasis<- mSpline(JDayForPlotting, 
                        knots=quantile(thisSite$JulianDay, probs=c(0.275,0.5,0.725)),
                        Boundary.knots=c(1,365),
                        periodic=T) 
      NTForPlotting<- seq(min(thisSite$NormTime), max(thisSite$NormTime), length=1000)
      NTBasis<- mSpline(NTForPlotting,  
                        knots=quantile(thisSite$NormTime, probs=c(0.275,0.5,0.725)),
                        Boundary.knots=c(-1,1),
                        periodic=T)
      MPhForPlotting<- seq(min(thisSite$MoonPhase), max(thisSite$MoonPhase), length=1000)
      MPhBasis<- mSpline(MPhForPlotting,
                         knots=quantile(thisSite$MoonPhase, probs=c(0.275,0.5,0.725)),
                         Boundary.knots=c(0,1),
                         periodic=T)
      
      # Calculate each smooth term at its mean value
      if (!isempty(JDInd)){JDmean = as.vector(apply(JDBasis,MARGIN=2,mean)%*%coef(tempMod)[JDInd])
      } else {JDmean = 0}
      if (!isempty(NTInd)){NTmean = as.vector(apply(NTBasis,MARGIN=2,mean)%*%coef(tempMod)[NTInd])
      } else {NTmean = 0}
      if (!isempty(MPhInd)) {MPhmean = as.vector(apply(MPhBasis,MARGIN=2,mean)%*%coef(tempMod)[MPhInd])
      } else {MPhmean=0}
      
      
      ### Generate GEEGLM partial residual plots ---------------------------
      # Julian Day ---------------------------
      if (!isempty(JDInd) & isempty(JDNTInd)){
        Fit<- JDBasis%*%(coef(tempMod)[JDInd]) # multiply basis functions by model coefficients to get values of spline at each X
        RealFit<- Fit+coef(tempMod)[1]+NTmean+MPhmean # add intercept and other terms at their mean values
        # RealFit<- inv.logit(RealFit)
        BootstrapFits<- (JDBasis%*%t(JDayBootstrapCoefs))+coef(tempMod)[1]+NTmean+MPhmean # get spread of spline values at each X based on distributions of each coefficient
        quant.func<- function(x){quantile(x, probs=c(0.0275,0.975))}
        cis<-apply(BootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
        # cil<-inv.logit(cis[1,]) # lowerCI bound
        # ciu<-inv.logit(cis[2,]) # upper CI bound
        cil<-cis[1,] # lowerCI bound
        ciu<-cis[2,] # upper CI bound
        
        plotDF = data.frame(JDayForPlotting,RealFit)
        colnames(plotDF) = c("Jday","Fit")
        
        PresStats[[CTname]]$JDPeak[siteInd] = JDayForPlotting[which(RealFit==max(RealFit))]
        
        JD = ggplot(plotDF, aes(Jday, Fit),
        ) + geom_smooth(aes(ymin=cil, ymax=ciu),
                        color="#16A7CA",
                        fill="#16A7CA",
                        alpha=0.2,
                        stat ="identity"
        ) + labs(x = "Julian Day",
                 # y = "Probability"
                 y = "s(JulianDay)"
        ) + scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),
                               label=c("J","F","M","A","M","J","J","A","S","O","N","D")
        ) + theme(axis.line = element_line(size=0.2),
                  panel.background = element_blank()
        )
        
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDPlot.png",sep="")
        ggsave(saveName,device="png", width=2, scale=3, height=1, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDPlot.pdf",sep="")
        ggsave(saveName,device="pdf", width=2, scale=3, height=1, units="in",dpi=600)
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
        ggsave(saveName,device="png", width=2, scale=4, height=1, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDDataDensity.pdf",sep="")
        ggsave(saveName,device="pdf", width=2, scale=4, height=1, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
      }
      
      # Moon Phase ---------------------------
      if (!isempty(MPhInd) & isempty(NTMPhInd)){
        Fit<- MPhBasis%*%(coef(tempMod)[MPhInd]) # multiply basis functions by model coefficients to get values of spline at each X
        RealFit<- Fit+coef(tempMod)[1]+JDmean+NTmean # add intercept and other terms at their mean values
        # RealFit<- inv.logit(RealFit)
        BootstrapFits<- (MPhBasis%*%t(MPhBootstrapCoefs))+coef(tempMod)[1]+JDmean+NTmean # get spread of spline values at each X based on distributions of each coefficient
        quant.func<- function(x){quantile(x, probs=c(0.0275,0.975))}
        cis<-apply(BootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
        # cil<-inv.logit(cis[1,]) # lowerCI bound
        # ciu<-inv.logit(cis[2,]) # upper CI bound
        cil<-cis[1,] # lowerCI bound
        ciu<-cis[2,] # upper CI bound
        
        plotDF = data.frame(MPhForPlotting,RealFit)
        colnames(plotDF) = c("MoonPhase","Fit")
        
        MPh = ggplot(plotDF, aes(MoonPhase, Fit),
        ) + geom_smooth(aes(ymin=cil, ymax=ciu),
                        color="#16A7CA",
                        fill="#16A7CA",
                        alpha=0.2,
                        stat ="identity"
        ) + labs(x = "MoonPhase",
                 # y = "Probability"
                 y = "s(MoonPhase)"
        ) + scale_x_continuous(breaks=c(0,0.5,1),
                               labels=c("New Moon","Full Moon","New Moon")
        ) + theme(axis.line = element_line(size=0.2),
                  panel.background = element_blank()
        )
        
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_MPhPlot.png",sep="")
        ggsave(saveName,device="png", width=2, scale=3, height=1, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_MPhPlot.pdf",sep="")
        ggsave(saveName,device="pdf", width=2, scale=3, height=1, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
        
        MPhdens = ggplot(thisSite,aes(x=MoonPhase)
        )+geom_histogram(aes(y=..ncount..),
                         fill='#66B2FF',
                         binwidth = 0.01,
                         alpha=0.5
        ) + scale_x_continuous(breaks=c(0,0.5,1),
                               labels=c("New Moon","Full Moon","New Moon")
        ) + labs(y="Count",x=NULL
        ) + theme_minimal()
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_MphDataDensity.png",sep="")
        ggsave(saveName,device="png", width=2, scale=4, height=1, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_MPhDataDensity.pdf",sep="")
        ggsave(saveName,device="pdf", width=2, scale=4, height=1, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
      }
      
      # Norm Time ---------------------------
      if (!isempty(NTInd) & isempty(JDNTInd) & isempty(NTMPhInd)){
        Fit<- NTBasis%*%(coef(tempMod)[NTInd]) # multiply basis functions by model coefficients to get values of spline at each X
        RealFit<- Fit+coef(tempMod)[1]+JDmean+MPhmean # add intercept and other terms at their mean values
        # RealFit<- inv.logit(RealFit)
        BootstrapFits<- (NTBasis%*%t(NTBootstrapCoefs))+coef(tempMod)[1]+JDmean+MPhmean # get spread of spline values at each X based on distributions of each coefficient
        quant.func<- function(x){quantile(x, probs=c(0.0275,0.975))}
        cis<-apply(BootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
        # cil<-inv.logit(cis[1,]) # lowerCI bound
        # ciu<-inv.logit(cis[2,]) # upper CI bound
        cil<-cis[1,] # lowerCI bound
        ciu<-cis[2,] # upper CI bound
        
        plotDF = data.frame(NTForPlotting,RealFit)
        colnames(plotDF) = c("NormTime","Fit")
        
        # PresStats[[CTname]]$DielPeak[siteInd] = NTForPlotting[which(RealFit==max(RealFit))]
        
        NT = ggplot(plotDF, aes(NormTime, Fit),
        ) + geom_smooth(aes(ymin=cil, ymax=ciu),
                        color="#16A7CA",
                        fill="#16A7CA",
                        alpha=0.2,
                        stat ="identity"
        ) + labs(x = NULL,
                 # y = "Probability"
                 y = "s(NormTime)"
        ) + coord_cartesian(xlim=c(-1,1)
        ) + scale_x_continuous(breaks=c(-1,0,1),
                               labels=c("Sunrise","Sunset","Sunrise")
        ) + theme(axis.line = element_line(size=0.2),
                  panel.background = element_blank()
        )
        
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NTPlot.png",sep="")
        ggsave(saveName,device="png", width=2, scale=3, height=1, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NTPlot.pdf",sep="")
        ggsave(saveName,device="pdf", width=2, scale=3, height=1, units="in",dpi=600)
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
        ggsave(saveName,device="png", width=1, scale=4, height=1, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NTDataDensity.pdf",sep="")
        ggsave(saveName,device="pdf", width=1, scale=4, height=1, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
      }
      
      # JD:NT Interaction ------------------------------
      if (!isempty(JDNTInd)){
        FitJD<- JDBasis%*%(coef(tempMod)[JDInd]) # multiply basis functions by model coefficients to get values of spline at each X
        FitNT<- NTBasis%*%(coef(tempMod)[NTInd]) # multiply basis functions by model coefficients to get values of spline at each X
        
        ### Plot NT smooth (if not also modified by MPh) ------------------
        if (isempty(NTMPhInd)){
          ## NT at different values of JD
          IntPlotDF = numeric()
          IntCI = numeric()
          testRows = c(1,250,500,750)
          for (l in 1:length(testRows)){
            # Recreate interaction basis
            IntBasis = numeric()
            for (k in 1:3){
              IntBasis = cbind(IntBasis,NTBasis[,1]*JDBasis[testRows[l],k],NTBasis[,2]*JDBasis[testRows[l],k],NTBasis[,3]*JDBasis[testRows[l],k])
            }
            
            # Calculate NT smooth at several different values of JD
            IntAdjust = IntBasis%*%coef(tempMod)[JDNTInd]
            RealFitInt<- FitNT+coef(tempMod)[1]+IntAdjust+MPhmean # add intercept, interaction term, and other terms at their mean values
            # RealFitInt<- inv.logit(RealFitInt)
            IntBootstrapFits<- (NTBasis%*%t(NTBootstrapCoefs))+coef(tempMod)[1]+(IntBasis%*%t(JDNTBootstrapCoefs))+MPhmean # get spread of spline values at each X based on distributions of each coefficient
            quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
            cisInt<-apply(IntBootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
            # IntCI = rbind(IntCI,cbind(NTForPlotting,inv.logit(t(cisInt)),rep(l,times=1000))) # get lower and upper CI bounds
            IntCI = rbind(IntCI,cbind(NTForPlotting,t(cisInt),rep(l,times=1000)))
            
            IntPlotDF = rbind(IntPlotDF,cbind(NTForPlotting,RealFitInt,rep(l,times=1000)))
            
          }
          
          IntPlotDF = as.data.frame(IntPlotDF)
          colnames(IntPlotDF) = c("NormTime","Data","Fit")
          IntCI = as.data.frame(IntCI)
          colnames(IntCI) = c("NormTime","Ymin","Ymax","Fit2")
          IntPlotDF$Fit = as.factor(IntPlotDF$Fit)
          IntCI$Fit2 = as.factor(IntCI$Fit2)
          
          NTJD = ggplot(
          ) + geom_line(data=IntPlotDF,aes(x=NormTime,y=Data,group=Fit,color=Fit),size=0.5
          ) + scale_color_manual("Fit",values=c("#33FFE0","#335CFF","#D933FF","#FF334B"),
                                 labels=c("Winter","Spring","Summer","Fall"),
                                 name=NULL
          ) + geom_ribbon(data=IntCI,aes(x=NormTime,ymin=Ymin, ymax=Ymax,group=Fit2,fill=Fit2),
                          alpha=0.2,
                          stat ="identity"
          ) + scale_fill_manual("Fit2",values=c("#33FFE0","#335CFF","#D933FF","#FF334B"),
                                guide=NULL
          ) + labs(x = NULL,
                   # y = "Probability"
                   y = "s(NormTime)"
          ) + coord_cartesian(xlim=c(-1,1)
          ) + scale_x_continuous(breaks=c(-1,0,1),
                                 labels=c("Sunrise","Sunset","Sunrise")
          ) + theme(legend.position="right",
                    axis.line = element_line(size=0.2),
                    panel.background = element_blank()
          )
          
          saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NTJDPlot_Partial.png",sep="")
          ggsave(saveName,device="png", width=2, scale=4, height=1, units="in",dpi=600)
          while (dev.cur()>1) {dev.off()}
          saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NTJDPlot_Partial.pdf",sep="")
          ggsave(saveName,device="pdf", width=2, scale=4, height=1, units="in",dpi=600)
          while (dev.cur()>1) {dev.off()}
          
          
          ## NT across values of JD
          IntPlotDF = numeric()
          IntCI = numeric()
          # Recreate interaction basis
          rowMat = numeric()
          IntBasis = numeric()
          for (l in 1:dim(NTBasis)[1]){ # for each row in NT basis
            for (k in 1:3){ # multiply by all values of JD
              rowMat = cbind(rowMat,NTBasis[l,1]*JDBasis[,k],NTBasis[l,2]*JDBasis[,k],NTBasis[l,3]*JDBasis[,k])
            }
            IntBasis = rbind(IntBasis,apply(rowMat,MARGIN=2,mean)) # average this NT across JD
            rowMat = numeric()
          }
          
          # Calculate NT smooth averaged across JD
          IntAdjust = IntBasis%*%coef(tempMod)[JDNTInd]
          RealFitInt<- FitNT+coef(tempMod)[1]+IntAdjust+MPhmean # add intercept, interaction term, and other terms at their mean values
          # RealFitInt<- inv.logit(RealFitInt)
          IntBootstrapFits<- (NTBasis%*%t(NTBootstrapCoefs))+coef(tempMod)[1]+(IntBasis%*%t(JDNTBootstrapCoefs))+MPhmean # get spread of spline values at each X based on distributions of each coefficient
          quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
          cisInt<-apply(IntBootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
          # IntCI = cbind(NTForPlotting,inv.logit(t(cisInt))) # get lower and upper CI bounds
          IntCI = cbind(NTForPlotting,t(cisInt)) # get lower and upper CI bounds
          
          
          IntPlotDF = cbind(NTForPlotting,RealFitInt)
          
          IntPlotDF = as.data.frame(IntPlotDF)
          colnames(IntPlotDF) = c("NormTime","Data")
          IntCI = as.data.frame(IntCI)
          colnames(IntCI) = c("NormTime","Ymin","Ymax")
          
          NTJDInt = ggplot(
          ) + geom_line(data=IntPlotDF,aes(x=NormTime,y=Data),size=0.5,color="51D4A5"
          ) + geom_ribbon(data=IntCI,aes(x=NormTime,ymin=Ymin, ymax=Ymax),
                          fill="51D4A5",
                          alpha=0.2,
                          stat ="identity"
          ) + labs(x = NULL,
                   # y = "Probability"
                   y = "s(NormTime)"
          ) + coord_cartesian(xlim=c(-1,1)
          ) + scale_x_continuous(breaks=c(-1,0,1),
                                 labels=c("Sunrise","Sunset","Sunrise")
          ) + theme(legend.position="right",
                    axis.line = element_line(size=0.2),
                    panel.background = element_blank()
          )
          
          saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NTJDPlot_Overall.png",sep="")
          ggsave(saveName,device="png", width=2, scale=4, height=1, units="in",dpi=600)
          while (dev.cur()>1) {dev.off()}
          saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NTJDPlot_Overall.pdf",sep="")
          ggsave(saveName,device="pdf", width=2, scale=4, height=1, units="in",dpi=600)
          while (dev.cur()>1) {dev.off()}
        }
        
        
        ### Plot JD smooth ------------------------------
        # if (as.numeric(PV$'p-value'[PV$'Variable'=="JDs"])<0.05){
          # ## JD at different values of NT
          # IntPlotDF = numeric()
          # IntCI = numeric()
          # testRows = c(1,250,500,750)
          # for (l in 1:length(testRows)){
          #   # Recreate interaction basis
          #   IntBasis = numeric()
          #   for (k in 1:3){
          #     IntBasis = cbind(IntBasis,NTBasis[testRows[l],1]*JDBasis[,k],NTBasis[testRows[l],2]*JDBasis[,k],NTBasis[testRows[l],3]*JDBasis[,k])
          #   }
          #   
          #   # Calculate JD smooth at several different values of NT
          #   IntAdjust = IntBasis%*%coef(tempMod)[JDNTInd]
          #   RealFitInt<- FitJD+coef(tempMod)[1]+IntAdjust+MPhmean # add intercept, interaction term, and other terms at their mean values
          #   RealFitInt<- inv.logit(RealFitInt)
          #   IntBootstrapFits<- (JDBasis%*%t(JDayBootstrapCoefs))+coef(tempMod)[1]+(IntBasis%*%t(JDNTBootstrapCoefs))+MPhmean # get spread of spline values at each X based on distributions of each coefficient
          #   quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
          #   cisInt<-apply(IntBootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
          #   IntCI = rbind(IntCI,cbind(JDayForPlotting,inv.logit(t(cisInt)),rep(l,times=1000))) # get lower and upper CI bounds
          #   
          #   IntPlotDF = rbind(IntPlotDF,cbind(JDayForPlotting,RealFitInt,rep(l,times=1000)))
          #   
          # }
          # 
          # IntPlotDF = as.data.frame(IntPlotDF)
          # colnames(IntPlotDF) = c("JulianDay","Data","Fit")
          # IntCI = as.data.frame(IntCI)
          # colnames(IntCI) = c("JulianDay","Ymin","Ymax","Fit2")
          # IntPlotDF$Fit = as.factor(IntPlotDF$Fit)
          # IntCI$Fit2 = as.factor(IntCI$Fit2)
          # 
          # JDNTInt = ggplot(
          # ) + geom_line(data=IntPlotDF,aes(x=JulianDay,y=Data,group=Fit,color=Fit),size=1
          # ) + scale_color_manual("Fit",values=c("#33FFE0","#335CFF","#D933FF","#FF334B"),
          #                        labels=c("Sunrise","Midday","Sunset","Midnight"),
          #                        name=NULL
          # ) + geom_ribbon(data=IntCI,aes(x=JulianDay,ymin=Ymin, ymax=1.05*Ymax,group=Fit2,fill=Fit2),
          #                 alpha=0.2,
          #                 stat ="identity"
          # ) + scale_fill_manual("Fit2",values=c("#33FFE0","#335CFF","#D933FF","#FF334B"),
          #                       guide=NULL
          # ) + labs(x = "Julian Day",
          #          y = "Probability"
          # ) + coord_cartesian(xlim=c(1,365)
          # ) + scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),
          #                        label=c("J","F","M","A","M","J","J","A","S","O","N","D")
          # ) + theme(legend.position="right",
          #           axis.line = element_line(size=0.2),
          #           panel.background = element_blank()
          # )
          # 
          # saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDNTPlot_Partial.png",sep="")
          # ggsave(saveName,device="png", width=2, scale=4, height=1, units="in",dpi=600)
          # while (dev.cur()>1) {dev.off()}
          # saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDNTPlot_Partial.pdf",sep="")
          # ggsave(saveName,device="pdf", width=2, scale=4, height=1, units="in",dpi=600)
          # while (dev.cur()>1) {dev.off()}
          
          ## JD across value of NT
          IntPlotDF = numeric()
          IntCI = numeric()
          # Recreate interaction basis
          rowMat = numeric()
          IntBasis = numeric()
          for (l in 1:dim(JDBasis)[1]){ # for each row in JD basis
            for (k in 1:3){ # multiply by all values of NT
              rowMat = cbind(rowMat,NTBasis[,1]*JDBasis[l,k],NTBasis[,2]*JDBasis[l,k],NTBasis[,3]*JDBasis[l,k])
            }
            IntBasis = rbind(IntBasis,apply(rowMat,MARGIN=2,mean)) # average this JD across NT
            rowMat = numeric()
          }
          
          # Calculate JD smooth averaged across NT values
          IntAdjust = IntBasis%*%coef(tempMod)[JDNTInd]
          RealFitInt<- FitJD+coef(tempMod)[1]+IntAdjust+MPhmean # add intercept, interaction term, and other terms at their mean values
          # RealFitInt<- inv.logit(RealFitInt)
          IntBootstrapFits<- (JDBasis%*%t(JDayBootstrapCoefs))+coef(tempMod)[1]+(IntBasis%*%t(JDNTBootstrapCoefs))+MPhmean # get spread of spline values at each X based on distributions of each coefficient
          quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
          cisInt<-apply(IntBootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
          # IntCI = cbind(JDayForPlotting,inv.logit(t(cisInt))) # get lower and upper CI bounds
          IntCI = cbind(JDayForPlotting,t(cisInt)) # get lower and upper CI bounds
          
          PresStats[[CTname]]$JDPeak[siteInd] = JDayForPlotting[which(RealFitInt==max(RealFitInt))]
          
          IntPlotDF = cbind(JDayForPlotting,RealFitInt)
          
          IntPlotDF = as.data.frame(IntPlotDF)
          colnames(IntPlotDF) = c("JulianDay","Data")
          IntCI = as.data.frame(IntCI)
          colnames(IntCI) = c("JulianDay","Ymin","Ymax")
          
          JDNTInt = ggplot(
          ) + geom_line(data=IntPlotDF,aes(x=JulianDay,y=Data,),size=0.5,color="51D4A5"
          ) + geom_ribbon(data=IntCI,aes(x=JulianDay,ymin=Ymin, ymax=Ymax),
                          fill="51D4A5",
                          alpha=0.2,
                          stat ="identity"
          ) + labs(x = "Julian Day",
                   y = "s(JulianDay)"
          ) + coord_cartesian(xlim=c(1,365)
          ) + scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),
                                 label=c("J","F","M","A","M","J","J","A","S","O","N","D")
          ) + theme(legend.position="right",
                    axis.line = element_line(size=0.2),
                    panel.background = element_blank()
          )
          
          saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDNTPlot_Overall.png",sep="")
          ggsave(saveName,device="png", width=2, scale=4, height=1, units="in",dpi=600)
          while (dev.cur()>1) {dev.off()}
          saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_JDNTPlot_Overall.pdf",sep="")
          ggsave(saveName,device="pdf", width=2, scale=4, height=1, units="in",dpi=600)
          while (dev.cur()>1) {dev.off()}
        # }
          
          # # calculate NT smooth at JD of peak predicted presence (based on JD averaged across all NT)
          # peakInd = which(FitJD==max(FitJD))
          # # Recreate interaction basis
          # IntBasis = numeric()
          # for (k in 1:3){
          #   IntBasis = cbind(IntBasis,NTBasis[,1]*JDBasis[peakInd,k],NTBasis[,2]*JDBasis[peakInd,k],NTBasis[,3]*JDBasis[peakInd,k])
          # }
          # IntAdjust = IntBasis%*%coef(tempMod)[JDNTInd]
          # RealFitInt<- FitNT+coef(tempMod)[1]+IntAdjust+MPhmean # add intercept, interaction term, and other terms at their mean values
          # PresStats[[CTname]]$DielPeak[siteInd] = NTForPlotting[which(RealFitInt==max(RealFitInt))]
        
      }
      
      
      # NT:MPh Interaction ------------------------------
      if (!isempty(NTMPhInd)){

        FitMPh<- MPhBasis%*%(coef(tempMod)[MPhInd]) # multiply basis functions by model coefficients to get values of spline at each X
        FitNT<- NTBasis%*%(coef(tempMod)[NTInd]) # multiply basis functions by model coefficients to get values of spline at each X

        ### Plot NT smooth (if not further modified by JD) ------------------
        if (isempty(JDNTInd)){
          ## NT at different values of MPh
          IntPlotDF = numeric()
          IntCI = numeric()
          testRows = c(1,250,500,750)
          for (l in 1:length(testRows)){
            # Recreate interaction basis
            IntBasis = numeric()
            for (k in 1:3){
              IntBasis = cbind(IntBasis,NTBasis[,1]*MPhBasis[testRows[l],k],NTBasis[,2]*MPhBasis[testRows[l],k],NTBasis[,3]*MPhBasis[testRows[l],k])
            }

            # Calculate NT smooth at several different values of MPh
            IntAdjust = IntBasis%*%coef(tempMod)[NTMPhInd]
            RealFitInt<- FitNT+coef(tempMod)[1]+IntAdjust+JDmean # adjust offset
            # RealFitInt<- inv.logit(RealFitInt)
            IntBootstrapFits<- (NTBasis%*%t(NTBootstrapCoefs))+coef(tempMod)[1]+(IntBasis%*%t(NTMPhBootstrapCoefs))+JDmean # get spread of spline values at each X based on distributions of each coefficient
            quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
            cisInt<-apply(IntBootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
            # IntCI = rbind(IntCI,cbind(NTForPlotting,inv.logit(t(cisInt)),rep(l,times=1000))) # get lower and upper CI bounds
            IntCI = rbind(IntCI,cbind(NTForPlotting,t(cisInt),rep(l,times=1000))) # get lower and upper CI bounds
            

            IntPlotDF = rbind(IntPlotDF,cbind(NTForPlotting,RealFitInt,rep(l,times=1000)))

          }

          IntPlotDF = as.data.frame(IntPlotDF)
          colnames(IntPlotDF) = c("NormTime","Data","Fit")
          IntCI = as.data.frame(IntCI)
          colnames(IntCI) = c("NormTime","Ymin","Ymax","Fit2")
          IntPlotDF$Fit = as.factor(IntPlotDF$Fit)
          IntCI$Fit2 = as.factor(IntCI$Fit2)

          NTMPh = ggplot(
          ) + geom_line(data=IntPlotDF,aes(x=NormTime,y=Data,group=Fit,color=Fit),size=0.5
          ) + scale_color_manual("Fit",values=c("#33FFE0","#335CFF","#D933FF","#FF334B"),
                                 labels=c("New Moon","First Quarter","Full Moon","Third Quarter"),
                                 name=NULL
          ) + geom_ribbon(data=IntCI,aes(x=NormTime,ymin=Ymin, ymax=Ymax,group=Fit2,fill=Fit2),
                          alpha=0.2,
                          stat ="identity"
          ) + scale_fill_manual("Fit2",values=c("#33FFE0","#335CFF","#D933FF","#FF334B"),
                                guide=NULL
          ) + labs(x = NULL,
                   # y = "Probability"
                   y = "s(NormTime)"
          ) + coord_cartesian(xlim=c(-1,1)
          ) + scale_x_continuous(breaks=c(-1,0,1),
                                 labels=c("Sunrise","Sunset","Sunrise")
          ) + theme(legend.position="right",
                    axis.line = element_line(size=0.2),
                    panel.background = element_blank()
          )

          saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NTMPhPlot_Partial.png",sep="")
          ggsave(saveName,device="png", width=2, scale=5, height=1, units="in",dpi=600)
          while (dev.cur()>1) {dev.off()}
          saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NTMPhPlot_Partial.pdf",sep="")
          ggsave(saveName,device="pdf", width=2, scale=5, height=1, units="in",dpi=600)
          while (dev.cur()>1) {dev.off()}


          ## NT across values of MPh
          IntPlotDF = numeric()
          IntCI = numeric()
          # Recreate interaction basis
          rowMat = numeric()
          IntBasis = numeric()
          for (l in 1:dim(NTBasis)[1]){ # for each row in NT basis
            for (k in 1:3){ # multiply by all values of MPh
              rowMat = cbind(rowMat,NTBasis[l,1]*MPhBasis[,k],NTBasis[l,2]*MPhBasis[,k],NTBasis[l,3]*MPhBasis[,k])
            }
            IntBasis = rbind(IntBasis,apply(rowMat,MARGIN=2,mean)) # average this NT across MPh
            rowMat = numeric()
          }

          # Calculate NT smooth averaged across MPh
          IntAdjust = IntBasis%*%coef(tempMod)[NTMPhInd]
          RealFitInt<- FitNT+coef(tempMod)[1]+IntAdjust+JDmean # adjust offset
          # RealFitInt<- inv.logit(RealFitInt)
          IntBootstrapFits<- (NTBasis%*%t(NTBootstrapCoefs))+coef(tempMod)[1]+(IntBasis%*%t(NTMPhBootstrapCoefs))+JDmean # get spread of spline values at each X based on distributions of each coefficient
          quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
          cisInt<-apply(IntBootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
          # IntCI = cbind(NTForPlotting,inv.logit(t(cisInt))) # get lower and upper CI bounds
          IntCI = cbind(NTForPlotting,t(cisInt)) # get lower and upper CI bounds
          

          IntPlotDF = cbind(NTForPlotting,RealFitInt)

          IntPlotDF = as.data.frame(IntPlotDF)
          colnames(IntPlotDF) = c("NormTime","Data")
          IntCI = as.data.frame(IntCI)
          colnames(IntCI) = c("NormTime","Ymin","Ymax")

          NTMPh = ggplot(
          ) + geom_line(data=IntPlotDF,aes(x=NormTime,y=Data),size=0.5,color="51D4A5"
          ) + geom_ribbon(data=IntCI,aes(x=NormTime,ymin=Ymin, ymax=Ymax),
                          fill="51D4A5",
                          alpha=0.2,
                          stat ="identity"
          ) + labs(x = NULL,
                   # y = "Probability"
                   y = "s(NormTime)"
          ) + coord_cartesian(xlim=c(-1,1)
          ) + scale_x_continuous(breaks=c(-1,0,1),
                                 labels=c("Sunrise","Sunset","Sunrise")
          ) + theme(legend.position="right",
                    axis.line = element_line(size=0.2),
                    panel.background = element_blank()
          )

          saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NTMPhPlot_Overall.png",sep="")
          ggsave(saveName,device="png", width=2, scale=5, height=1, units="in",dpi=600)
          while (dev.cur()>1) {dev.off()}
          saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NTMPhPlot_Overall.pdf",sep="")
          ggsave(saveName,device="pdf", width=2, scale=5, height=1, units="in",dpi=600)
          while (dev.cur()>1) {dev.off()}
        }

        # # calculate NT smooth at MPh of peak predicted presence
        # peakInd = which(FitMPh==max(FitMPh))
        # # Recreate interaction basis
        # IntBasis = numeric()
        # for (k in 1:3){
        #   IntBasis = cbind(IntBasis,NTBasis[,1]*JDBasis[peakInd,k],NTBasis[,2]*JDBasis[peakInd,k],NTBasis[,3]*JDBasis[peakInd,k])
        # }
        # IntAdjust = IntBasis%*%coef(tempMod)[JDNTInd]
        # RealFitInt<- FitNT+coef(tempMod)[1]+IntAdjust+MPhmean # add intercept, interaction term, and other terms at their mean values
        # PresStats[[CTname]]$DielPeak[siteInd] = NTForPlotting[which(RealFitInt==max(RealFitInt))]
        
        ### Plot MPh smooth  ------------------------------
        # if (as.numeric(PV$'p-value'[PV$'Variable'=="MPhs"])<0.05){
          ## MPh at during daytime and nighttime
          IntPlotDF = numeric()
          IntCI = numeric()
          dayRows = which(NTForPlotting<0)
          nightRows = which(NTForPlotting>0)
          testRows = list(dayRows,nightRows)

          for (l in 1:length(testRows)){
            IntBasis = numeric()
            # Recreate interaction basis
            for (m in 1:dim(MPhBasis)[1]){ # for each row in MPh basis
              rowMat = numeric()
              for (k in 1:3){ # multiply by all daytime OR nighttime values of NT
                rowMat = cbind(rowMat,NTBasis[unlist(testRows[l]),1]*MPhBasis[m,k],NTBasis[unlist(testRows[l]),2]*MPhBasis[m,k],NTBasis[unlist(testRows[l]),3]*MPhBasis[m,k])
              }
              IntBasis = rbind(IntBasis,apply(rowMat,MARGIN=2,mean)) # average this MPh across NT
            }

            # Calculate MPh smooth in daytime/nighttime
            IntAdjust = IntBasis%*%coef(tempMod)[NTMPhInd]
            RealFitInt<- FitMPh+coef(tempMod)[1]+IntAdjust+JDmean # adjust offset
            # RealFitInt<- inv.logit(RealFitInt)
            IntBootstrapFits<- (MPhBasis%*%t(MPhBootstrapCoefs))+coef(tempMod)[1]+(IntBasis%*%t(NTMPhBootstrapCoefs))+JDmean # get spread of spline values at each X based on distributions of each coefficient
            quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
            cisInt<-apply(IntBootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
            # IntCI = rbind(IntCI,cbind(MPhForPlotting,inv.logit(t(cisInt)),rep(l,times=1000))) # get lower and upper CI bounds
            IntCI = rbind(IntCI,cbind(MPhForPlotting,t(cisInt),rep(l,times=1000))) # get lower and upper CI bounds
            
            IntPlotDF = rbind(IntPlotDF,cbind(MPhForPlotting,RealFitInt,rep(l,times=1000)))

          }

          IntPlotDF = as.data.frame(IntPlotDF)
          colnames(IntPlotDF) = c("MoonPhase","Data","Fit")
          IntCI = as.data.frame(IntCI)
          colnames(IntCI) = c("MoonPhase","Ymin","Ymax","Fit2")
          IntPlotDF$Fit = as.factor(IntPlotDF$Fit)
          IntCI$Fit2 = as.factor(IntCI$Fit2)

          MPhNT = ggplot(
          ) + geom_line(data=IntPlotDF,aes(x=MoonPhase,y=Data,group=Fit,color=Fit),size=0.5
          ) + scale_color_manual("Fit",values=c("#335CFF","#D933FF"),
                                 labels=c("Day","Night"),
                                 name=NULL
          ) + geom_ribbon(data=IntCI,aes(x=MoonPhase,ymin=Ymin, ymax=Ymax,group=Fit2,fill=Fit2),
                          alpha=0.2,
                          stat ="identity"
          ) + scale_fill_manual("Fit2",values=c("#335CFF","#D933FF"),
                                guide=NULL
          ) + labs(x = NULL,
                   # y = "Probability"
                   y = "s(MoonPhase)"
          ) + scale_x_continuous(breaks=c(0,0.5,1),
                                 label=c("New Moon","Full Moon","New Moon")
          ) + theme(legend.position="right",
                    axis.line = element_line(size=0.2),
                    panel.background = element_blank()
          )

          saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_MPhNTPlot_Partial.png",sep="")
          ggsave(saveName,device="png", width=2, scale=5, height=1, units="in",dpi=600)
          while (dev.cur()>1) {dev.off()}
          saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_MPhNTPlot_Partial.pdf",sep="")
          ggsave(saveName,device="pdf", width=2, scale=5, height=1, units="in",dpi=600)
          while (dev.cur()>1) {dev.off()}

        # }

      }

      # NT:JD & NT:MPh Interaction -----------------------------
      if (!isempty(JDNTInd) & !isempty(NTMPhInd)){
        # Plot NT as modified by a few different values of JD and MPh
        
        FitNT<- NTBasis%*%(coef(tempMod)[NTInd]) # multiply basis functions by model coefficients to get values of spline at each X
        
        Titles = c("Winter","Spring","Summer","Fall")
        ylabs = c("s(NormTime)","","","")
        plotList = list(Winter=list(),Spring=list(),Summer=list(),Fall=list(),Legend=list())
        
        JDIntAdjust = numeric()
        MPhIntAdjust = numeric()
        JDtestRows = c(1,250,500,750)
        MPhtestRows = c(1,250,500,750)
        
        extract_legend <- function(my_ggp) {
          step1 <- ggplot_gtable(ggplot_build(my_ggp))
          step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
          step3 <- step1$grobs[[step2]]
          return(step3)
        }
        
        for (l in 1:length(JDtestRows)){
          # Reconstruct basis functions for each combination of NT and JD and NT and MPh of interest
          JDIntBasis = numeric()
          for (k in 1:3){
            JDIntBasis = cbind(JDIntBasis,NTBasis[,1]*JDBasis[JDtestRows[l],k],NTBasis[,2]*JDBasis[JDtestRows[l],k],NTBasis[,3]*JDBasis[JDtestRows[l],k])
          }
          JDIntAdjust = JDIntBasis%*%coef(tempMod)[JDNTInd]
          
          IntPlotDF = numeric()
          IntCI = numeric()
          for (m in 1:length(MPhtestRows)){
            MPhIntBasis = numeric()
            for (k in 1:3){
            MPhIntBasis = cbind(MPhIntBasis,NTBasis[,1]*MPhBasis[MPhtestRows[m],k],NTBasis[,2]*MPhBasis[MPhtestRows[m],k],NTBasis[,3]*MPhBasis[MPhtestRows[l],k])
            }
            MPhIntAdjust = MPhIntBasis%*%coef(tempMod)[NTMPhInd]
            RealFitInt<- FitNT+coef(tempMod)[1]+JDIntAdjust+MPhIntAdjust
            # RealFitInt<- inv.logit(RealFitInt)
            IntPlotDF = rbind(IntPlotDF,cbind(NTForPlotting,RealFitInt,rep(m,times=1000)))
            
            IntBootstrapFits<- (NTBasis%*%t(NTBootstrapCoefs))+coef(tempMod)[1]+(JDIntBasis%*%t(JDNTBootstrapCoefs))+(MPhIntBasis%*%t(NTMPhBootstrapCoefs)) # get spread of spline values at each X based on distributions of each coefficient
            quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
            cisInt<-apply(IntBootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
            # IntCI = rbind(IntCI,cbind(NTForPlotting,inv.logit(t(cisInt)),rep(m,times=1000))) # get lower and upper CI bounds
            IntCI = rbind(IntCI,cbind(NTForPlotting,t(cisInt),rep(m,times=1000)))
            }
          
          IntPlotDF = as.data.frame(IntPlotDF)
          colnames(IntPlotDF) = c("NormTime","Data","Fit")
          IntCI = as.data.frame(IntCI)
          colnames(IntCI) = c("NormTime","Ymin","Ymax","Fit2")
          IntPlotDF$Fit = as.factor(IntPlotDF$Fit)
          IntCI$Fit2 = as.factor(IntCI$Fit2)
          
          plotList[[Titles[l]]] = ggplot(
          ) + geom_line(data=IntPlotDF,aes(x=NormTime,y=Data,group=Fit,color=Fit),size=0.5
          ) + scale_color_manual("Fit",values=c("#33FFE0","#335CFF","#D933FF","#FF334B"),
                                 # labels=c("New Moon","First Quarter","Full Moon","Third Quarter"),
                                 # name=NULL
                                 guide=NULL
          ) + geom_ribbon(data=IntCI,aes(x=NormTime,ymin=Ymin, ymax=Ymax,group=Fit2,fill=Fit2),
                          alpha=0.2,
                          stat ="identity"
          ) + scale_fill_manual("Fit2",values=c("#33FFE0","#335CFF","#D933FF","#FF334B"),
                                guide=NULL
          ) + labs(x = NULL,
                   y = ylabs[l],
                   title = Titles[l]
          ) + coord_cartesian(xlim=c(-1,1)
          ) + scale_x_continuous(breaks=c(-1,0,1),
                                 labels=c("Sunrise","Sunset","Sunrise")
          ) + theme(legend.position="right",
                    axis.line = element_line(size=0.2),
                    panel.background = element_blank()
          )
          
          if (l==4){ # grab legend to later plot once alongside multiple panels
          legendPlot = ggplot(
          ) + geom_line(data=IntPlotDF,aes(x=NormTime,y=Data,group=Fit,color=Fit),size=0.5
          ) + scale_color_manual("Fit",values=c("#33FFE0","#335CFF","#D933FF","#FF334B"),
                                 labels=c("New Moon","First Quarter","Full Moon","Third Quarter"),
                                 name=NULL)
          plotList[[l+1]] = extract_legend(legendPlot)}
          
        }
        
        png(file=paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NTJDMPhPlot_Partial.png",sep=""),width = 1800, height = 400, units = "px",res=125)
        grid.arrange(grobs=plotList,ncol=9,nrow=1,layout_matrix=rbind(c(1,1,2,2,3,3,4,4,5)),top=paste(CTname,'at',sites[j]))
        while (dev.cur()>1) {dev.off()}
        
        g=arrangeGrob(grobs=plotList,ncol=9,nrow=1,layout_matrix=rbind(c(1,1,2,2,3,3,4,4,5)),top=paste(CTname,'at',sites[j]))
        savename=paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_NTJDMPhPlot_Partial.pdf",sep="")
        ggsave(savename,g,device="pdf", width=1800, scale=5, height=400, units="px",dpi=600)
        while (dev.cur()>1) {dev.off()}
      }
      
      # Year (as boxplot) -----------------------------------
      if (!isempty(YrInd)){
        if (sites[j]=="HAT"){ # only 2 years of data from HATB
          AdjustedYearCoefs = data.frame(c(YearBootstrapCoefs[,1]-mean(YearBootstrapCoefs[,1]),
                                           YearBootstrapCoefs[,2]),
                                         as.factor(rep(1:2,each=10000)))
        } else {
          AdjustedYearCoefs = data.frame(c(YearBootstrapCoefs[,1]-mean(YearBootstrapCoefs[,1]),
                                           YearBootstrapCoefs[,2],
                                           YearBootstrapCoefs[,3]),
                                         as.factor(rep(1:3,each=10000)))}
        colnames(AdjustedYearCoefs) = c("Coefficient","Year")
        # AdjustedYearCoefs$Prob = inv.logit(AdjustedYearCoefs$Coefficient)
        
        # calculate quantiles for plotting limits
        quants = AdjustedYearCoefs %>%
          group_by(Year) %>%
          summarize(q25 = quantile(Coefficient,probs=0.25),
                    q75 = quantile(Coefficient,probs=0.75))
        iqr = quants$q75-quants$q25
        
        Yr = ggplot(AdjustedYearCoefs,aes(Year,Coefficient)
        ) + geom_boxplot(varwidth=TRUE,
                         outlier.shape=NA,
                         lwd=0.3
        ) + labs(x='Study Year',y="Coefficient"
        ) + coord_cartesian(ylim = c(min(quants$q25)-(1.5*max(iqr)),(1.5*max(iqr))+max(quants$q75))
        ) + theme(axis.line = element_line(size=0.2),
                  panel.background = element_blank())
        
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_YearPlot.png",sep="")
        ggsave(saveName,device="png", width=2, scale=3.5, height=1, units="in",dpi=600)
        while (dev.cur()>1) {dev.off()}
        saveName = paste(outDir,'/',CTname,'/',sites[j],"_",int,"_GEEGLM_YearPlot.pdf",sep="")
        ggsave(saveName,device="pdf", width=2, scale=3.5, height=1, units="in",dpi=600)
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
      
      # Calculate proportion of presence bins during day, during period of peak presence
      if (PresStats[[CTname]]$JDPeak[siteInd]!=0){ # If we know when seasonal peak presence is
        # grab 90 days centered around peak presence
        peakPres = seq(round(PresStats[[CTname]]$JDPeak[siteInd])-45,round(PresStats[[CTname]]$JDPeak[siteInd]+45))
        peakPres[peakPres<0] = peakPres[peakPres<0]+365
        peakPres[peakPres>365] = peakPres[peakPres>365]-365
        peakPresBins = which(thisSite$JulianDay%in%peakPres)
      } else { 
        # Otherwise look across the whole year
        peakPresBins = seq(1,length(thisSite$JulianDay))}
      
      # during time of interest, find presence bins
      presBins = which(thisSite$Presence[peakPresBins]==1)
      # identify which presence bins are during the day
      dayBins = which(thisSite$NormTime[peakPresBins[presBins]]<0)
      # calculate proportion of presence bins are during the day
      PresStats[[CTname]]$DielProp[siteInd] = length(dayBins)/length(presBins)
    }
  }
  
  gotMods = which(rownames(PresStats[[CTname]])%in%sites)
  noMods = setdiff(seq(1,11),gotMods)
  PresStats[[CTname]][noMods,] = rep(NA,ncol(PresStats[[CTname]]))
  
  
}

colnames(modPerf) = c("Species","Site","PropGoodResid")
save(modPerf,file=paste(outDir,'/TempLun3ModPerf.Rdata',sep=""))
save(PresStats,file=paste(outDir,'/PresStats.Rdata',sep=""))


## Old model fit plots ---------------------
# modFit = data.frame(Date=as.Date(thisSite$TimeStamp),
#                     Pres=thisSite$Presence,
#                     Fits=tempMod$fitted.values,
#                     Res=tempMod$residuals)
# 
# PresPlot = ggplot(modFit
# ) + geom_point(aes(x=Date,y=Pres),
#                color="#000000",
#                size=2
# )+scale_x_continuous(breaks=c(as.Date("2016-05-01"),
#                               as.Date("2016-11-01"),
#                               as.Date("2017-05-01"),
#                               as.Date("2017-11-01"),
#                               as.Date("2018-05-01"),
#                               as.Date("2018-11-01"),
#                               as.Date("2019-04-30"))
# ) + labs(x="",y="Presence")
# FitPlot = ggplot(modFit
# )+geom_point(aes(x=Date,y=Fits),
#              color="#1976D2",
#              size=2
# )+scale_x_continuous(breaks=c(as.Date("2016-05-01"),
#                               as.Date("2016-11-01"),
#                               as.Date("2017-05-01"),
#                               as.Date("2017-11-01"),
#                               as.Date("2018-05-01"),
#                               as.Date("2018-11-01"),
#                               as.Date("2019-04-30"))
# ) + labs(x="",y="Fitted Values")
# ModRes = ggplot(modFit
# ) + geom_point(aes(x=Date,y=Res),
#                size=2
# ) +scale_x_continuous(breaks=c(as.Date("2016-05-01"),
#                                as.Date("2016-11-01"),
#                                as.Date("2017-05-01"),
#                                as.Date("2017-11-01"),
#                                as.Date("2018-05-01"),
#                                as.Date("2018-11-01"),
#                                as.Date("2019-04-30"))
# ) + labs(x="",y="Residuals")
# 
# HL = hoslem.test(thisSite$Presence,fitted(tempMod))
# HLPlot = data.frame(Expected=HL$expected[,1],Observed=HL$observed[,1])
# HLPlot$X = seq(min(HL$observed[,1]),max(HL$observed[,1]),length.out=length(HLPlot$Expected))
# HLPlot$Y = seq(min(HL$observed[,1]),max(HL$observed[,1]),length.out=length(HLPlot$Expected))
#
# ObsPred = ggplot(HLPlot
# ) + geom_line(aes(x=X,y=Y)
# ) + geom_point(aes(y=Expected,x=Observed)
# ) + labs(x="Observed",y="Expected"
# ) + annotate("text",
#              label=paste('p-value: ',as.character(HL$p.value),sep=""),
#              size=4,
#              x=1.005*min(HL$observed[,1]),
#              y=0.999*max(HL$expected[,1]))
#
# png(file=paste(outDir,'/',CTname,'/',sites[j],"_ModelFit.png",sep=""),width = 1400, height = 800, units = "px")
# grid.arrange(PresPlot,FitPlot,ModRes,ncol=1,nrow=3,top=paste(CTname,'at',sites[j]))
# while (dev.cur()>1) {dev.off()}
