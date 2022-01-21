library(tidyverse)
library(lubridate)
library(geepack)
library(splines2)
library(SimDesign)
library(pracma)
library(mgcv)
library(car)
library(corrplot)
library(scales)
library(patchwork)
library(boot)
source("getPvalues.r")


modelDFDir = 'I:/TimeSeries_ScaledByEffortError'
outDir = 'I:/ModelOutput'
int = "5minBin"

dfList = list.files(path=modelDFDir,pattern=paste('*Master.csv',sep=""),
                    full.names=TRUE,recursive=FALSE,
                    include.dirs=FALSE,no..=TRUE)


#### EXPLORE SPLINE FITS ---------------------------------

# numKnots_QIC_Comp = data.frame(LunLin=as.numeric(),
#                                Lun4=as.numeric(),
#                                Lun5=as.numeric(),
#                                JD4=as.numeric(),
#                                JD5=as.numeric())
# for (i in 1:numel(dfList)){ 
#   
#   # load file
#   thisSite = data.frame(read.csv(dfList[i]))
#   CTname = str_remove(dfList[i],paste(modelDFDir,'/',sep="")) # get the species/CT name
#   site = str_remove(CTname,paste("_",int,"_Master.csv",sep=""))
#   site = sub(".*_","",site)
#   CTname = sub("_.*","",CTname)
#   if (str_detect(CTname,"Atl")){
#     CTname = "Gervais"
#   }
#   
#   # if it doesn't already exist, create directory to save figures
#   if (!dir.exists(paste(outDir,'/',CTname,sep=""))){
#     dir.create(paste(outDir,'/',CTname,sep=""))
#   }
#   
#   # # Fit GAM
#   # if (j==7){
#   #   JDayGAM = gam(thisSite~s(Jday,bs="cc",k=6)+as.factor(yearGroup)+as.factor(hatSite),
#   #                 family=tw)
#   # } else {JDayGAM = gam(thisSite~s(Jday,bs="cc",k=6)+as.factor(yearGroup),
#   #               family=tw)}
#   # sinkName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GAMSummary.txt",sep="")
#   # sink(sinkName)
#   # print(summary(JDayGAM))
#   # sink()
#   # 
#   # # Plot GAM partial residuals
#   # saveName = paste(seasDir,'/',CTname,'/',sites[j],"_",int,"_GAM.png",sep="")
#   # if (j==7){
#   #   png(saveName,width=800,height=500)
#   # } else {
#   #   png(saveName,width=800,height=400)
#   # }
#   # 
#   # plot.gam(JDayGAM,all.terms=TRUE,pages=1,main=paste(CTname,'at',sites[j]))
#   # while (dev.cur()>1) {dev.off()}
#   
#   
#   knots = list(c(0.333,0.666),c(0.275,0.5,0.725))
#   # Test whether to include Lunar Illuminance as a linear or a smooth term:
#   mod01 = geeglm(Presence~LunarIllum,family=binomial,data=thisSite,id=GroupID,corstr="ar1")
#   mod02 = geeglm(Presence~mSpline(LunarIllum,
#                                   knots=quantile(thisSite$LunarIllum,probs=unlist(knots[1])),
#                                   Boundary.knots=c(0,1)),
#                  family=binomial,data=thisSite,id=GroupID,corstr="ar1")
#   mod03 = geeglm(Presence~mSpline(LunarIllum,
#                                   knots=quantile(thisSite$LunarIllum,probs=unlist(knots[2])),
#                                   Boundary.knots=c(0,1)),
#                  family=binomial,data=thisSite,id=GroupID,corstr="ar1")
#   QICLun = c(QIC(mod01)[[1]],QIC(mod02)[[1]],QIC(mod03)[[1]])
#   # QICLun = c(QIC(mod01)[[1]],QIC(mod02)[[1]])
#   
#   # Test how many knots to use for Julian Day (not doing linear cause it should be cyclic)
#   mod04 = geeglm(Presence~mSpline(JulianDay,
#                                   knots=quantile(thisSite$JulianDay, probs=unlist(knots[1])),
#                                   Boundary.knots=c(1,365),
#                                   periodic=T),
#                  family=binomial,data=thisSite,id=GroupID,corstr="ar1")
#   mod05 = geeglm(Presence~mSpline(JulianDay,
#                                   knots=quantile(thisSite$JulianDay, probs=unlist(knots[2])),
#                                   Boundary.knots=c(1,365),
#                                   periodic=T),
#                  family=binomial,data=thisSite,id=GroupID,corstr="ar1")
#   QICJD = c(QIC(mod04)[[1]],QIC(mod05)[[1]])
#   
#   
#   numKnots_QIC_Comp[i,] = cbind(QIC(mod01)[[1]],QIC(mod02)[[1]],QIC(mod03)[[1]],
#                        QIC(mod04)[[1]],QIC(mod05)[[1]])
#   
# }
# 
# BestLun = data.frame(colnames(numKnots)[unlist(data.frame(apply(numKnots[,1:3],1,which.min)))])
# colnames(BestLun) = "BestModel"
# BestJD = data.frame(colnames(numKnots)[unlist(data.frame(apply(numKnots[,4:5],1,which.min)))+3])
# colnames(BestJD) = "BestModel"
# save(dfList,numKnots_QIC_Comp,BestLun,BestJD,file=paste(outDir,"/SmoothFitEval.Rdata",sep=""))


# ##### CREATE MODELS (JD, LunIllum, Year) --------------------------
## Best # knots for each continuous covar, selected based on which most commonly had lower QIC
# smoothKnots = data.frame(Species = c("Blainville","Cuvier","Gervais","Kogia",
#                                      "Risso","Sowerby","Sperm Whale","True",
#                                      "UD19","UD26","UD28","UD36"),
#                          Lunar = c(4,4,4,4,5,4,5,5,5,4,4,5),
#                          JD = c(5,5,5,5,5,5,5,5,5,5,5,5))
# 
# for (i in 1:numel(dfList)){
# 
#     # load file
#     thisSite = data.frame(read.csv(dfList[i]))
#     CTname = str_remove(dfList[i],paste(modelDFDir,'/',sep="")) # get the species/CT name
#     site = str_remove(CTname,paste("_",int,"_Master.csv",sep=""))
#     site = sub(".*_","",site)
#     CTname = sub("_.*","",CTname)
#     if (str_detect(CTname,"Atl")){
#       CTname = "Gervais"
#     }
# 
# # Fit GEEGLM:
#     knots = list(c(NaN),c(0.333,0.666),c(0.275,0.5,0.725))
#     knotInd = which(smoothKnots$Species==CTname)
#     LunKnots = unlist(knots[smoothKnots$Lunar[knotInd]-2])
#     JDKnots = unlist(knots[smoothKnots$JD[knotInd]-2])
# 
#   if (is.na(LunKnots)){ # include Lunar Illuminance as a linear
#   tempMod = geeglm(Presence~mSpline(JulianDay,
#                                     knots=quantile(JulianDay,probs=JDKnots),
#                                     Boundary.knots=c(1,365),
#                                     periodic=T)
#                    +LunarIllum
#                    +as.factor(StudyYear),
#                    #+as.factor(DayPhase),
#                    family=binomial,
#                    data=thisSite,
#                    id=GroupID,
#                    corstr="ar1")
#   LunInd = 4
#   } else {
#     tempMod = geeglm(Presence~mSpline(JulianDay,
#                                       knots=quantile(JulianDay,probs=JDKnots),
#                                       Boundary.knots=c(1,365),
#                                       periodic=T)
#                      +mSpline(LunarIllum,
#                               knots=quantile(LunarIllum,probs=LunKnots),
#                               Boundary.knots=c(0,1))
#                      +as.factor(StudyYear),
#                      #+as.factor(DayPhase),
#                      family=binomial,
#                      data=thisSite,
#                      id=GroupID,
#                      corstr="ar1")
#   }
# 
#   PV = getPvalues(tempMod)
#   sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_4knots.txt",sep="")
#   sink(sinkName)
#   print(PV)
#   sink()
# 
#   save(tempMod,PV,JDKnots,LunKnots,file=paste(outDir,'/',CTname,'/',site,"_",int,"_Model_4knots.Rdata",sep=""))
# }

##### CREATE MODELS (JD & Year) --------------------------

for (i in 1:numel(dfList)){
  
  # load file
  thisSite = data.frame(read.csv(dfList[i]))
  CTname = str_remove(dfList[i],paste(modelDFDir,'/',sep="")) # get the species/CT name
  site = str_remove(CTname,paste("_",int,"_Master.csv",sep=""))
  site = sub(".*_","",site)
  CTname = sub("_.*","",CTname)
  if (str_detect(CTname,"Atl")){
    CTname = "Gervais"
  }
  
  # Find nighttime data
  nightInd = which(thisSite$DayPhase=='Night')
  gaps = which(diff(nightInd,lag=1)>1)
  nightEnd = nightInd[gaps]
  nightSt = nightInd[gaps+1]
  nightEnd = nightEnd[2:length(nightEnd)]
  nightEnd = c(nightEnd,nightInd[length(nightInd)])
  
  # calculate proportion of each night w clicks
  
  # calculate proportion of each night w moon
  
  # calculate average lunar illuminance each night
  
  # make data frame for modeling

  
  # Fit GEEGLM:
  
  tempMod = geeglm(Presence~mSpline(JulianDay,
                                    knots=quantile(JulianDay,probs=c(0.275,0.5,0.725)),
                                    Boundary.knots=c(1,365),
                                    periodic=T)
                   +as.factor(StudyYear),
                   family=binomial,
                   data=thisSite,
                   id=GroupID,
                   corstr="ar1")
  
  PV = getPvalues(tempMod)
  sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_JDYr.txt",sep="")
  sink(sinkName)
  print(PV)
  sink()
  
  save(tempMod,PV,file=paste(outDir,'/',CTname,'/',site,"_",int,"_Model_JDYr.Rdata",sep=""))
}

##### CREATE MODELS (LunIllum) ------------------------------------

for (i in 1:numel(dfList)){
  
  # load file
  thisSite = data.frame(read.csv(dfList[i]))
  CTname = str_remove(dfList[i],paste(modelDFDir,'/',sep="")) # get the species/CT name
  site = str_remove(CTname,paste("_",int,"_Master.csv",sep=""))
  site = sub(".*_","",site)
  CTname = sub("_.*","",CTname)
  if (str_detect(CTname,"Atl")){
    CTname = "Gervais"
  }
  
  
  
  # Fit GEEGLM:
  
  tempMod = geeglm(Presence~mSpline(JulianDay,
                                    knots=quantile(JulianDay,probs=c(0.275,0.5,0.725)),
                                    Boundary.knots=c(1,365),
                                    periodic=T)
                   +as.factor(StudyYear),
                   family=binomial,
                   data=thisSite,
                   id=GroupID,
                   corstr="ar1")
  
  PV = getPvalues(tempMod)
  sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_JDYr.txt",sep="")
  sink(sinkName)
  print(PV)
  sink()
  
  save(tempMod,PV,file=paste(outDir,'/',CTname,'/',site,"_",int,"_Model_JDYr.Rdata",sep=""))
}


# #### PLOT ALL MODELS FOR A GIVEN SPECIES --------------------------------------

# species = list.dirs(outDir,recursive=FALSE)
# 
# for ( i in 1:numel(species)){
# 
# 
#   modFiles = list.files(path=species[i],pattern="*Model.Rdata",
#                         full.names=TRUE,recursive=FALSE,include.dirs=FALSE,no..=TRUE)
#   CTname = str_remove(species[i],paste(outDir,'/',sep=""))
#   sites = list()
# 
#   # JDList = list()
#   # LunIllumList = list()
#   # YrList = list()
# 
#   for (j in 1:numel(modFiles)){
# 
#     site = str_remove(modFiles[j],paste(outDir,"/",CTname,"/",sep=""))
#     sites = c(sites,str_remove(site,"_5minBin_Model.Rdata"))
# 
#     load(modFiles[j]) # load model
#     # find associated master dataframe
#     thisSpec = which(str_detect(dfList,CTname))
#     thisSite = which(str_detect(dfList,unlist(sites[j])))
#     thisModInd = intersect(thisSpec,thisSite)
#     thisSite = data.frame(read.csv(dfList[thisModInd]))
# 
#     # Get indices of coefficients for each covar
#     if (numel(JDKnots)==2){
#       JDind = c(2,3)
#     } else {JDind = c(2:4)}
#     if (numel(LunKnots)==2){
#       LunInd = JDind[length(JDind)]+c(1:5)
#     } else if (numel(LunKnots)==3){
#       LunInd = JDind[length(JDind)]+c(1:6)
#     }
#     YrInd = c(1,LunInd[length(LunInd)]+c(1:2))
#     #PhsInd = c(1,YrInd[length(YrInd)]+c(1:3))
# 
#     ## Bootstrap GEEGLM parameter estimates for later construction of confidence intervals ----------------
#     BootstrapParameters<-rmvnorm(10000, coef(tempMod), summary(tempMod)$cov.unscaled)
#     JDayBootstrapCoefs<- BootstrapParameters[,JDind]
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
#     JxJD = model.matrix(tempMod)[,JDind]%*%coef(tempMod)[JDind]
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
#     RealFitJ<- JBasis%*%coef(tempMod)[JDind] # multiply basis functions by model coefficients to get values of spline at each X
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
