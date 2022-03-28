library(tidyverse)
library(lubridate)
library(geepack)
library(splines2)
library(pracma)
library(mgcv)
library(car)
library(stringr)
# library(corrplot)
# library(scales)
# library(patchwork)
source("getPvalues.r")


modelDFDir = 'J:/Chpt_2/TimeSeries_ScaledByEffortError'
outDir = 'J:/Chpt_2/ModelOutput'
int = "5minBin"

dfList = list.files(path=modelDFDir,pattern=paste('*',int,'_Master.csv',sep=""),
                    full.names=TRUE,recursive=FALSE,
                    include.dirs=FALSE,no..=TRUE)
lunList = list.files(path=modelDFDir,pattern=paste('*',int,'_MasterLun.csv',sep=""),
                     full.names=TRUE,recursive=FALSE,
                     include.dirs=FALSE,no..=TRUE)

#### EXPLORE SPLINE FITS ---------------------------------

## Julian Day + Normalized Time of Day --------------

# QIC_Comp = data.frame(NT4=as.numeric(),
#                       NT5=as.numeric(),
#                       JD4=as.numeric(),
#                       JD5=as.numeric())
# 
# for (i in 1:numel(dfList)){
# 
#   # load file
#   thisSite = data.frame(read.csv(dfList[i]))
#   # CTname = str_remove(dfList[i],paste(modelDFDir,'/',sep="")) # get the species/CT name
#   # site = str_remove(CTname,paste("_",int,"_Master.csv",sep=""))
#   # site = sub(".*_","",site)
#   # CTname = sub("_.*","",CTname)
#   # if (str_detect(CTname,"Atl")){
#   #   CTname = "Gervais"
#   # }
# 
#   # if it doesn't already exist, create directory to save figures
#   if (!dir.exists(paste(outDir,'/',CTname,sep=""))){
#     dir.create(paste(outDir,'/',CTname,sep=""))
#   }
# 
#   knots = list(c(0.333,0.666),c(0.275,0.5,0.725))
# 
#   # Test how many knots to use for NormToD (not doing linear cause it should be cyclic)
#   mod01 = geeglm(Presence~mSpline(NormTime,
#                                    knots=quantile(thisSite$NormTime,probs=unlist(knots[1])),
#                                    Boundary.knots=c(min(thisSite$NormTime),max(thisSite$NormTime)),
#                                    periodic=T),
#                  family=binomial,data=thisSite,id=GroupID,corstr="ar1")
#   mod02 = geeglm(Presence~mSpline(NormTime,
#                                    knots=quantile(thisSite$NormTime,probs=unlist(knots[2])),
#                                    Boundary.knots=c(min(thisSite$NormTime),max(thisSite$NormTime)),
#                                    periodic=T),
#                  family=binomial,data=thisSite,id=GroupID,corstr="ar1")
#   QICNormTime = c(QIC(mod01)[[1]],QIC(mod02)[[1]])
# 
#   # Test how many knots to use for Julian Day (not doing linear cause it should be cyclic)
#   mod03 = geeglm(Presence~mSpline(JulianDay,
#                                   knots=quantile(thisSite$JulianDay, probs=unlist(knots[1])),
#                                   Boundary.knots=c(1,365),
#                                   periodic=T),
#                  family=binomial,data=thisSite,id=GroupID,corstr="ar1")
#   mod04 = geeglm(Presence~mSpline(JulianDay,
#                                   knots=quantile(thisSite$JulianDay, probs=unlist(knots[2])),
#                                   Boundary.knots=c(1,365),
#                                   periodic=T),
#                  family=binomial,data=thisSite,id=GroupID,corstr="ar1")
#   QICJD = c(QIC(mod03)[[1]],QIC(mod04)[[1]])
# 
# 
#   QIC_Comp[i,] = cbind(QIC(mod01)[[1]],QIC(mod02)[[1]],QIC(mod03)[[1]],
#                        QIC(mod04)[[1]])
# 
# 
# }
# 
# BestNormTime = data.frame(colnames(QIC_Comp)[unlist(data.frame(apply(QIC_Comp[,1:2],1,which.min)))])
# BestJD = data.frame(colnames(QIC_Comp)[unlist(data.frame(apply(QIC_Comp[,3:4],1,which.min)))+2])
# colnames(BestNormTime) = "BestModel"
# colnames(BestJD) = "BestModel"
# save(dfList,QIC_Comp,BestNormTime,BestJD,file=paste(outDir,"/JDNT_SmoothFitEval.Rdata",sep=""))

## Lunar ----------------

# Lun_QIC_Comp = data.frame(Phase4=as.numeric(),
#                                Phase5=as.numeric(),
#                                AltLin=as.numeric(),
#                                Alt4=as.numeric(),
#                                Alt5=as.numeric())

# for (i in 1:numel(lunList)){
for (i in 65:numel(lunList)){
  

  # load file
  masterLun = data.frame(read.csv(lunList[i]))

  masterLun$MoonAltitude = masterLun$MoonAltitude*(180/pi) # convert altitude from radians to degrees

  knots = list(c(0.333,0.666),c(0.275,0.5,0.725))
  # Test whether to include Moon Phase with 4 or 5 knots (not doing linear cause it should be cyclic)
  mod01 = geeglm(NightPres~mSpline(MoonPhase,
                                   knots=quantile(masterLun$MoonPhase,probs=unlist(knots[1])),
                                   Boundary.knots=c(min(masterLun$MoonPhase),max(masterLun$MoonPhase)),
                                   periodic=T),
                 family=binomial,data=masterLun,id=GroupID,corstr="ar1")
  mod02 = geeglm(NightPres~mSpline(MoonPhase,
                                   knots=quantile(masterLun$MoonPhase,probs=unlist(knots[2])),
                                   Boundary.knots=c(min(masterLun$MoonPhase),max(masterLun$MoonPhase)),
                                   periodic=T),
                 family=binomial,data=masterLun,id=GroupID,corstr="ar1")
  QICMoonPhase = c(QIC(mod01)[[1]],QIC(mod02)[[1]])


  # Test whether to include Moon Altitude as a linear or a smooth term with 4 or 5 knots:
  mod03 = geeglm(NightPres~MoonAltitude,family=binomial,data=masterLun,id=GroupID,corstr="ar1")
  mod04 = geeglm(NightPres~mSpline(MoonAltitude,
                                   knots=quantile(masterLun$MoonAltitude,probs=unlist(knots[1])),
                                   Boundary.knots=c(min(masterLun$MoonAltitude),max(masterLun$MoonAltitude))),
                 family=binomial,data=masterLun,id=GroupID,corstr="ar1")
  mod05 = geeglm(NightPres~mSpline(MoonAltitude,
                                   knots=quantile(masterLun$MoonAltitude,probs=unlist(knots[2])),
                                   Boundary.knots=c(min(masterLun$MoonAltitude),max(masterLun$MoonAltitude))),
                 family=binomial,data=masterLun,id=GroupID,corstr="ar1")
  QICMoonAltitude = c(QIC(mod03)[[1]],QIC(mod04)[[1]],QIC(mod05)[[1]])


  Lun_QIC_Comp[i,] = cbind(QIC(mod01)[[1]],QIC(mod02)[[1]],QIC(mod03)[[1]],
                       QIC(mod04)[[1]],QIC(mod05)[[1]])

}

BestPhase = data.frame(colnames(Lun_QIC_Comp)[unlist(data.frame(apply(Lun_QIC_Comp[,1:2],1,which.min)))])
colnames(BestPhase) = "BestModel"
BestAlt = data.frame(colnames(Lun_QIC_Comp)[(unlist(data.frame(apply(Lun_QIC_Comp[,3:5],1,which.min))))+2])
colnames(BestAlt) = "BestModel"
save(lunList,Lun_QIC_Comp,BestPhase,BestAlt,file=paste(outDir,"/MPMA_SmoothFitEval.Rdata",sep=""))

