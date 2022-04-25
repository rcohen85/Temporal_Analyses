library(tidyverse)
library(lubridate)
library(geepack)
library(splines2)
library(pracma)
library(mgcv)
library(car)
library(stringr)
source("getPvalues.r")

start = Sys.time()
modelDFDir = 'J:/Chpt_2/TimeSeries_ScaledByEffortError'
outDir = 'J:/Chpt_2/ModelOutput'
int = "5minBin"

dfList = list.files(path=modelDFDir,pattern=paste('*',int,'_MasterTempLun.csv',sep=""),
                    full.names=TRUE,recursive=FALSE,
                    include.dirs=FALSE,no..=TRUE)
# lunList = list.files(path=modelDFDir,pattern=paste('*',int,'_MasterLun.csv',sep=""),
#                      full.names=TRUE,recursive=FALSE,
#                      include.dirs=FALSE,no..=TRUE)

##### CREATE MODELS USING FIXED SPLINES --------------------------

### JD, MPh, NT, MPr, Yr models ----------------------------------
# for (i in 1:length(dfList)){
for (i in 76:length(dfList)){
  
# load file
thisSite = data.frame(read.csv(dfList[i]))
CTname = str_remove(dfList[i],paste(modelDFDir,'/',sep="")) # get the species/CT name
site = str_remove(CTname,paste("_",int,"_MasterTempLun.csv",sep=""))
site = sub(".*_","",site)
CTname = sub("_.*","",CTname)
if (str_detect(CTname,"Atl")){
  CTname = "Gervais"
}

# Prepare covars:

if (site=="HAT"){ # for HAT only model 2017-2019
  keepInd = which(thisSite$StudyYear>1)
  Pres = thisSite$Presence[keepInd]
  GroupID = thisSite$GroupID[keepInd]
  NTs = mSpline(thisSite$NormTime[keepInd],
                knots=quantile(thisSite$NormTime[keepInd],probs=c(0.275,0.5,0.725)),
                Boundary.knots=c(-1,1),
                periodic=T)
  MPhs = mSpline(thisSite$MoonPhase[keepInd],
                 knots=quantile(thisSite$MoonPhase[keepInd],probs=c(0.275,0.5,0.725)),
                 Boundary.knots=c(0,1),
                 periodic=T)
  JDs = mSpline(thisSite$JulianDay[keepInd],
                knots=quantile(thisSite$JulianDay[keepInd],probs=c(0.275,0.5,0.725)),
                Boundary.knots=c(1,366),
                periodic=T)
  MPr = as.factor(thisSite$MoonPres[keepInd])
  # LF = as.factor(thisSite$LunFact[keepInd])
  YrF = as.factor(thisSite$StudyYear[keepInd])

} else { # otherwise use all data
  Pres = thisSite$Presence
  GroupID = thisSite$GroupID
  NTs = mSpline(thisSite$NormTime,
                knots=quantile(thisSite$NormTime,probs=c(0.275,0.5,0.725)),
                Boundary.knots=c(-1,1),
                periodic=T)
  MPhs = mSpline(thisSite$MoonPhase,
                 knots=quantile(thisSite$MoonPhase,probs=c(0.275,0.5,0.725)),
                 Boundary.knots=c(0,1),
                 periodic=T)
  JDs = mSpline(thisSite$JulianDay,
                knots=quantile(thisSite$JulianDay,probs=c(0.275,0.5,0.725)),
                Boundary.knots=c(1,366),
                periodic=T)
  MPr = as.factor(thisSite$MoonPres)
  # LF = as.factor(thisSite$LunFact)
  YrF = as.factor(thisSite$StudyYear)
}

# Fit full model
tempMod = geeglm(Pres~NTs
                 +JDs
                 +MPhs
                 # +LF
                 +YrF
                 +NTs:JDs
                 +NTs:MPhs,
                 # +MPhs:LF,
                 family=binomial,
                 id=GroupID,
                 corstr="ar1")


PV = getPvalues(tempMod)

# Iteratively remove non-signif covars, re-run model to arrive at final model
if (sum(PV$'p-value'<0.05)==length(PV$'p-value')){ # if all variables are significant, save
  if (tempMod$geese$error==0){ # if model converged, save model and summary
    sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_TempLun3.txt",sep="")
    sink(sinkName)
    print(PV)
    sink()
    save(tempMod,PV,file=paste(outDir,'/',CTname,'/',site,"_",int,"_Model_TempLun3.Rdata",sep=""))
  } else { # if model did not converge, don't save model, only summary indicating non-convergence
    sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_TempLun3.txt",sep="")
    sink(sinkName)
    print('Model did not converge')
    sink()
  }
} else if (sum(PV$'p-value'<0.05)==0){ # if no significant variables, only save summary indicating non-significance
  sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_TempLun3.txt",sep="")
  sink(sinkName)
  print('No significance')
  sink()
} else {
  PV2 = PV
  while (sum(PV2$'p-value'<0.05) < length(PV2$'p-value')){
    # if interaction term is significant, keep both contributing terms in model
    if (any(PV2$'Variable'=="NTs:JDs") & length(PV2$'p-value'[PV2$'Variable'=="NTs:JDs"]<0.05)>0){ 
      PV2$'p-value'[1:2] = 0.01
    }
    if (any(PV2$'Variable'=="NTs:MPhs") & length(PV2$'p-value'[PV2$'Variable'=="NTs:MPhs"]<0.05)>0){
      PV2$'p-value'[c(1,3)] = 0.01
    }
    # if (any(PV2$'Variable'=="MPhs:LF") & length(PV2$'p-value'[PV2$'Variable'=="MPhs:LF"]<0.05)>0){
    #   PV2$'p-value'[3:4] = 0.01
    # }
    whichVars = which(PV2$'p-value'<0.05)
    str1 = "tempMod = geeglm(Pres~"
    str2 = ""
    str3 = ",family=binomial,id=GroupID,corstr='ar1')"
    for (j in 1:length(whichVars)){
      if(j==1){
        str2 = paste(str2,PV2$Variable[whichVars[j]],sep="")
      } else {
        str2 = paste(str2,PV2$Variable[whichVars[j]],sep="+")
      }
    }
    modelCall = paste(str1,str2,str3,sep="")
    eval(parse(text=modelCall))
    PV = getPvalues(tempMod)
    PV2 = PV
    # if interaction term is significant, keep both contributing terms in model
    if (any(PV2$'Variable'=="NTs:JDs") & length(PV2$'p-value'[PV2$'Variable'=="NTs:JDs"]<0.05)>0){ 
      PV2$'p-value'[1:2] = 0.01
    }
    if (any(PV2$'Variable'=="NTs:MPhs") & length(PV2$'p-value'[PV2$'Variable'=="NTs:MPhs"]<0.05)>0){
      PV2$'p-value'[c(1,3)] = 0.01
    }
    # if (any(PV2$'Variable'=="MPhs:LF") & length(PV2$'p-value'[PV2$'Variable'=="MPhs:LF"]<0.05)>0){
    #   PV2$'p-value'[3:4] = 0.01
    # }
  }
  if (tempMod$geese$error==0){ # if model converged, save model and summary
    sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_TempLun3.txt",sep="")
    sink(sinkName)
    print(PV)
    sink()
    save(tempMod,PV,file=paste(outDir,'/',CTname,'/',site,"_",int,"_Model_TempLun3.Rdata",sep=""))
  } else { # if model did not converge, don't save model, only summary indicating non-convergence
    sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_TempLun3.txt",sep="")
    sink(sinkName)
    print('Model did not converge')
    sink()
  }
}

end = Sys.time()
end-start
}
### JD, Yr, NormTime Models --------------------------
# for (i in 1:numel(dfList)){
#   # 
#   # # load file
#   # thisSite = data.frame(read.csv(dfList[i]))
#   # CTname = str_remove(dfList[i],paste(modelDFDir,'/',sep="")) # get the species/CT name
#   # site = str_remove(CTname,paste("_",int,"_Master.csv",sep=""))
#   # site = sub(".*_","",site)
#   # CTname = sub("_.*","",CTname)
#   # if (str_detect(CTname,"Atl")){
#   #   CTname = "Gervais"
#   # }
#   # 
#   # # Fit JD + Year GEEGLM:
#   # 
#   # if (site=="HAT"){ # for HAT only model 2017-2019
#   #   startInd = which(thisSite$StudyYear>1)
#   #   JDs = mSpline(thisSite$JulianDay[startInd],
#   #                 knots=quantile(thisSite$JulianDay[startInd],probs=c(0.275,0.5,0.725)),
#   #                 Boundary.knots=c(1,365),
#   #                 periodic=T)
#   #   NTs = mSpline(thisSite$NormTime[startInd],
#   #                 knots=quantile(thisSite$NormTime[startInd],probs=c(0.275,0.5,0.725)),
#   #                 Boundary.knots=c(-1,1),
#   #                 periodic=T)
#   #   YrF = as.factor(thisSite$StudyYear[startInd])
#   #   Pres = thisSite$Presence[startInd]
#   #   GroupID = thisSite$GroupID[startInd]
#   #   
#   # } else { # otherwise use all data
#   #   Pres = thisSite$Presence
#   #   GroupID = thisSite$GroupID
#   #   JDs = mSpline(thisSite$JulianDay,
#   #                 knots=quantile(thisSite$JulianDay,probs=c(0.275,0.5,0.725)),
#   #                 Boundary.knots=c(1,365),
#   #                 periodic=T)
#   #   NTs = mSpline(thisSite$NormTime,
#   #                 knots=quantile(thisSite$NormTime,probs=c(0.275,0.5,0.725)),
#   #                 Boundary.knots=c(-1,1),
#   #                 periodic=T)
#   #   YrF = as.factor(thisSite$StudyYear)
#   # }
#   # 
#   # tempMod = geeglm(Pres~JDs
#   #                  +NTs
#   #                  +YrF
#   #                  +JDs:NTs,
#   #                  family=binomial,
#   #                  id=GroupID,
#   #                  corstr="ar1")
#   # 
#   # 
#   # PV = getPvalues(tempMod)
#   # 
#   # if (sum(PV$'p-value'<0.05)==length(PV$'p-value')){ # if all variables are significant, save
#   #   if (tempMod$geese$error==0){ # if model converged, save model and summary
#   #     sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_JDNTYr_5I.txt",sep="")
#   #     sink(sinkName)
#   #     print(PV)
#   #     sink()
#   #     save(tempMod,PV,file=paste(outDir,'/',CTname,'/',site,"_",int,"_Model_JDNTYr_5I.Rdata",sep=""))
#   #   } else { # if model did not converge, don't save model, only summary indicating non-convergence
#   #     sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_JDNTYr_5I.txt",sep="")
#   #     sink(sinkName)
#   #     print('Model did not converge')
#   #     sink()
#   #   }
#   #   
#   # } else {
#   #   PV2 = PV
#   #   while (sum(PV2$'p-value'<0.05) < length(PV2$'p-value')){
#   #     if (PV2$'p-value'[PV2$'Variable'=="JDs:NTs"]<0.05){ # if interaction term is significant, keep both contributing terms in model
#   #       PV2$'p-value'[1:2] = 0.01
#   #     }
#   #     whichVars = which(PV2$'p-value'<0.05)
#   #     str1 = "tempMod = geeglm(Pres~"
#   #     str2 = ""
#   #     str3 = ",family=binomial,id=GroupID,corstr='ar1')"
#   #     for (j in 1:length(whichVars)){
#   #       if(j==1){
#   #         str2 = paste(str2,PV2$Variable[whichVars[j]],sep="")
#   #       } else {
#   #         str2 = paste(str2,PV2$Variable[whichVars[j]],sep="+")
#   #       }
#   #     }
#   #     modelCall = paste(str1,str2,str3,sep="")
#   #     eval(parse(text=modelCall))
#   #     PV = getPvalues(tempMod)
#   #     PV2 = PV
#   #     if (any(PV2$'Variable'=="JDs:NTs")){
#   #       if (PV2$'p-value'[PV2$'Variable'=="JDs:NTs"]<0.05){ # if interaction term is significant, keep both contributing terms in model
#   #         PV2$'p-value'[1:2] = 0.01
#   #       }
#   #     }
#   #   }
#   #   if (tempMod$geese$error==0){ # if model converged, save model and summary
#   #     sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_JDNTYr_5I.txt",sep="")
#   #     sink(sinkName)
#   #     print(PV)
#   #     sink()
#   #     save(tempMod,PV,file=paste(outDir,'/',CTname,'/',site,"_",int,"_Model_JDNTYr_5I.Rdata",sep=""))
#   #   } else { # if model did not converge, don't save model, only summary indicating non-convergence
#   #     sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_JDNTYr_5I.txt",sep="")
#   #     sink(sinkName)
#   #     print('Model did not converge')
#   #     sink()
#   #   }
#   # }
# }

### Lunar Models ------------------
# for (i in 1:numel(lunList)){
#   # 
#   # masterLun = data.frame(read.csv(lunList[i]))
#   # CTname = str_remove(lunList[i],paste(modelDFDir,'/',sep="")) # get the species/CT name
#   # site = str_remove(CTname,paste("_",int,"_MasterLun.csv",sep=""))
#   # site = sub(".*_","",site)
#   # CTname = sub("_.*","",CTname)
#   # if (str_detect(CTname,"Atl")){
#   #   CTname = "Gervais"
#   # }
#   # 
#   # masterLun$MoonAltitude = masterLun$MoonAltitude*(180/pi) # convert altitude from radians to degrees
#   # 
#   # # Fit GEEGLM:
#   # if (site=="HAT"){
#   #   
#   #   startInd = which(masterLun$NightBinTimes>=as.POSIXct('2017-05-01 00:00:00',format="%Y-%m-%d %H:%M:%S",tz="GMT"))
#   #   NightPres = masterLun$NightPres[startInd]
#   #   GroupID = masterLun$GroupID[startInd]
#   #   MPhs = mSpline(masterLun$MoonPhase[startInd],
#   #                  knots=quantile(masterLun$MoonPhase[startInd],probs=c(0.275,0.5,0.725)),
#   #                  Boundary.knots=c(0,1),
#   #                  periodic=T)
#   #   MAs = mSpline(masterLun$MoonAltitude[startInd],
#   #                 knots=quantile(masterLun$MoonAltitude[startInd],probs=c(0.333,0.666)),
#   #                 Boundary.knots=c(min(masterLun$MoonAltitude[startInd]),max(masterLun$MoonAltitude[startInd])))
#   #   MPrF = as.factor(masterLun$MoonPres[startInd])
#   #   
#   # } else {
#   #   NightPres = masterLun$NightPres
#   #   GroupID = masterLun$GroupID
#   #   MPhs = mSpline(masterLun$MoonPhase,
#   #                  knots=quantile(masterLun$MoonPhase,probs=c(0.275,0.5,0.725)),
#   #                  Boundary.knots=c(0,1),
#   #                  periodic=T)
#   #   MAs = mSpline(masterLun$MoonAltitude,
#   #                 knots=quantile(masterLun$MoonAltitude,probs=c(0.333,0.666)),
#   #                 Boundary.knots=c(min(masterLun$MoonAltitude),max(masterLun$MoonAltitude)))
#   #   MPrF = as.factor(masterLun$MoonPres)
#   #   
#   # }
#   # 
#   # lunMod = geeglm(NightPres
#   #                 ~ MPhs
#   #                 + MAs
#   #                 + MPrF
#   #                 + MPhs:MAs
#   #                 + MPhs:MPrF,
#   #                 family=binomial(link="logit"),
#   #                 id=GroupID,
#   #                 corstr="ar1")
#   # 
#   # 
#   # PV = getPvalues(lunMod)
#   # if (sum(PV$'p-value'<0.05)==length(PV$'p-value')){ # if all variables are significant, save
#   #   if (lunMod$geese$error==0){ # if model converged, save model and summary
#   #     sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_LunPhAPr_5I.txt",sep="")
#   #     sink(sinkName)
#   #     print(PV)
#   #     sink()
#   #     save(lunMod,PV,file=paste(outDir,'/',CTname,'/',site,"_",int,"_Model_LunPhAPr_5I.Rdata",sep=""))
#   #   } else { # if model did not converge, don't save model, only summary indicating non-convergence
#   #     sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_LunPhAPr_5I.txt",sep="")
#   #     sink(sinkName)
#   #     print('Model did not converge')
#   #     sink()
#   #   }
#   #   
#   # } else if (sum(PV$'p-value'<0.05)==0){ # if no significant variables
#   #   sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_LunPhAPr_5I.txt",sep="")
#   #   sink(sinkName)
#   #   print('No significance')
#   #   sink()
#   # }
#   #   else {
#   #   PV2 = PV
#   #   while (sum(PV2$'p-value'<0.05) < length(PV2$'p-value')){
#   #     if (PV2$'p-value'[PV2$'Variable'=="MPhs:MAs"]<0.05){ # if interaction term is significant, keep both contributing terms in model
#   #       PV2$'p-value'[1:2] = 0.01
#   #     }
#   #     if (any(PV2$'Variable'=="MPhs:MPrF")){
#   #       if (PV2$'p-value'[PV2$'Variable'=="MPhs:MPrF"]<0.05){ # if interaction term is significant, keep both contributing terms in model
#   #         PV2$'p-value'[c(1,3)] = 0.01
#   #       }
#   #     }
#   #     whichVars = which(PV2$'p-value'<0.05)
#   #     str1 = "lunMod = geeglm(NightPres~"
#   #     str2 = ""
#   #     str3 = ",family=binomial,id=GroupID,corstr='ar1')"
#   #     for (j in 1:length(whichVars)){
#   #       if(j==1){
#   #         str2 = paste(str2,PV2$Variable[whichVars[j]],sep="")
#   #       } else {
#   #         str2 = paste(str2,PV2$Variable[whichVars[j]],sep="+")
#   #       }
#   #     }
#   #     modelCall = paste(str1,str2,str3,sep="")
#   #     eval(parse(text=modelCall))
#   #     PV = getPvalues(lunMod)
#   #     PV2 = PV
#   #     if (any(PV2$'Variable'=="MPhs:MAs")){
#   #       if (PV2$'p-value'[PV2$'Variable'=="MPhs:MAs"]<0.05){ # if interaction term is significant, keep both contributing terms in model
#   #         PV2$'p-value'[1:2] = 0.01
#   #       }
#   #     }
#   #     if (any(PV2$'Variable'=="MPhs:MPrF")){
#   #       if (PV2$'p-value'[PV2$'Variable'=="MPhs:MPrF"]<0.05){ # if interaction term is significant, keep both contributing terms in model
#   #         PV2$'p-value'[c(1,3)] = 0.01
#   #       }
#   #     }
#   #   }
#   #   if (lunMod$geese$error==0){ # if model converged, save model and summary
#   #     sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_LunPhAPr_5I.txt",sep="")
#   #     sink(sinkName)
#   #     print(PV)
#   #     sink()
#   #     save(lunMod,PV,file=paste(outDir,'/',CTname,'/',site,"_",int,"_Model_LunPhAPr_5I.Rdata",sep=""))
#   #   } else { # if model did not converge, don't save model, only summary indicating non-convergence
#   #     sinkName = paste(outDir,'/',CTname,'/',site,"_",int,"_GEEGLMSummary_LunPhAPr_5I.txt",sep="")
#   #     sink(sinkName)
#   #     print('Model did not converge')
#   #     sink()
#   #   }
#   # }
#   # 
# }
# end = Sys.time()
# end-start

# ##### CREATE MODELS USING BEST FIT SPLINES  --------------------------
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


### Old Code ---------------------------------------------------------------------------
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
