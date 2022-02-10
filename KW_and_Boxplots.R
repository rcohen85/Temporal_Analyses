library(pracma)
library(lubridate)
library(tidyverse)
library(conover.test)
library(dunn.test)
library(suncalc)
library(rstatix)
library(expss)
library(ggsignif)
library(dplyr)

# KW test & boxplots to compare months ------------------------------------
# 
# datDir = 'I:/TimeSeries_ScaledByEffortError'
# seasDir = 'I:/Seasonal_Plots'
# int = "Daily"
# 
# fileList = list.files(path=datDir,pattern=paste('*_',int,'.csv',sep=""),
#                       full.names=TRUE,recursive=FALSE,
#                       include.dirs=FALSE,no..=TRUE)
# goodIdx = c(1,2,4,8,11,13:19)
# 
# for (i in goodIdx){  # for each species' file
# 
#   thisCT = data.frame(read.csv(fileList[i]))  # load file
#   attribs = attributes(thisCT)  # find names of sites (all columns but first)
#   sites = attribs$names[2:numel(attribs$names)]
# 
#   CTname = str_remove(fileList[i],paste(datDir,'/',sep="")) # get the species/CT name
#   CTname = str_remove(CTname,paste('_',int,'.csv',sep=""))
# 
#   dateVec = as.Date(thisCT$Date,format="%d-%b-%Y")
# 
#   # if it doesn't already exist, create directory to save figures
#   if (!dir.exists(paste(seasDir,'/',CTname,sep=""))){
#     dir.create(paste(seasDir,'/',CTname,sep=""))
#   }
# 
#   for (j in 1:numel(sites)){ # for each site
# 
#     thisSite = as.numeric(thisCT[,j+1]) # get presence data for this site
#     noDat = which(is.na(thisSite)) # remove days with no effort
#     if (!numel(noDat) == 0){
#       thisSite = thisSite[-noDat]
#       reducedDateVec = dateVec[-noDat]
#       monthGroup = month(reducedDateVec) # find which month each day of data falls in
#     } else if (numel(noDat) == 0) {
#       reducedDateVec = dateVec
#       monthGroup = month(reducedDateVec) # find which month each day of data falls in
#     }
# 
# 
#     pres = which(thisSite!=0)
# 
#     if (!numel(pres)==0){ # if there is any presence
# 
#       # test for significant differences between months
#       kruskal.test(thisSite,monthGroup,na.action=na.pass)
#       pairCompCI = conover.test(thisSite,monthGroup,method='bonferroni',label=TRUE,
#                                 wrap=TRUE,table=TRUE,alpha=0.05)
#       # pairCompDunn = dunn.test(thisSite,monthGroup,method='bonferroni',wrap=TRUE,
#       #                      table=TRUE,alpha=0.05)
#       CITable = matrix(NA,11,11)
#       CITable[upper.tri(CITable,diag=TRUE)] = t(pairCompCI$P.adjusted)
#       CITable = t(CITable)
#       rownames(CITable) = c("10","11","12","2","3","4","5","6","7","8","9")
#       colnames(CITable) = c("1","10","11","12","2","3","4","5","6","7","8")
#       saveName = paste(seasDir,'/',CTname,'/',sites[j],"_Monthly_CI_Comparisons.csv",sep="")
#       write.csv(CITable,file=saveName,row.names=TRUE)
# 
#       # Identify which pairwise comparisons were significant
#       # sig = which(pairCompCI$P.adjusted<0.025)
#       # sigPairs = pairCompCI$comparisons[sig]
#       # sigNumbers = regmatches(sigPairs, gregexpr("[[:digit:]]+", sigPairs))
#       # sigNumbers = as.numeric(unlist(sigNumbers))
# 
#       # Plot data by month
#       saveName = paste(seasDir,'/',CTname,'/',sites[j],"_MonthlyBoxplot.png",sep="")
#       png(saveName)
#       boxplot(thisSite~monthGroup,main=(paste(CTname,'at',sites[j])),
#               xlab=c("Month"),ylab=c("Daily Counts"))
#       while (dev.cur()>1) {dev.off()}
# 
#     }
# 
#   }
# 
# }
# 

# KW test & boxplots to compare years, phases of the day -------------------------

datDir = 'J:/Chpt_2/TimeSeries_ScaledByEffortError'
illumDir = 'J:/Chpt_2/IlluminationFiles'
normTimeDir = 'J:/Chpt_2/NormTimes'
dielDir = 'J:/Chpt_2/Diel Plots'
int = "5minBin"
phases = 4 # 4 for Dawn/Day/Dusk/Night or 2 for Day/Night
type = "B" # "V" for violin plots or "B" for boxplots
goodIdx = c(1,2,4,8,11,13:19)
goodSite = list(c(8:10),c(8:10),c(),c(1:8),c(),c(),c(),c(4,6:10),c(),c(),c(1:11),
                c(),c(1:5),c(1:11),c(1:5),
                c(1:11),c(1:8),c(1:11),c(1:7))
# goodIdx = c(14)
# goodSite = list(c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),
#                 c(),c(),c(1:11),c(),
#                 c(),c(),c(),c())


dateStart = as.POSIXct('2016-05-01 00:00:00',format="%Y-%m-%d %H:%M:%S",tz="GMT");
dateEnd = as.POSIXct('2019-04-30 23:59:59',format="%Y-%m-%d %H:%M:%S",tz="GMT");
lats = as.numeric(c(41.06165,40.22999,39.83295,
                    39.19192,38.37337,37.16452,
                    35.30183,33.66992,32.10527,
                    30.58295,30.27818))
lons = as.numeric(c(-66.35155,-67.97798,-69.98194,
                    -72.22735,-73.36985,-74.46585,
                    -74.87895,-75.9977,-77.09067,
                    -77.39002,-80.22085))

datFileList = list.files(path=datDir,pattern=paste('*_',int,'.csv',sep=""),
                         full.names=TRUE,recursive=FALSE,
                         include.dirs=FALSE,no..=TRUE)
effFileList = list.files(path=datDir,pattern=paste('*_',int,'_Effort.csv',sep=""),
                         full.names=TRUE,recursive=FALSE,
                         include.dirs=FALSE,no..=TRUE)
normTimeFileList = list.files(path=normTimeDir,pattern='*.csv',
                              full.names=TRUE,recursive=FALSE,
                              include.dirs=FALSE,no..=TRUE)


for (i in goodIdx){  # for each species' file

  thisCT = data.frame(read.csv(datFileList[i]))  # load file
  attribs = attributes(thisCT)  # find names of sites (all columns but first)
  sites = attribs$names[2:numel(attribs$names)]

  CTname = str_remove(datFileList[i],paste(datDir,'/',sep="")) # get the species/CT name
  CTname = str_remove(CTname,paste('_',int,'.csv',sep=""))
  if (str_detect(CTname,"Atl")){
    CTname = "Gervais"
  }

  binVec = as.POSIXct(thisCT$Bin,format="%d-%b-%Y %H:%M:%S", tz="GMT")

  # if it doesn't already exist, create directory to save figures
  if (!dir.exists(paste(dielDir,'/',CTname,sep=""))){
    dir.create(paste(dielDir,'/',CTname,sep=""))
  }

  for (j in unlist(goodSite[i])){ # for each site

    if (sites[j]=="HAT"){ # for HAT, only look at site B
      goodInd = which(binVec>=as.POSIXct('2017-05-01 00:00:00',format="%Y-%m-%d %H:%M:%S",tz="GMT"))
      thiCT = thisCT[goodInd,]
    }
    thisSite = as.numeric(thisCT[,j+1]) # get presence data for this site
    noDat = which(is.na(thisSite)) # remove bins with no effort
    if (!numel(noDat) == 0){
      thisSite = thisSite[-noDat]
      reducedBinVec = binVec[-noDat]
      yearGroup = year(reducedBinVec) # find year each bin falls in
    } else if (numel(noDat) == 0) {
      reducedBinVec = binVec
      yearGroup = year(reducedBinVec) # find year each bin falls in
    }

    # presBins = thisSite!=0
    if (any(i==c(11,16:19))){
      presBins = thisSite>=50 # consider dolphins "present" if more than 50 clicks in bin
    } else if (any(i==c(1,2,4,8,13:15))){
      presBins = thisSite>=20 # consider BWs/Kogia/Pm "present" if more than 20 clicks in bin
    }
    pres = which(presBins)

    if (!numel(pres)==0){ # if there is any presence

      # YEAR ANALYSIS ---------------------
      # # quantify effort in each year
      # effInd = which(!is.na(str_match(effFileList,sites[j])))
      # effort = data.frame(read.csv(effFileList[effInd]))
      # effDate = as.POSIXlt(effort$Bin,format="%d-%b-%Y %H:%M:%S", tz="GMT")
      # effYear = year(effDate)
      #
      # # sum presence in each year & normalize by recording effort in that year
      # yearly_pres = aggregate(as.numeric(presBins),list(yearGroup),FUN=sum)
      # yearly_eff = aggregate(effort$Effort,list(effYear),FUN=sum)
      # yearly_eff$Norm = yearly_eff$x/max(yearly_eff$x)
      # yearly_pres$Norm = yearly_pres$x*(1/yearly_eff$Norm)
      #
      # YrP = ggplot(yearly_pres,
      #       )+geom_col(aes(Group.1,Norm)
      #       )+labs(x="",y="Number of 5-min Bins with Presence")


      # DAY/NIGHT/DAWN/DUSK ANALYSIS ---------------------
      if (phases==4){
        # Get day phase data
        dayData = getSunlightTimes(date=seq.Date(as.Date(dateStart),as.Date(dateEnd),by=1),
                                   lat=lats[j],lon=lons[j],
                                   keep=c("nauticalDawn","sunrise","sunset","nauticalDusk"),
                                   tz="UTC")
        dayData = dayData[,4:7]
        phaseBins = as.POSIXct(c(t(dayData)),format="%Y-%m-%d %H:%M:%S",tz="GMT")
        tooLate = which(phaseBins>=dateEnd)
        phaseBins = phaseBins[-tooLate] # cutting out some good bins??
        # phaseBins = c(dateStart,phaseBins,dateEnd)
        # phaseVec = unlist(rep(list("Night","Dawn","Day","Dusk"),length.out=length(phaseBins)-1))
        phaseVec = unlist(rep(list("Dawn","Day","Dusk","Night"),length.out=length(phaseBins)-1))

        # find which day phase each presence bin falls into
        whichBin = histc(as.numeric(reducedBinVec[pres]),as.numeric(phaseBins))
        whichBin$cnt = whichBin$cnt[-length(whichBin$cnt)]

        # normalize counts to phase length
        phaseDurs = diff(as.numeric(phaseBins),lag=1) # length of each phase in seconds
        normPhaseLength = phaseDurs/max(phaseDurs)
        normCount = (whichBin$cnt)*(1/normPhaseLength)
        # normCount = (whichBin$cnt/(60/5))*(1/phaseDurs)

        # # quantify proportion of each phase with clicking
        # phaseDurs = diff(as.numeric(phaseBins),lag=1)/(60*5) # length of each phase bin in 5-min bins
        # clickProp = whichBin$cnt/phaseDurs

        # find phases with presence
        # presPhases = which(clickProp>0)
        presPhases = which(normCount>0)

        DielDF = data.frame(normCount[presPhases],as.factor(phaseVec[presPhases]))
        colnames(DielDF) = c("NormPres","DayPhase")
        # DielDF = data.frame(clickProp[presPhases],as.factor(phaseVec[presPhases]))
        # colnames(DielDF) = c("ClickProp","DayPhase")

        # Plot click prop each day phase & test for differences
        if (length(levels(DielDF$DayPhase))>1){ # if there is presence in more than one day phase
          # Test significance of differences btwn day phases
          KW = kruskal.test(DielDF$NormPres,DielDF$DayPhase,na.action=na.pass)
          pairCompCI = conover.test(DielDF$NormPres,DielDF$DayPhase,method='bonferroni',label=TRUE,
                                    wrap=TRUE,table=TRUE,alpha=0.05)

          ndim=numel(levels(DielDF$DayPhase))-1
          CITable = matrix(NA,ndim,ndim)
          CITable[upper.tri(CITable,diag=TRUE)] = t(pairCompCI$P.adjusted)
          CITable = t(CITable)
          CITable = signif(CITable,digits=3)
          Comparisons = matrix(NA,ndim,ndim)
          Comparisons[upper.tri(Comparisons,diag=TRUE)] = t(pairCompCI$comparisons)
          Comparisons = t(Comparisons)
          CITable = rbind(CITable,Comparisons)
          saveName = paste(dielDir,'/',CTname,'/minClicks50/',sites[j],"_4Phase_CI_Comparisons.csv",sep="")
          write.csv(CITable,file=saveName,row.names=TRUE)

          # calculate quantiles for plotting limits
          quants = DielDF %>%
            group_by(DayPhase) %>%
            summarize(q25 = quantile(NormPres,probs=0.25),
                      q50 = quantile(NormPres,probs=0.50),
                      q75 = quantile(NormPres,probs=0.75),
                      q90 = quantile(NormPres,probs=0.9))
          iqr = quants$q75-quants$q25

          # Find and remove outliers for plotting (only if doing violin plots)
          if (type=="V"){
            out = DielDF %>%
              group_by(DayPhase) %>%
              identify_outliers("NormPres")
            if (numel(out!=0)){
              out = out[,-c(3,4)]
              outInd = which(!is.na(match_row(out,DielDF)))
              DielDF = DielDF[-outInd,]
            }
          }

          my_comparisons <- list( c("Dawn", "Day"), c("Dawn", "Dusk"),c("Day","Dusk"),
                                  c("Dawn", "Night"), c("Day","Night"),c("Dusk","Night"))
          sigComps = which(pairCompCI$P.adjusted<0.025)
          nonSigComps = which(pairCompCI$P.adjusted>=0.025)

          # Plot presence per day phase
          if(type=="B"){ # boxplot for if there are clicks in more than one day phase
            ggplot(DielDF,aes(DayPhase,NormPres)
            ) + geom_boxplot(varwidth=FALSE,
                            outlier.shape=NA
            ) + scale_y_continuous(limits=c(0,(3.2*max(iqr))+max(quants$q75))
            ) + scale_x_discrete(limits=c("Dawn","Day","Dusk","Night") # make sure NaNs aren't displayed, if they remain in the data
            ) + labs(title = paste(CTname, 'at',sites[j]),
                     x=(""),y="Proportion of Phase with Presence"
            ) + geom_signif(comparisons=my_comparisons[sigComps],
                            y_position=(max(iqr)*c(1.5,1.75,2,2.25,2.5,2.75))+max(quants$q75),
                            test=wilcox.test,
                            textsize=6,
                            # step_increase=0.07,
                            tip_length=0.01,
                            vjust=0.5,
                            # annotations = sprintf("p = %.3g",pairCompCI$P.adjusted[sigComps]))
                            map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.025))


            saveName = paste(dielDir,'/',CTname,'/minClicks50/',sites[j],"_4PhaseBoxplot_normPres.png",sep="")
            ggsave(saveName,device="png")

          } else if (type=="V") { # violin plot for if there are clicks in more than one day phase
            ggplot(DielDF,aes(DayPhase,NormPres)
            ) + geom_violin(scale="count",
                            draw_quantiles=c(0.25, 0.5, 0.75)
            ) + scale_x_discrete(limits=c("Dawn","Day","Dusk","Night") # make sure NaNs aren't displayed, if they remain in the data
            ) + labs(title = paste(CTname, 'at',sites[j]),
                     x=(""),y="Proportion of Phase with Presence"
            ) + geom_signif(comparisons=my_comparisons[sigComps],
                            y_position=(max(iqr)*c(1.5,1.7,1.9))+max(quants$q75),
                            test=wilcox.test,
                            # step_increase=0.07,
                            tip_length=0.01,
                            annotations = sprintf("p = %.3g",pairCompCI$P.adjusted))
                            # map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.025))

            saveName = paste(dielDir,'/',CTname,'/minClicks50/',sites[j],"_4PhaseViolinplot_normPres.png",sep="")
            ggsave(saveName,device="png")
          }

        } else if (length(levels(DielDF$DayPhase))==1) {
          if (type=="B"){ # boxplot for clicks in only one day phase
            ggplot(DielDF,aes(DayPhase,NormPres)
            ) + geom_boxplot(varwidth=FALSE,
                             outlier.shape=NA
            ) + scale_y_continuous(limits=c(0,(1.65*max(iqr))+max(quants$q75))
            ) + annotate("text",size=8,
                         x=levels(DielDF$DayPhase),y=1.05*max(DielDF$NormPres),
                         label=paste("Presence only during ",levels(DielDF$DayPhase),sep="")
            ) + labs(title = paste(CTname, 'at',sites[j]),
                     x=(""),y="Proportion of Phase with Presence")

            saveName = paste(dielDir,'/',CTname,'/',sites[j],"_4PhaseBoxplot_normPres.png",sep="")
            ggsave(saveName,device="png")

          } else if (type=="V") { # violin plot for clicks in only one day phase
            ggplot(DielDF,aes(DayPhase,NormPres)
            ) + geom_violin(scale="count",
                            draw_quantiles=c(0.25, 0.5, 0.25)
            ) + annotate("text",size=8,
                         x=levels(DielDF$DayPhase),y=1.05*max(DielDF$NormPres),
                         label=paste("Presence only during ",levels(DielDF$DayPhase),sep="")
            ) + labs(title = paste(CTname, 'at',sites[j]),
                     x=(""),y="Proportion of Phase with Presence")
            saveName = paste(dielDir,'/',CTname,'/',sites[j],"_4PhaseViolinplot_normPres.png",sep="")
            ggsave(saveName,device="png")
          }
        }
      }

      # DAY/NIGHT ANALYSIS -------------------------
      if (phases==2){

        # # Load normalized bin times
        # whichNormTimeFile = str_detect(normTimeFileList,sites[j])
        # normTimes = data.frame(read.csv(normTimeFileList[whichNormTimeFile]))  # load file
        # normBinTimes = normTimes$Bin[-noDat]
        #
        # # create vector coding for day/night
        # dayNight = rep(0,length(normBinTimes))
        # dayNight[normBinTimes>0] = "Night"
        # dayNight[normBinTimes<0] = "Day"
        # # find and remove presence which falls right on sunrise/sunset
        # neither = which(dayNight[pres]==0)
        # if (length(neither)!=0){
        #   keep = pres[-neither]
        # } else { keep = pres }

        # Get day phase data
        dayData = getSunlightTimes(date=seq.Date(as.Date(dateStart),as.Date(dateEnd),by=1),
                                   lat=lats[j],lon=lons[j],
                                   keep=c("sunrise","sunset"),
                                   tz="UTC")
        dayData = dayData[,4:5]
        phaseBins = as.POSIXlt(c(t(dayData)),format="%Y-%m-%d %H:%M:%S",tz="GMT")
        tooLate = which(phaseBins>=dateEnd)
        if (numel(tooLate)!=0){
        phaseBins = phaseBins[-tooLate]} # cutting out some good bins??
        phaseBins = c(dateStart,phaseBins,dateEnd)
        phaseVec = unlist(rep(list("Night","Day"),length.out=length(phaseBins)-1))

        # find which day phase each presence bin falls into
        whichBin = histc(as.numeric(reducedBinVec[pres]),as.numeric(phaseBins))
        whichBin$cnt = whichBin$cnt[-length(whichBin$cnt)]

        # quantify proportion of each phase with clicking
        phaseDurs = diff(as.numeric(phaseBins),lag=1)/(60*5) # length of each phase bins
        clickProp = whichBin$cnt/phaseDurs

        # find phases with presence
        presPhases = which(clickProp>0)

        # DielDF = data.frame(thisSite[keep],as.factor(dayNight[keep]))
        # colnames(DielDF) = c("Presence","DayNight")
        DielDF = data.frame(clickProp[presPhases],as.factor(phaseVec[presPhases]))
        colnames(DielDF) = c("ClickProp","DayPhase")

        # calculate quantiles for plotting limits
        quants = DielDF %>%
          group_by(DayPhase) %>%
          summarize(q25 = quantile(ClickProp,probs=0.25),
                    q50 = quantile(ClickProp,probs=0.50),
                    q75 = quantile(ClickProp,probs=0.75),
                    q90 = quantile(ClickProp,probs=0.9))
        iqr = quants$q75 - quants$q25

        # # Find and remove outliers for plotting (only for violin plots)
        # if (type=="V") {
        #   out = DielDF %>%
        #     group_by(DayPhase) %>%
        #     identify_outliers("Presence")
        #   out = out[,-c(3,4)]
        #   outInd = which(!is.na(match_row(out,DielDF)))
        #   DielDF = DielDF[-outInd,]
        # }

        # Plot click prop in day/night & test for differences
        if (length(levels(DielDF$DayPhase))==2){ # if there are clicks in more than one day phase
          # Test significance of differences btwn day & night clicking
          KW = kruskal.test(DielDF$ClickProp,DielDF$DayPhase,na.action=na.pass)

          # Plot day vs night presence
          if (type=="B"){ # boxplot
            ggplot(DielDF,aes(DayPhase,ClickProp)
            ) + geom_boxplot(varwidth=FALSE,
                             outlier.shape=NA
            ) + scale_y_continuous(limits=c(0,(1.65*max(iqr))+max(quants$q75))
            ) + annotate("text",size=8,
                         x=levels(DielDF$DayPhase)[1],y=(1.6*max(iqr))+max(quants$q75),
                         label=paste("P Value: ",signif(KW$p.value,digits=3),sep="")
            ) + labs(title = paste(CTname, 'at',sites[j]),
                     x=(""),y="Proportion of Phase with Presence")
            saveName = paste(dielDir,'/',CTname,'/minClicks50/',sites[j],"_DielBoxplot_prop.png",sep="")
          } else if (type=="V") { # violin plot
            ggplot(DielDF,aes(DayPhase,ClickProp)
            ) + geom_violin(scale="count",
                            draw_quantiles=c(0.25, 0.5, 0.25)
            ) + annotate("text",size=8,
                         x=levels(DielDF$DayPhase)[1],y=1.2*max(quants$q75),
                         label=paste("P Value: ",signif(KW$p.value,digits=3),sep="")
            ) + labs(title = paste(CTname, 'at',sites[j]),
                     x=(""),y="Proportion of Phase with Presence")
            saveName = paste(dielDir,'/',CTname,'/minClicks50/',sites[j],"_DielViolinplot_prop.png",sep="")
          }

          ggsave(saveName,device="png")

        } else { # if there are clicks in only one day phase, just plot
          if (type=="B"){ # boxplot
            ggplot(DielDF,aes(DayPhase,ClickProp)
            ) + geom_boxplot(varwidth=FALSE,
                             outlier.shape=NA
            ) + annotate("text",size=8,
                         x=levels(DielDF$DayPhase),y=(1.6*max(iqr))+max(quants$q75),
                         label=paste("Presence only during ",levels(DielDF$DayPhase),sep="")
            ) + labs(title = paste(CTname, 'at',sites[j]),
                     x=(""),y="Proportion of Phase with Presence")
            saveName = paste(dielDir,'/',CTname,'/minClicks50/',sites[j],"_DielBoxplot_prop.png",sep="")
          } else if (type=="V") { # violin plot
            ggplot(DielDF,aes(DayPhase,ClickProp)
            ) + geom_violin(scale="count",
                            draw_quantiles=c(0.25, 0.5, 0.25)
            ) + annotate("text",size=8,
                         x=levels(DielDF$DayPhase),y=1.05*max(DielDF$ClickProp),
                         label=paste("Presence only during ",levels(DielDF$DayPhase),sep="")
            ) + labs(title = paste(CTname, 'at',sites[j]),
                     x=(""),y="Proportion of Phase with Presence")
            saveName = paste(dielDir,'/',CTname,'/minClicks50/',sites[j],"_DielViolinplot_prop.png",sep="")
          }

          ggsave(saveName,device="png")
        }
      }
    }
  }
}





# # KW Tests and Boxplots to Compare Before/During/After Moon-Up Period -------------------
# 
# datDir = 'J:/Chpt_2/TimeSeries_ScaledByEffortError'
# dielDir = 'J:/Chpt_2/Diel Plots'
# int = "5minBin"
# dateStart = as.POSIXct('2016-05-01 00:00:00',format="%Y-%m-%d %H:%M:%S",tz="GMT");
# dateEnd = as.POSIXct('2019-04-30 23:59:59',format="%Y-%m-%d %H:%M:%S",tz="GMT");
# sites = c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS','JAX')
# lats = as.numeric(c(41.06165,40.22999,39.83295,
#                     39.19192,38.37337,37.16452,
#                     35.30183,33.66992,32.10527,
#                     30.58295,30.27818))
# lons = as.numeric(c(-66.35155,-67.97798,-69.98194,
#                     -72.22735,-73.36985,-74.46585,
#                     -74.87895,-75.9977,-77.09067,
#                     -77.39002,-80.22085))
# 
# datFileList = list.files(path=datDir,pattern=paste('*_',int,'_MasterLun.csv',sep=""),
#                          full.names=TRUE,recursive=FALSE,
#                          include.dirs=FALSE,no..=TRUE)
# 
# for (i in 1:length(datFileList)){
# 
#   masterLun = data.frame(read.csv(datFileList[i]))  # load file
#   masterLun$NightBinTimes = as.POSIXct(masterLun$NightBinTimes,format="%Y-%m-%d %H:%M:%S",tz="GMT")
#   CTname = str_remove(datFileList[i],paste(datDir,'/',sep="")) # get the species/CT name
#   site = str_remove(CTname,paste("_",int,"_MasterLun.csv",sep=""))
#   site = sub(".*_","",site)
#   CTname = sub("_.*","",CTname)
#   if (str_detect(CTname,"Atl")){
#     CTname = "Gervais"
#   }
# 
#   thislat = lats[which(str_detect(sites,site))]
#   thislon = lons[which(str_detect(sites,site))]
# 
#   if (site=="HAT"){
#     goodInd = which(masterLun$NightBinTimes>=as.POSIXct('2017-05-01 00:00:00',format="%Y-%m-%d %H:%M:%S",tz="GMT"))
#     masterLun = masterLun[goodInd,]
#   }
# 
#   pres = which(masterLun$NightPres>0)
# 
#   # get start/end of each night
#   dayData = getSunlightTimes(date=seq.Date(as.Date(dateStart),as.Date(dateEnd),by=1),
#                              lat=thislat,lon=thislon,
#                              keep=c("sunrise","sunset"),
#                              tz="UTC")
#   sunR = dayData$sunrise[!is.na(dayData$sunrise)]
#   sunS = dayData$sunset[!is.na(dayData$sunset)]
# 
#   # get moon rise/set times
#   moonData = getMoonTimes(date=seq.Date(as.Date(dateStart),as.Date(dateEnd+(60*60*24)),by=1),
#                           lat=thislat,lon=thislon,
#                           keep=c("rise","set"),
#                           tz="UTC")
# 
#   moonR = moonData$rise[!is.na(moonData$rise)]
#   moonS = moonData$set[!is.na(moonData$set)]
# 
#   # name each nighttime moon bin (pre-moon, moon up, post-moon)
#   moonBins = matrix(nrow=length(c(sunS,sunR,moonS,moonR)),ncol=4)
#   ind = 1
#   for (j in 1:length(sunS)){ # for each night
# 
#     # find the subsequent sunrise
#     sunRind = which(sunR>sunS[j])[1]
# 
#     if (!is.na(sunRind)){ # if we have start & end info for this night
# 
#       # find moonrise and/or moonset on this night
#       moonRind = which(moonR>=sunS[j] & moonR<sunR[sunRind])
#       moonSind = which(moonS>=sunS[j] & moonS<sunR[sunRind])
# 
#       moonBins[ind,1] = sunS[j]
#       moonBins[ind,2] = "SunS"
# 
#       if (!is_empty(moonRind) & is_empty(moonSind)){ # if the moon rises this night, but doesn't set
# 
#         moonBins[(ind+1):(ind+2),1] = c(moonR[moonRind][1],sunR[sunRind])
#         moonBins[(ind+1):(ind+2),2] = c("MoonR","SunR")
#         moonBins[ind:(ind+2),3] = c("Pre","MoonUp","Day")
# 
#         ind = ind+3
# 
#       } else if (is_empty(moonRind) & !is_empty(moonSind)){ # if the moon sets tonight, but doesn't rise
# 
#         moonBins[(ind+1):(ind+2),1] = c(moonS[moonSind],sunR[sunRind])
#         moonBins[(ind+1):(ind+2),2] = c("MoonS","SunR")
#         moonBins[ind:(ind+2),3] = c("MoonUp","Post","Day")
#         ind = ind+3
# 
#       } else if (is_empty(moonRind) & is_empty(moonSind)){ # if moon doesn't rise or set this night
# 
#         moonBins[ind+1,1] = sunR[sunRind] # just add sunrise info
#         moonBins[ind+1,2] = "SunR"
#         moonBins[ind:(ind+1),3] = c("Pre","Day")
#         ind = ind+2
# 
#       } else if (!is_empty(moonRind) & !is_empty(moonSind)){ # if the moon both rises and sets this night
# 
#         whichFirst = which.min(c(moonR[moonRind],moonS[moonSind]))
#         if (whichFirst==1){
# 
#           moonBins[(ind+1):(ind+3),1] = c(moonR[moonRind],moonS[moonSind],sunR[sunRind])
#           moonBins[(ind+1):(ind+3),2] = c("MoonR","MoonS","SunR")
#           moonBins[ind:(ind+2),3] = c("Pre","MoonUp","Post","Day")
#           ind = ind+4
# 
#         } else if (whichFirst==2){
# 
#           moonBins[(ind+1):(ind+3),1] = c(moonS[moonSind],moonR[moonRind],sunR[sunRind])
#           moonBins[(ind+1):(ind+3),2] = c("MoonS","MoonR","SunR")
#           moonBins[ind:(ind+3),3] = c("MoonUp","Pre","MoonUp","Day")
#           ind = ind+4
# 
#         }
# 
#       }
#     }
#   }
# 
#   # get rid of extra rows
#   moonBins = moonBins[which(!is.na(moonBins[,1])),]
# 
#   # calculate duration of each moon pres bin, in units of 5-min bins
#   moonBins[1:dim(moonBins)[1]-1,4] = diff(as.numeric(moonBins[,1]),lag=1)/(60*5)
# 
#   # find which moon bin each nighttime presence bin falls into
#   whichBin = histc(as.numeric(masterLun$NightBinTimes[pres]),as.numeric(moonBins[,1]))
#   whichBin$cnt = whichBin$cnt[-length(whichBin$cnt)]
# 
#   # quantify proportion of each moon bin with clicking
#   clickProp = whichBin$cnt/as.numeric(moonBins[dim(moonBins)[1]-1,4])
# 
#   # find moon bins with presence
#   presPhases = which(clickProp>0)
# 
#   LunDF = data.frame(clickProp[presPhases],as.factor(moonBins[presPhases,3]))
#   colnames(LunDF) = c("ClickProp","MoonPres")
# 
#   # calculate quantiles for plotting limits
#   quants = LunDF %>%
#     group_by(MoonPres) %>%
#     summarize(q25 = quantile(ClickProp,probs=0.25),
#               q50 = quantile(ClickProp,probs=0.50),
#               q75 = quantile(ClickProp,probs=0.75),
#               q90 = quantile(ClickProp,probs=0.9))
#   iqr = quants$q75 - quants$q25
# 
#   # Plot # clicks per each moon pres bin & test for differences
# 
#   if (length(levels(LunDF$MoonPres))>1){ # if there are clicks in more than one moon pres state
# 
#     # Test significance of differences btwn moon pres states
#     KW = kruskal.test(LunDF$ClickProp,LunDF$MoonPres,na.action=na.pass)
# 
#     if (length(levels(LunDF$MoonPres))==2){
#       # if only comparing 2 moon pres states, save KW output
#       sinkName = paste(dielDir,'/',CTname,'/minClicks50/',site,"_MoonPres_CI_Comparisons.txt",sep="")
#       sink(sinkName)
#       print(KW)
#       sink()
# 
#       # Plot click prop per moon pres state with p-value shown
#       ggplot(LunDF,aes(MoonPres,ClickProp)
#       ) + geom_boxplot(varwidth=FALSE,
#                        outlier.shape=NA
#       ) + scale_y_continuous(limits=c(0,(1.65*max(iqr))+max(quants$q75))
#       ) + annotate("text",size=8,
#                    x=levels(LunDF$MoonPres)[1],y=(1.6*max(iqr))+max(quants$q75),
#                    label=paste("P Value: ",signif(KW$p.value,digits=3),sep="")
#       ) + labs(title = paste(CTname, 'at',site),
#                x=(""),y="Proportion of Phase with Presence")
# 
#       saveName = paste(dielDir,'/',CTname,'/minClicks50/',site,"_MoonPresBoxplot_prop.pdf",sep="")
#       ggsave(saveName,device="pdf")
#       saveName = paste(dielDir,'/',CTname,'/minClicks50/',site,"_MoonPresBoxplot_prop.png",sep="")
#       ggsave(saveName,device="png")
# 
#     } else if (length(levels(LunDF$MoonPres))==3) {
# 
#       # if comparing 3 moon pres states, carry out pairwise comparisons
#       pairCompCI = conover.test(LunDF$ClickProp,LunDF$MoonPres,method='bonferroni',label=TRUE,
#                                 wrap=TRUE,table=TRUE,alpha=0.05)
# 
#       ndim=numel(levels(LunDF$MoonPres))-1
#       CITable = matrix(NA,ndim,ndim)
#       CITable[upper.tri(CITable,diag=TRUE)] = t(pairCompCI$P.adjusted)
#       CITable = t(CITable)
#       CITable = signif(CITable,digits=3)
#       Comparisons = matrix(NA,ndim,ndim)
#       Comparisons[upper.tri(Comparisons,diag=TRUE)] = t(pairCompCI$comparisons)
#       Comparisons = t(Comparisons)
#       CITable = rbind(CITable,Comparisons)
#       saveName = paste(dielDir,'/',CTname,'/minClicks50/',site,"_MoonPres_CI_Comparisons.csv",sep="")
#       write.csv(CITable,file=saveName,row.names=TRUE)
# 
#       my_comparisons <- list( c("MoonUp","Post"),c("MoonUp","Pre"), c("Post","Pre"))
#       # sigComps = which(pairCompCI$P.adjusted<0.025)
#       # nonSigComps = which(pairCompCI$P.adjusted>=0.025)
# 
#       # Plot click prop per moon pres state with brackets indicating significance
#       ggplot(LunDF,aes(MoonPres,ClickProp)
#       ) + geom_boxplot(varwidth=FALSE,
#                        outlier.shape=NA
#       ) + scale_y_continuous(limits=c(0,(2.5*max(iqr))+max(quants$q75))
#       ) + scale_x_discrete(limits=c("Pre","MoonUp","Post") # make sure NaNs aren't displayed, if they remain in the data
#       ) + labs(title = paste(CTname, 'at',site),
#                x=(""),y="Proportion of Phase with Presence"
#       ) + geom_signif(comparisons=my_comparisons,
#                       y_position=(max(iqr)*c(1.5,1.7,1.9))+max(quants$q75),
#                       test=wilcox.test,
#                       # step_increase=0.03,
#                       tip_length = 0.01,
#                       annotations = sprintf("p = %.3g",pairCompCI$P.adjusted))
#                       # map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.025))
# 
#       saveName = paste(dielDir,'/',CTname,'/minClicks50/',site,"_MoonPresBoxplot_prop.pdf",sep="")
#       ggsave(saveName,device="pdf")
#       saveName = paste(dielDir,'/',CTname,'/minClicks50/',site,"_MoonPresBoxplot_prop.png",sep="")
#       ggsave(saveName,device="png")
# 
#     }
# 
# 
#   } else if (length(levels(LunDF$MoonPres))==1) { # if there are clicks in only 1 moon pres state
#     # just plot variability of click prop in that state
#     ggplot(LunDF,aes(MoonPres,ClickProp)
#     ) + geom_boxplot(varwidth=FALSE,
#                      outlier.shape=NA
#     ) + scale_y_continuous(limits=c(0,(1.65*max(iqr))+max(quants$q75))
#     ) + annotate("text",size=8
#                  x=levels(LunDF$MoonPres),y=1.05*max(LunDF$ClickProp),
#                  label=paste("Presence only during ",levels(LunDF$MoonPres),sep="")
#     ) + labs(title = paste(CTname, 'at',site),
#              x=(""),y="Proportion of Phase with Presence")
# 
#     saveName = paste(dielDir,'/',CTname,'/minClicks50/',site,"_MoonPresBoxplot_prop.pdf",sep="")
#     ggsave(saveName,device="pdf")
#     saveName = paste(dielDir,'/',CTname,'/minClicks50/',site,"_MoonPresBoxplot_prop.png",sep="")
#     ggsave(saveName,device="png")
# 
# 
#   }
# }
# 
# 
# 
# 
# 
# 
# 
