library(pracma)
library(tidyverse)
library(conover.test)
library(dunn.test)
library(suncalc)
library(rstatix)
library(expss)
library(ggsignif)

# KW test & boxplots to compare months ------------------------------------

# datDir = 'I:/TimeSeries'
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
#   CTname = str_remove(fileList[i],paste(inDir,'/',sep="")) # get the species/CT name
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
#       dev.off()
#       
#     }
# 
#   }
#   
# }



# KW test & boxplots to compare phases of the day -------------------------

datDir = 'I:/TimeSeries_ScaledByEffortError'
illumDir = 'I:/IlluminationFiles'
normTimeDir = 'I:/NormTimes'
dielDir = 'I:/Diel Plots'
int = "5minBin"
phases = 4 # 4 for Dawn/Day/Dusk/Night or 2 for Day/Night
type = "V" # "V" for violin plots or "B" for boxplots
goodIdx = c(1,2,4,8,11,13:19)

dateStart = as.POSIXlt('2016-05-01 00:00:00',format="%Y-%m-%d %H:%M:%S",tz="GMT");
dateEnd = as.POSIXlt('2019-04-30 23:59:59',format="%Y-%m-%d %H:%M:%S",tz="GMT");
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
normTimeFileList = list.files(path=normTimeDir,pattern='*.csv',
                              full.names=TRUE,recursive=FALSE,
                              include.dirs=FALSE,no..=TRUE)


for (i in goodIdx){  # for each species' file
  
  thisCT = data.frame(read.csv(datFileList[i]))  # load file
  attribs = attributes(thisCT)  # find names of sites (all columns but first)
  sites = attribs$names[2:numel(attribs$names)]
  
  CTname = str_remove(datFileList[i],paste(datDir,'/',sep="")) # get the species/CT name
  CTname = str_remove(CTname,paste('_',int,'.csv',sep=""))
  
  binVec = as.POSIXlt(thisCT$Bin,format="%d-%b-%Y %H:%M:%S", tz="GMT")
  
  # if it doesn't already exist, create directory to save figures
  if (!dir.exists(paste(dielDir,'/',CTname,sep=""))){
    dir.create(paste(dielDir,'/',CTname,sep=""))
  }
  
  for (j in 1:numel(sites)){ # for each site
    
    thisSite = as.numeric(thisCT[,j+1]) # get presence data for this site
    noDat = which(is.na(thisSite)) # remove bins with no effort
    if (!numel(noDat) == 0){
      thisSite = thisSite[-noDat]
      reducedBinVec = binVec[-noDat]
    } else if (numel(noDat) == 0) {
      reducedBinVec = binVec
    }
    
    pres = which(thisSite!=0)
    
    if (!numel(pres)==0){ # if there is any presence
      if (phases==4){  # DAY/NIGHT/DAWN/DUSK ANALYSIS
        # Get day phase data
        dayData = getSunlightTimes(date=seq.Date(as.Date(dateStart),as.Date(dateEnd),by=1),
                                   lat=lats[j],lon=lons[j],
                                   keep=c("nauticalDawn","sunrise","sunset","nauticalDusk"),
                                   tz="UTC")
        nextNauticalDawn = c(dayData$nauticalDawn[2:length(dayData$date)],as.POSIXct(dateEnd))
        dayData$nextNauticalDawn = nextNauticalDawn
        
        # find which day phase each presence bin falls into 
        # REDO THIS WITH HISTC
        dayPhase = rep(NaN,length(pres))
        for (b in 1:length(pres)){
          if (numel(which(reducedBinVec[pres[b]]>=dayData$nauticalDawn & reducedBinVec[pres[b]]<dayData$sunrise))!=0){
            dayPhase[b] = "Dawn"
          } else if (numel(which(reducedBinVec[pres[b]]>=dayData$sunrise & reducedBinVec[pres[b]]<dayData$sunset))!=0) {
            dayPhase[b] = "Day"
          } else if (numel(which(reducedBinVec[pres[b]]>=dayData$sunset & reducedBinVec[pres[b]]<dayData$nauticalDusk))!=0) {
            dayPhase[b] = "Dusk"
          } else if (numel(which(reducedBinVec[pres[b]]>=dayData$nauticalDusk & reducedBinVec[pres[b]]<dayData$nextNauticalDawn))!=0) {
            dayPhase[b] = "Night"
          } else if (numel(which(reducedBinVec[pres[b]]>=dateStart & reducedBinVec[pres[b]]<dayData$nauticalDawn[1]))!=0) {
            dayPhase[b] = "Night"
          }
        }
        
        # TO DO: sum presence in each day phase & normalize to phase length (normalized # bins per phase, instead of # clicks per bin per phase)
        
        # # find and remove presence which isn't in any day phase (shouldn't be any)
        # noGood = which(dayPhase=="NaN")
        # if (length(noGood)!=0){
        #   keep = pres[-noGood]
        # } else { keep = pres }
        
        DielDF = data.frame(thisSite[pres],as.factor(dayPhase))
        colnames(DielDF) = c("Presence","DayPhase")
        
        # Plot # clicks per bin each day phase & test for differences
        if (length(levels(DielDF$DayPhase))>1){ 
          # Test significance of differences btwn day phases
          KW = kruskal.test(DielDF$Presence,DielDF$DayPhase,na.action=na.pass)
          pairCompCI = conover.test(DielDF$Presence,DielDF$DayPhase,method='bonferroni',label=TRUE,
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
            summarize(q25 = quantile(Presence,probs=0.25),
                      q50 = quantile(Presence,probs=0.50),
                      q75 = quantile(Presence,probs=0.75),
                      q90 = quantile(Presence,probs=0.9))
          iqr = quants$q75-quants$q25
          
          # Find and remove outliers for plotting (only if doing violin plots)
          if (type=="V"){
            out = DielDF %>% 
              group_by(DayPhase) %>%
              identify_outliers("Presence")
            if (numel(out!=0)){
              out = out[,-c(3,4)]
              outInd = which(!is.na(match_row(out,DielDF)))
              DielDF = DielDF[-outInd,]
            }
          }
          
          my_comparisons <- list( c("Dawn", "Day"), c("Dawn", "Dusk"), 
                                  c("Dawn", "Night"),c("Day","Dusk"), 
                                  c("Day","Night"),c("Dusk","Night"))
          
          # Plot presence per day phase
          if(type=="B"){ # boxplot for if there are clicks in more than one day phase
            ggplot(DielDF,aes(DayPhase,Presence)
            ) + geom_boxplot(varwidth=TRUE,
                             outlier.shape=NA
            ) + coord_cartesian(ylim = (1.75*max(iqr))+max(quants$q75)
            ) + scale_x_discrete(limits=c("Dawn","Day","Dusk","Night") # make sure NaNs aren't displayed, if they remain in the data
            ) + labs(title = paste(CTname, 'at',sites[j]),
                     x=(""),y="Clicks per 5-min Bin"
            ) + geom_signif(comparisons=my_comparisons,
                            y_position=(1.55*max(iqr))+max(quants$q75),
                            test=wilcox.test,
                            step_increase=0.07,
                            map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05))
            
            saveName = paste(dielDir,'/',CTname,'/minClicks50/',sites[j],"_4PhaseBoxplot.png",sep="")
          } else if (type=="V") { # violin plot for if there are clicks in more than one day phase
            ggplot(DielDF,aes(DayPhase,Presence)
            ) + geom_violin(scale="count",
                            draw_quantiles=c(0.25, 0.5, 0.75)
            ) + scale_x_discrete(limits=c("Dawn","Day","Dusk","Night") # make sure NaNs aren't displayed, if they remain in the data
            ) + labs(title = paste(CTname, 'at',sites[j]),
                     x=(""),y="Clicks per 5-min Bin"
            ) + geom_signif(comparisons=list(c("Day","Night")),
                            y_position=(1.55*(max(iqr))+max(quants$q75)),
                            test=wilcox.test,
                            step_increase=0.07,
                            map_signif_level=c("***"=0.001, "**"=0.01, "*"=0.05))
            
            saveName = paste(dielDir,'/',CTname,'/minClicks50/',sites[j],"_4PhaseViolinplot.png",sep="")
          }
          
          ggsave(saveName,device="png")
          
        } else if (length(levels(DielDF$DayPhase))==1) {
          if (type=="B"){ # boxplot for clicks in only one day phase
            ggplot(DielDF,aes(DayPhase,Presence)
            ) + geom_boxplot(varwidth=TRUE
            ) + annotate("text",
                         x=levels(DielDF$DayPhase),y=1.05*max(DielDF$Presence),
                         label=paste("Presence only during ",levels(DielDF$DayPhase),sep="")
            ) + labs(title = paste(CTname, 'at',sites[j]),
                     x=(""),y="Click per 5-min Bin")
            
            saveName = paste(dielDir,'/',CTname,'/',sites[j],"_DielBoxplot.png",sep="")
          } else if (type=="V") { # violin plot for clicks in only one day phase
            ggplot(DielDF,aes(DayPhase,Presence)
            ) + geom_violin(scale="count",
                            draw_quantiles=c(0.25, 0.5, 0.25)
            ) + annotate("text",
                         x=levels(DielDF$DayPhase),y=1.05*max(DielDF$Presence),
                         label=paste("Presence only during ",levels(DielDF$DayPhase),sep="")
            ) + labs(title = paste(CTname, 'at',sites[j]),
                     x=(""),y="Click per 5-min Bin")
            saveName = paste(dielDir,'/',CTname,'/',sites[j],"_DielViolinplot.png",sep="")
          }

          ggsave(saveName,device="png")
        }
        
      } else if (phases==2){   #  DAY/NIGHT ANALYSIS
        
        # Load normalized bin times
        whichNormTimeFile = str_detect(normTimeFileList,sites[j]) 
        normTimes = data.frame(read.csv(normTimeFileList[whichNormTimeFile]))  # load file
        normBinTimes = normTimes$Bin[-noDat]
        
        # create vector coding for day/night
        dayNight = rep(0,length(normBinTimes))
        dayNight[normBinTimes>0] = "Night"
        dayNight[normBinTimes<0] = "Day"
        # find and remove presence which falls right on sunrise/sunset
        neither = which(dayNight[pres]==0)
        if (length(neither)!=0){
          keep = pres[-neither]
        } else { keep = pres }
        
        DielDF = data.frame(thisSite[keep],as.factor(dayNight[keep]))
        colnames(DielDF) = c("Presence","DayNight")
        
        # calculate quantiles for plotting limits
        quants = DielDF %>%
          group_by(DayPhase) %>%
          summarize(q25 = quantile(Presence,probs=0.25),
                    q50 = quantile(Presence,probs=0.50),
                    q75 = quantile(Presence,probs=0.75),
                    q90 = quantile(Presence,probs=0.9))
        iqr = quants&q75 - quants$q25
        
        # Find and remove outliers for plotting (only for violin plots)
        if (type=="V") {
          out = DielDF %>% 
            group_by(DayPhase) %>%
            identify_outliers("Presence")
          out = out[,-c(3,4)]
          outInd = which(!is.na(match_row(out,DielDF)))
          DielDF = DielDF[-outInd,]
        }
        
        # Plot # clicks per bin in day/night & test for differences
        if (length(levels(DielDF$DayNight))==2){
          # Test significance of differences btwn day & night clicking
          KW = kruskal.test(DielDF$Presence,DielDF$DayNight,na.action=na.pass)
          
          # Plot day vs night presence
          if (type=="B"){ # boxplot
            ggplot(DielDF,aes(DayNight,Presence)
            ) + geom_boxplot(varwidth=TRUE,
                             outlier.shape=NA
            ) + coord_cartesian(ylim = (1.75*max(iqr))+max(quants$q75)
            ) + annotate("text",
                         x=levels(DielDF$DayNight)[1],y=1.2*max(quants$q75),
                         label=paste("P Value: ",KW$p.value,sep="")
            ) + labs(title = paste(CTname, 'at',sites[j]),
                     x=(""),y="Clicks per 5-min Bin")
            saveName = paste(dielDir,'/',CTname,'/minClicks50/',sites[j],"_DielBoxplot.png",sep="")
          } else if (type=="V") { # violin plot
            ggplot(DielDF,aes(DayNight,Presence)
            ) + geom_violin(scale="count",
                            draw_quantiles=c(0.25, 0.5, 0.25)
            ) + annotate("text",
                         x=levels(DielDF$DayNight)[1],y=1.2*max(quants$q75),
                         label=paste("P Value: ",KW$p.value,sep="")
            ) + labs(title = paste(CTname, 'at',sites[j]),
                     x=(""),y="Clicks per 5-min Bin")
            saveName = paste(dielDir,'/',CTname,'/minClicks50/',sites[j],"_DielViolinplot.png",sep="")
          }
          
          ggsave(saveName,device="png")
          
        } else {
          if (type=="B"){ # boxplot
            ggplot(DielDF,aes(DayNight,Presence)
            ) + geom_boxplot(varwidth=TRUE

            ) + annotate("text",
                         x=levels(DielDF$DayNight),y=1.05*max(DielDF$Presence),
                         label=paste("Presence only during ",levels(DielDF$DayNight),sep="")
            ) + labs(title = paste(CTname, 'at',sites[j]),
                     x=(""),y="Click per 5-min Bin")
            saveName = paste(dielDir,'/',CTname,'/minClicks50/',sites[j],"_DielBoxplot.png",sep="")
          } else if (type=="V") { # violin plot
            ggplot(DielDF,aes(DayNight,Presence)
            ) + geom_violin(scale="count",
                            draw_quantiles=c(0.25, 0.5, 0.25)
            ) + annotate("text",
                         x=levels(DielDF$DayNight),y=1.05*max(DielDF$Presence),
                         label=paste("Presence only during ",levels(DielDF$DayNight),sep="")
            ) + labs(title = paste(CTname, 'at',sites[j]),
                     x=(""),y="Click per 5-min Bin")
            saveName = paste(dielDir,'/',CTname,'/minClicks50/',sites[j],"_DielViolinplot.png",sep="")
          }
          
          ggsave(saveName,device="png")
        }
      }
    }
    
  }
}



