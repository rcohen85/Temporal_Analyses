# Plot histograms of temporal and lunar activity patterns
library(gridExtra)
library(stringr)
library(pracma)
library(ggplot2)
library(lubridate)

DFDir = 'J:/Chpt_2/TimeSeries_ScaledByEffortError'
tempDir = 'J:/Chpt_2/DataDists'
# lunDir = 'J:/Chpt_2/Lunar Plots'
int = "5minBin"

dfList = list.files(path=DFDir,pattern=paste('*',int,'_MasterTempLun.csv',sep=""),
                    full.names=TRUE,recursive=FALSE,
                    include.dirs=FALSE,no..=TRUE)
lunList = list.files(path=DFDir,pattern=paste('*',int,'_MasterLun.csv',sep=""),
                     full.names=TRUE,recursive=FALSE,
                     include.dirs=FALSE,no..=TRUE)

for (i in 1:length(dfList)){
  thisSite = data.frame(read.csv(dfList[i]))
  CTname = str_remove(dfList[i],paste(DFDir,'/',sep="")) # get the species/CT name
  site = str_remove(CTname,paste("_",int,"_MasterTempLun.csv",sep=""))
  site = sub(".*_","",site)
  CTname = sub("_.*","",CTname)
  if (str_detect(CTname,"Atl")){
    CTname = "Gervais"
  }
  thisSite$TimeStamp = as_datetime(thisSite$TimeStamp)

  # if it doesn't already exist, create directory to save figures
  if (!dir.exists(paste(tempDir,'/',CTname,sep=""))){
    dir.create(paste(tempDir,'/',CTname,sep=""))
  }

  # thisSite$DayPhase = as.numeric(as.factor(thisSite$DayPhase))
  pres = which(thisSite$Presence==1)
  presDF = thisSite[pres,]

  # Plot presence bins across all data for each covar
  NTedges = seq(-1,1,by=0.1)
  NTDF = data.frame(Obs=histc(thisSite$NormTime,edges=NTedges)$cnt[1:20]/dim(thisSite)[1],
                    Pres=histc(presDF$NormTime,edges=NTedges)$cnt[1:20]/dim(presDF)[1],
                    BinCenter=seq(-.95,0.95,by=0.1))

  NT = ggplot(NTDF
  )+geom_col(aes(x=BinCenter,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=BinCenter,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(-1,0,1),
                           label=c("Sunrise","Sunset","Sunrise")
  ) + coord_cartesian(xlim=c(-1,1)
  ) + labs(y="Normalized Counts",x=NULL,title="Normalized Time of Day"
  ) + theme_minimal()

  # DPedges = seq(0.5,4.5,by=1)
  # DPDF = data.frame(Obs=histc(thisSite$DayPhase,edges=DPedges)$cnt[1:4]/dim(thisSite)[1],
  #                   Pres=histc(presDF$DayPhase,edges=DPedges)$cnt[1:4]/dim(presDF)[1],
  #                   BinCenter=seq(1,4,by=1))
  # 
  # DP = ggplot(DPDF
  # )+geom_col(aes(x=BinCenter,y=Obs),
  #            fill='#bfbfbf',
  #            alpha=0.6
  # )+geom_col(aes(x=BinCenter,y=Pres),
  #            fill='#66B2FF',
  #            alpha=0.5
  # ) + scale_x_continuous(breaks=c(1,2,3,4),
  #                        label=c("Dawn","Day","Dusk","Night")
  # ) + coord_cartesian(xlim=c(0.5,4.5)
  # ) + labs(y="Normalized Counts",x=NULL,title="Day Phase"
  # ) + theme_minimal()
  
  PHedges = seq(0,1,by=0.05)
  MPHDF = data.frame(Obs=histc(thisSite$MoonPhase,edges=PHedges)$cnt[1:20]/dim(thisSite)[1],
                    Pres=histc(presDF$MoonPhase,edges=PHedges)$cnt[1:20]/dim(presDF)[1],
                    BinCenter=seq(0.025,0.975,by=0.05))

  MPh = ggplot(MPHDF
  )+geom_col(aes(x=BinCenter,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=BinCenter,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(0,0.5,1),
                         label=c("New Moon","Full Moon","New Moon")
  ) + coord_cartesian(xlim=c(0,1)
  ) + labs(y="Normalized Counts",x=NULL,title="Moon Phase"
  ) + theme_minimal()

  JDedges = seq(1,365,by=14)
  JDDF = data.frame(Obs=histc(thisSite$JulianDay,edges=JDedges)$cnt[1:26]/dim(thisSite)[1],
                    Pres=histc(presDF$JulianDay,edges=JDedges)$cnt[1:26]/dim(presDF)[1],
                    BinCenter=seq(4.5,361.5,by=14))

  JD = ggplot(JDDF
  )+geom_col(aes(x=BinCenter,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=BinCenter,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),
                         label=c("J","F","M","A","M","J","J","A","S","O","N","D")
  ) + coord_cartesian(xlim=c(1,365)
  ) + labs(y="Normalized Counts",x=NULL,title="Julian Day"
  ) + theme_minimal()

  Yredges = seq(0.5,3.5,by=1)
  YrDF = data.frame(Obs=histc(thisSite$StudyYear,edges=Yredges)$cnt[1:3]/dim(thisSite)[1],
                    Pres=histc(presDF$StudyYear,edges=Yredges)$cnt[1:3]/dim(presDF)[1],
                    BinCenter=seq(1,3,by=1))

  Yr = ggplot(YrDF
  )+geom_col(aes(x=BinCenter,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=BinCenter,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(1,2,3),
  ) + coord_cartesian(xlim=c(0.5,3.5)
  ) + labs(y="Normalized Counts",x=NULL,title="Study Year"
  ) + theme_minimal()


  png(file=paste(tempDir,'/',CTname,'/',site,"_TempHist.png",sep=""),width = 1000, height = 250, units = "px")
  grid.arrange(NT,MPh,JD,Yr,ncol=4,nrow=1,top=paste(CTname,'at',site))
  while (dev.cur()>1) {dev.off()}
  
  
  # Plot NT presence bins in different seasons --------------------------------
  NTedges = seq(-1,1,by=0.1)
  winInd1 = which(month(thisSite$TimeStamp)==1)
  winInd2 = which(month(thisSite$TimeStamp[pres])==1)
  NTWinter = data.frame(Obs=histc(thisSite$NormTime[winInd1],edges=NTedges)$cnt[1:20]/dim(thisSite)[1],
                    Pres=histc(presDF$NormTime[winInd2],edges=NTedges)$cnt[1:20]/dim(presDF)[1],
                    BinCenter=seq(-.95,0.95,by=0.1))
  

  
  sprInd1 = which(month(thisSite$TimeStamp)==4)
  sprInd2 = which(month(thisSite$TimeStamp[pres])==4)
  NTSpring = data.frame(Obs=histc(thisSite$NormTime[sprInd1],edges=NTedges)$cnt[1:20]/dim(thisSite)[1],
                        Pres=histc(presDF$NormTime[sprInd2],edges=NTedges)$cnt[1:20]/dim(presDF)[1],
                        BinCenter=seq(-.95,0.95,by=0.1))
  
  sumInd1 = which(month(thisSite$TimeStamp)==7)
  sumInd2 = which(month(thisSite$TimeStamp[pres])==7)
  NTSummer = data.frame(Obs=histc(thisSite$NormTime[sumInd1],edges=NTedges)$cnt[1:20]/dim(thisSite)[1],
                        Pres=histc(presDF$NormTime[sumInd2],edges=NTedges)$cnt[1:20]/dim(presDF)[1],
                        BinCenter=seq(-.95,0.95,by=0.1))
  
  fallInd1 = which(month(thisSite$TimeStamp)==10)
  fallInd2 = which(month(thisSite$TimeStamp[pres])==10)
  NTFall = data.frame(Obs=histc(thisSite$NormTime[fallInd1],edges=NTedges)$cnt[1:20]/dim(thisSite)[1],
                        Pres=histc(presDF$NormTime[fallInd2],edges=NTedges)$cnt[1:20]/dim(presDF)[1],
                        BinCenter=seq(-.95,0.95,by=0.1))
  
  ymax = 1.05*max(c(max(NTWinter[,1:2]),max(NTSpring[,1:2]),max(NTSummer[,1:2]),max(NTFall[,1:2])))
  yrange = c(0,ymax)
  
  Win = ggplot(NTWinter
  )+geom_col(aes(x=BinCenter,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=BinCenter,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(-1,0,1),
                         label=c("Sunrise","Sunset","Sunrise")
  ) + coord_cartesian(xlim=c(-1,1),ylim=yrange
  ) + labs(y="Normalized Counts",x=NULL,title="January"
  ) + theme_minimal()
  
  Spr = ggplot(NTSpring
  )+geom_col(aes(x=BinCenter,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=BinCenter,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(-1,0,1),
                         label=c("Sunrise","Sunset","Sunrise")
  ) + coord_cartesian(xlim=c(-1,1),ylim=yrange
  ) + labs(y="Normalized Counts",x=NULL,title="April"
  ) + theme_minimal()
  
  Sum = ggplot(NTSummer
  )+geom_col(aes(x=BinCenter,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=BinCenter,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(-1,0,1),
                         label=c("Sunrise","Sunset","Sunrise")
  ) + coord_cartesian(xlim=c(-1,1),ylim=yrange
  ) + labs(y="Normalized Counts",x=NULL,title="July"
  ) + theme_minimal()
  
  Fall = ggplot(NTFall
  )+geom_col(aes(x=BinCenter,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=BinCenter,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(-1,0,1),
                         label=c("Sunrise","Sunset","Sunrise")
  ) + coord_cartesian(xlim=c(-1,1),ylim=yrange
  ) + labs(y="Normalized Counts",x=NULL,title="October"
  ) + theme_minimal()
  
  png(file=paste(tempDir,'/',CTname,'/',site,"_NTJDHist.png",sep=""),width = 1000, height = 250, units = "px")
  grid.arrange(Win,Spr,Sum,Fall,ncol=4,nrow=1,top=paste(CTname,'at',site,'NormTime'))
  while (dev.cur()>1) {dev.off()}
  
  
  ## Plot NT presence bins at new, 1st quarter, full, & 3rd quarter ------------------------
  NTedges = seq(-1,1,by=0.1)
  newInd1 = c(which(thisSite$MoonPhase<0.01), which(thisSite$MoonPhase>0.99))
  newInd2 = c(which(thisSite$MoonPhase[pres]<0.01), which(thisSite$MoonPhase[pres]>0.99))
  NTNew = data.frame(Obs=histc(thisSite$NormTime[newInd1],edges=NTedges)$cnt[1:20]/dim(thisSite)[1],
                        Pres=histc(presDF$NormTime[newInd2],edges=NTedges)$cnt[1:20]/dim(presDF)[1],
                        BinCenter=seq(-.95,0.95,by=0.1))
  
  
  
  firstInd1 = which(thisSite$MoonPhase<0.26 & thisSite$MoonPhase>0.24)
  firstInd2 = which(thisSite$MoonPhase[pres]<0.26 & thisSite$MoonPhase[pres]>0.24)
  NTFirst = data.frame(Obs=histc(thisSite$NormTime[firstInd1],edges=NTedges)$cnt[1:20]/dim(thisSite)[1],
                        Pres=histc(presDF$NormTime[firstInd2],edges=NTedges)$cnt[1:20]/dim(presDF)[1],
                        BinCenter=seq(-.95,0.95,by=0.1))
  
  fullInd1 = which(thisSite$MoonPhase<0.51 & thisSite$MoonPhase>0.49)
  fullInd2 = which(thisSite$MoonPhase[pres]<0.51 & thisSite$MoonPhase[pres]>0.49)
  NTFull = data.frame(Obs=histc(thisSite$NormTime[fullInd1],edges=NTedges)$cnt[1:20]/dim(thisSite)[1],
                        Pres=histc(presDF$NormTime[fullInd2],edges=NTedges)$cnt[1:20]/dim(presDF)[1],
                        BinCenter=seq(-.95,0.95,by=0.1))
  
  thirdInd1 = which(thisSite$MoonPhase<0.76 & thisSite$MoonPhase>0.74)
  thirdInd2 = which(thisSite$MoonPhase[pres]<0.76 & thisSite$MoonPhase[pres]>0.74)
  NTThird = data.frame(Obs=histc(thisSite$NormTime[thirdInd1],edges=NTedges)$cnt[1:20]/dim(thisSite)[1],
                      Pres=histc(presDF$NormTime[thirdInd2],edges=NTedges)$cnt[1:20]/dim(presDF)[1],
                      BinCenter=seq(-.95,0.95,by=0.1))
  
  ymax = 1.05*max(c(max(NTNew[,1:2]),max(NTFirst[,1:2]),max(NTFull[,1:2]),max(NTThird[,1:2])))
  yrange = c(0,ymax)
  
  New = ggplot(NTNew
  )+geom_col(aes(x=BinCenter,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=BinCenter,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(-1,0,1),
                         label=c("Sunrise","Sunset","Sunrise")
  ) + coord_cartesian(xlim=c(-1,1),ylim=yrange
  ) + labs(y="Normalized Counts",x=NULL,title="New Moon"
  ) + theme_minimal()
  
  First = ggplot(NTFirst
  )+geom_col(aes(x=BinCenter,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=BinCenter,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(-1,0,1),
                         label=c("Sunrise","Sunset","Sunrise")
  ) + coord_cartesian(xlim=c(-1,1),ylim=yrange
  ) + labs(y="Normalized Counts",x=NULL,title="First Quarter"
  ) + theme_minimal()
  
  Full = ggplot(NTFull
  )+geom_col(aes(x=BinCenter,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=BinCenter,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(-1,0,1),
                         label=c("Sunrise","Sunset","Sunrise")
  ) + coord_cartesian(xlim=c(-1,1),ylim=yrange
  ) + labs(y="Normalized Counts",x=NULL,title="Full Moon"
  ) + theme_minimal()
  
  Third = ggplot(NTThird
  )+geom_col(aes(x=BinCenter,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=BinCenter,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(-1,0,1),
                         label=c("Sunrise","Sunset","Sunrise")
  ) + coord_cartesian(xlim=c(-1,1),ylim=yrange
  ) + labs(y="Normalized Counts",x=NULL,title="Third Quarter"
  ) + theme_minimal()
  
  png(file=paste(tempDir,'/',CTname,'/',site,"_NTMPhHist.png",sep=""),width = 1000, height = 250, units = "px")
  grid.arrange(New,First,Full,Third,ncol=4,nrow=1,top=paste(CTname,'at',site,'NormTime'))
  while (dev.cur()>1) {dev.off()}
  
  ## Plot MPh presence bins in day vs night ----------------------
  PHedges = seq(0,1,by=0.05)
  dayInd1 = which(thisSite$NormTime<0)
  dayInd2 = which(thisSite$NormTime[pres]<0)
  MPhDay = data.frame(Obs=histc(thisSite$MoonPhase[dayInd1],edges=PHedges)$cnt[1:20]/dim(thisSite)[1],
                     Pres=histc(presDF$MoonPhase[dayInd2],edges=PHedges)$cnt[1:20]/dim(presDF)[1],
                     BinCenter=seq(0.025,0.975,by=0.05))
  
  nightInd1 = which(thisSite$NormTime>0)
  nightInd2 = which(thisSite$NormTime[pres]>0)
  MPhNight = data.frame(Obs=histc(thisSite$MoonPhase[nightInd1],edges=PHedges)$cnt[1:20]/dim(thisSite)[1],
                      Pres=histc(presDF$MoonPhase[nightInd2],edges=PHedges)$cnt[1:20]/dim(presDF)[1],
                      BinCenter=seq(0.025,0.975,by=0.05))
  
  ymax = 1.05*max(c(max(MPhDay[,1:2]),max(MPhNight[,1:2])))
  yrange = c(0,ymax)
  
  MPhDay = ggplot(MPhDay
  )+geom_col(aes(x=BinCenter,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=BinCenter,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(0,0.5,1),
                         label=c("New Moon","Full Moon","New Moon")
  ) + coord_cartesian(xlim=c(0,1),ylim=yrange
  ) + labs(y="Normalized Counts",x=NULL,title="Day"
  ) + theme_minimal()
  
  MPhNight = ggplot(MPhNight
  )+geom_col(aes(x=BinCenter,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=BinCenter,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(0,0.5,1),
                         label=c("New Moon","Full Moon","New Moon")
  ) + coord_cartesian(xlim=c(0,1),ylim=yrange
  ) + labs(y="Normalized Counts",x=NULL,title="Night"
  ) + theme_minimal()
  
  
  png(file=paste(tempDir,'/',CTname,'/',site,"_MPhNTHist.png",sep=""),width = 500, height = 250, units = "px")
  grid.arrange(MPhDay,MPhNight,ncol=2,nrow=1,top=paste(CTname,'at',site,'MoonPhase'))
  while (dev.cur()>1) {dev.off()}
  
}

for (i in 1:length(lunList)){
  # 
  # masterLun = data.frame(read.csv(lunList[i]))
  # CTname = str_remove(lunList[i],paste(DFDir,'/',sep="")) # get the species/CT name
  # site = str_remove(CTname,paste("_",int,"_MasterLun.csv",sep=""))
  # site = sub(".*_","",site)
  # CTname = sub("_.*","",CTname)
  # if (str_detect(CTname,"Atl")){
  #   CTname = "Gervais"
  # }
  # 
  # # if it doesn't already exist, create directory to save figures
  # if (!dir.exists(paste(lunDir,'/',CTname,sep=""))){
  #   dir.create(paste(lunDir,'/',CTname,sep=""))
  # }
  # 
  # masterLun$MoonPres = as.numeric(as.factor(masterLun$MoonPres))
  # masterLun$MoonAltitude = masterLun$MoonAltitude*(180/pi) # convert altitude from radians to degrees
  # pres = which(masterLun$NightPres>0)
  # nightPresDF = masterLun[pres,]
  # 
  # PHedges = seq(0,1,by=0.05)
  # PHDF = data.frame(Obs=histc(masterLun$MoonPhase,edges=PHedges)$cnt[1:20]/dim(masterLun)[1],
  #                   Pres=histc(nightPresDF$MoonPhase,edges=PHedges)$cnt[1:20]/dim(nightPresDF)[1],
  #                   BinCenter=seq(0.025,0.975,by=0.05))
  # 
  # Ph = ggplot(PHDF
  # )+geom_col(aes(x=BinCenter,y=Obs),
  #            fill='#bfbfbf',
  #            alpha=0.6
  # )+geom_col(aes(x=BinCenter,y=Pres),
  #            fill='#66B2FF',
  #            alpha=0.5
  # ) + scale_x_continuous(breaks=c(0,0.5,1),
  #                        label=c("New Moon","Full Moon","New Moon")
  # ) + coord_cartesian(xlim=c(0,1)
  # ) + labs(y="Normalized Counts",x=NULL,title="Moon Phase"
  # ) + theme_minimal()
  # 
  # Altedges = seq(-90,90,by=10)
  # AltDF = data.frame(Obs=histc(masterLun$MoonAltitude,edges=Altedges)$cnt[1:18]/dim(masterLun)[1],
  #                   Pres=histc(nightPresDF$MoonAltitude,edges=Altedges)$cnt[1:18]/dim(nightPresDF)[1],
  #                   BinCenter=seq(-85,85,by=10))
  # 
  # Alt = ggplot(AltDF
  # )+geom_col(aes(x=BinCenter,y=Obs),
  #            fill='#bfbfbf',
  #            alpha=0.6
  # )+geom_col(aes(x=BinCenter,y=Pres),
  #            fill='#66B2FF',
  #            alpha=0.5
  # ) + scale_x_continuous(breaks=c(-90,-45,0,45,90),
  #                        label=c(paste("-90",intToUtf8(176),sep=""),
  #                                paste("-45",intToUtf8(176),sep=""),
  #                                paste("0",intToUtf8(176),sep=""),
  #                                paste("45",intToUtf8(176),sep=""),
  #                                paste("90",intToUtf8(176),sep=""))
  # ) + coord_cartesian(xlim=c(-90,90)
  # ) + labs(y="Normalized Counts",x=NULL,title="Moon Altitude"
  # ) + theme_minimal()
  # 
  # Presedges = seq(0.5,3.5,by=1)
  # PresDF = data.frame(Obs=histc(masterLun$MoonPres,edges=Presedges)$cnt[1:3]/dim(masterLun)[1],
  #                   Pres=histc(nightPresDF$MoonPres,edges=Presedges)$cnt[1:3]/dim(nightPresDF)[1],
  #                   BinCenter=seq(1,3,by=1))
  # 
  # Pres = ggplot(PresDF
  # )+geom_col(aes(x=BinCenter,y=Obs),
  #            fill='#bfbfbf',
  #            alpha=0.6
  # )+geom_col(aes(x=BinCenter,y=Pres),
  #            fill='#66B2FF',
  #            alpha=0.5
  # ) + scale_x_continuous(breaks=c(1,2,3),
  #                        label=c("Before","MoonUp","After")
  # ) + coord_cartesian(xlim=c(0.5,3.5)
  # ) + labs(y="Normalized Counts",x=NULL,title="Moon Presence"
  # ) + theme_minimal()
  # 
  # png(file=paste(lunDir,'/',CTname,'/',site,"_LunHists.png",sep=""),width = 1000, height = 300, units = "px")
  # grid.arrange(Ph,Alt,Pres,ncol=3,nrow=1,top=paste(CTname,'at',site))
  # while (dev.cur()>1) {dev.off()}
  # 
}


## Plot presence and effort on separate histograms --------------
#
# NT = ggplot(presDF,aes(x=NormTime)
# )+geom_histogram(aes(y=..count../dim(presDF)[1]),
#                  fill='#66B2FF',
#                  binwidth = 0.1,
#                  alpha=0.5
# ) + scale_x_continuous(breaks=c(-1,0,1),
#                        label=c("Sunrise","Sunset","Sunrise")
# ) + coord_cartesian(xlim=c(-1,1)
# ) + labs(y="Normalized Presence Counts",x=NULL,title="Normalized Time of Day"
# ) + theme_minimal()

# NTdat = ggplot(thisSite,aes(x=NormTime)
# )+geom_histogram(aes(y=..count../dim(thisSite)[1]),
#                  fill='#bfbfbf',
#                  binwidth = 0.1,
#                  alpha=0.5
# ) + scale_x_continuous(breaks=c(-0.99,0,0.99),
#                        label=c("Sunrise","Sunset","Sunrise")
# ) + coord_cartesian(xlim=c(-1,1)
# ) + labs(y="Normalized Data Counts",x=NULL
# ) + theme_minimal()
#
# DP = ggplot(presDF,aes(x=DayPhase)
# )+geom_histogram(aes(y=..count../dim(presDF)[1]),
#                  fill='#66B2FF',
#                  binwidth = 0.5,
#                  alpha=0.5
# ) + scale_x_continuous(breaks=c(1,2,3,4),
#                      labels=c("Dawn","Day","Dusk","Night")
# ) + coord_cartesian(xlim=c(0,5)
# ) + labs(y=NULL,x=NULL,title="Day Phase"
# ) + theme_minimal()
# 
# DPdat = ggplot(thisSite,aes(x=DayPhase)
# )+geom_histogram(aes(y=..count../dim(thisSite)[1]),
#                  fill='#66B2FF',
#                  binwidth = 0.5,
#                  alpha=0.5
# ) + scale_x_continuous(breaks=c(1,2,3,4),
#                      label=c("Dawn","Day","Dusk","Night")
# ) + coord_cartesian(xlim=c(0,5)
# ) + labs(y=NULL,x=NULL
# ) + theme_minimal()
#
# JD = ggplot(presDF,aes(x=JulianDay)
# )+geom_histogram(aes(y=..count../dim(presDF)[1]),
#                  fill='#66B2FF',
#                  binwidth = 7,
#                  alpha=0.5
# ) + scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),
#                        label=c("J","F","M","A","M","J","J","A","S","O","N","D")
# ) + coord_cartesian(xlim=c(1,365)
# ) + labs(y=NULL,x=NULL,title="Julian Day"
# ) + theme_minimal()
#
# Yr = ggplot(presDF,aes(x=StudyYear)
# )+geom_histogram(aes(y=..count../dim(presDF)[1]),
#                  fill='#66B2FF',
#                  binwidth = 1,
#                  alpha=0.5
# ) + coord_cartesian(xlim=c(1,3)
# ) + labs(y=NULL,x=NULL,title="Study Year"
# ) + theme_minimal()
#
# Ph = ggplot(nightPresDF,aes(x=MoonPhase)
# )+geom_histogram(aes(y=..count../dim(nightPresDF)[1]),
#                  fill='#66B2FF',
#                  binwidth = 0.1,
#                  alpha=0.5
# # ) + scale_x_continuous(breaks=c(0,0.5,1),
# #                        label=c("New Moon","Full Moon","New Moon")
# ) + labs(y="Normalized Presence Counts",x=NULL,title="Moon Phase"
# ) + theme_minimal()
# 
# Phdat = ggplot(masterLun,aes(x=MoonPhase)
# )+geom_histogram(aes(y=..count../dim(masterLun)[1]),
#                  fill='#66B2FF',
#                  binwidth = 0.1,
#                  alpha=0.5
# ) + scale_x_continuous(breaks=c(0.01,0.5,0.99),
#                        label=c("New Moon","Full Moon","New Moon")
# ) + labs(y="Normalized Data Counts",x=NULL
# ) + theme_minimal()
#
# Alt = ggplot(nightPresDF,aes(x=MoonAltitude)
# )+geom_histogram(aes(y=..count../dim(nightPresDF)[1]),
#                  fill='#66B2FF',
#                  binwidth = 10,
#                  alpha=0.5
# ) + scale_x_continuous(breaks=c(-45,0,45),
#                        label=c(paste("-45",intToUtf8(176),sep=""),
#                                paste("0",intToUtf8(176),sep=""),
#                                paste("45",intToUtf8(176),sep=""))
# ) + labs(y=NULL,x=NULL,title="Moon Altitude"
# ) + theme_minimal()
# 
# Altdat = ggplot(masterLun,aes(x=MoonAltitude)
# )+geom_histogram(aes(y=..count../dim(masterLun)[1]),
#                  fill='#66B2FF',
#                  binwidth = 10,
#                  alpha=0.5
# ) + scale_x_continuous(breaks=c(-45,0,45),
#                        label=c(paste("-45",intToUtf8(176),sep=""),
#                                paste("0",intToUtf8(176),sep=""),
#                                paste("45",intToUtf8(176),sep=""))
# ) + labs(y=NULL,x=NULL
# ) + theme_minimal()
#
# Pres = ggplot(nightPresDF,aes(x=as.numeric(as.factor(MoonPres)))
# )+geom_histogram(aes(y=..count../dim(nightPresDF)[1]),
#                  fill='#66B2FF',
#                  binwidth = 0.5,
#                  alpha=0.5
# ) + scale_x_continuous(breaks=c(1,2,3),
#                      label=c("Before","MoonUp","After")
# ) + labs(y=NULL,x=NULL,title="Moon Presence"
# ) + theme_minimal()
# 
# Presdat = ggplot(masterLun,aes(x=as.numeric(as.factor(MoonPres)))
# )+geom_histogram(aes(y=..count../dim(masterLun)[1]),
#                  fill='#66B2FF',
#                  binwidth = 0.5,
#                  alpha=0.5
# ) + scale_x_continuous(breaks=c(1,2,3),
#                      label=c("Before","MoonUp","After")
# ) + labs(y=NULL,x=NULL
# ) + theme_minimal()

