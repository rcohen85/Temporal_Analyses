# Plot histograms of temporal and lunar activity patterns
library(gridExtra)

DFDir = 'J:/Chpt_2/TimeSeries_ScaledByEffortError'
dielDir = 'J:/Chpt_2/Diel Plots'
lunDir = 'J:/Chpt_2/Lunar Plots'
int = "5minBin"

dfList = list.files(path=DFDir,pattern=paste('*',int,'_Master.csv',sep=""),
                    full.names=TRUE,recursive=FALSE,
                    include.dirs=FALSE,no..=TRUE)
lunList = list.files(path=DFDir,pattern=paste('*',int,'_MasterLun.csv',sep=""),
                     full.names=TRUE,recursive=FALSE,
                     include.dirs=FALSE,no..=TRUE)

for (i in 1:length(dfList)){
  thisSite = data.frame(read.csv(dfList[i]))
  CTname = str_remove(dfList[i],paste(DFDir,'/',sep="")) # get the species/CT name
  site = str_remove(CTname,paste("_",int,"_Master.csv",sep=""))
  site = sub(".*_","",site)
  CTname = sub("_.*","",CTname)
  if (str_detect(CTname,"Atl")){
    CTname = "Gervais"
  }
  
  # if it doesn't already exist, create directory to save figures
  if (!dir.exists(paste(dielDir,'/',CTname,sep=""))){
    dir.create(paste(dielDir,'/',CTname,sep=""))
  }
  
  thisSite$DayPhase = as.numeric(as.factor(thisSite$DayPhase))
  pres = which(thisSite$Presence==1)
  presDF = thisSite[pres,]
  
  NTedges = seq(-1,1,by=0.1)
  NTDF = data.frame(Obs=histc(thisSite$NormTime,edges=NTedges)$cnt[1:20]/dim(thisSite)[1],
                    Pres=histc(presDF$NormTime,edges=NTedges)$cnt[1:20]/dim(presDF)[1],
                    Edges=seq(-.95,0.95,by=0.1))
  
  NT = ggplot(NTDF
  )+geom_col(aes(x=Edges,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=Edges,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(-1,0,1),
                           label=c("Sunrise","Sunset","Sunrise")
  ) + coord_cartesian(xlim=c(-1,1)
  ) + labs(y="Normalized Presence Counts",x=NULL,title="Normalized Time of Day"
  ) + theme_minimal()
  
  DPedges = seq(0.5,4.5,by=1)
  DPDF = data.frame(Obs=histc(thisSite$DayPhase,edges=DPedges)$cnt[1:4]/dim(thisSite)[1],
                    Pres=histc(presDF$DayPhase,edges=DPedges)$cnt[1:4]/dim(presDF)[1],
                    Edges=seq(1,4,by=1))
  
  DP = ggplot(DPDF
  )+geom_col(aes(x=Edges,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=Edges,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(1,2,3,4),
                         label=c("Dawn","Day","Dusk","Night")
  ) + coord_cartesian(xlim=c(0.5,4.5)
  ) + labs(y="Normalized Presence Counts",x=NULL,title="Day Phase"
  ) + theme_minimal()
  
  JDedges = seq(1,365,by=14)
  JDDF = data.frame(Obs=histc(thisSite$JulianDay,edges=JDedges)$cnt[1:26]/dim(thisSite)[1],
                    Pres=histc(presDF$JulianDay,edges=JDedges)$cnt[1:26]/dim(presDF)[1],
                    Edges=seq(4.5,361.5,by=14))
  
  JD = ggplot(JDDF
  )+geom_col(aes(x=Edges,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=Edges,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),
                         label=c("J","F","M","A","M","J","J","A","S","O","N","D")
  ) + coord_cartesian(xlim=c(1,365)
  ) + labs(y="Normalized Presence Counts",x=NULL,title="Julian Day"
  ) + theme_minimal()
  
  Yredges = seq(0.5,3.5,by=1)
  YrDF = data.frame(Obs=histc(thisSite$StudyYear,edges=Yredges)$cnt[1:3]/dim(thisSite)[1],
                    Pres=histc(presDF$StudyYear,edges=Yredges)$cnt[1:3]/dim(presDF)[1],
                    Edges=seq(1,3,by=1))
  
  Yr = ggplot(YrDF
  )+geom_col(aes(x=Edges,y=Obs),
             fill='#bfbfbf',
             alpha=0.6
  )+geom_col(aes(x=Edges,y=Pres),
             fill='#66B2FF',
             alpha=0.5
  ) + scale_x_continuous(breaks=c(1,2,3),
  ) + coord_cartesian(xlim=c(0.5,3.5)
  ) + labs(y="Normalized Presence Counts",x=NULL,title="Study Year"
  ) + theme_minimal()
  
  
  png(file=paste(dielDir,'/',CTname,'/',site,"_TempHist.png",sep=""),width = 800, height = 800, units = "px")
  grid.arrange(NT,DP,JD,Yr,ncol=2,nrow=2,top=paste(CTname,'at',site))
  while (dev.cur()>1) {dev.off()}
  
}

for (i in 1:length(lunList)){
  
  masterLun = data.frame(read.csv(lunList[i]))
  CTname = str_remove(lunList[i],paste(DFDir,'/',sep="")) # get the species/CT name
  site = str_remove(CTname,paste("_",int,"_MasterLun.csv",sep=""))
  site = sub(".*_","",site)
  CTname = sub("_.*","",CTname)
  if (str_detect(CTname,"Atl")){
    CTname = "Gervais"
  }
  
  # if it doesn't already exist, create directory to save figures
  if (!dir.exists(paste(lunDir,'/',CTname,sep=""))){
    dir.create(paste(lunDir,'/',CTname,sep=""))
  }
  
  masterLun$MoonAltitude = masterLun$MoonAltitude*(180/pi) # convert altitude from radians to degrees
  pres = which(masterLun$NightPres>0)
  nightPresDF = masterLun[pres,]
  
  
  Ph = ggplot(nightPresDF,aes(x=MoonPhase)
  )+geom_histogram(aes(y=..count../dim(nightPresDF)[1]),
                   fill='#66B2FF',
                   binwidth = 0.1,
                   alpha=0.5
  # ) + scale_x_continuous(breaks=c(0.01,0.5,0.99),
  #                        label=c("New Moon","Full Moon","New Moon")
  ) + labs(y="Normalized Presence Counts",x=NULL,title="Moon Phase"
  ) + theme_minimal()
  
  Phdat = ggplot(masterLun,aes(x=MoonPhase)
  )+geom_histogram(aes(y=..count../dim(masterLun)[1]),
                   fill='#66B2FF',
                   binwidth = 0.1,
                   alpha=0.5
  ) + scale_x_continuous(breaks=c(0.01,0.5,0.99),
                         label=c("New Moon","Full Moon","New Moon")
  ) + labs(y="Normalized Data Counts",x=NULL
  ) + theme_minimal()
  
  Alt = ggplot(nightPresDF,aes(x=MoonAltitude)
  )+geom_histogram(aes(y=..count../dim(nightPresDF)[1]),
                   fill='#66B2FF',
                   binwidth = 10,
                   alpha=0.5
  ) + scale_x_continuous(breaks=c(-45,0,45),
                         label=c(paste("-45",intToUtf8(176),sep=""),
                                 paste("0",intToUtf8(176),sep=""),
                                 paste("45",intToUtf8(176),sep=""))
  ) + labs(y=NULL,x=NULL,title="Moon Altitude"
  ) + theme_minimal()
  
  Altdat = ggplot(masterLun,aes(x=MoonAltitude)
  )+geom_histogram(aes(y=..count../dim(masterLun)[1]),
                   fill='#66B2FF',
                   binwidth = 10,
                   alpha=0.5
  ) + scale_x_continuous(breaks=c(-45,0,45),
                         label=c(paste("-45",intToUtf8(176),sep=""),
                                 paste("0",intToUtf8(176),sep=""),
                                 paste("45",intToUtf8(176),sep=""))
  ) + labs(y=NULL,x=NULL
  ) + theme_minimal()
  
  Pres = ggplot(nightPresDF,aes(x=as.numeric(as.factor(MoonPres)))
  )+geom_histogram(aes(y=..count../dim(nightPresDF)[1]),
                   fill='#66B2FF',
                   binwidth = 0.5,
                   alpha=0.5
  ) + scale_x_continuous(breaks=c(1,2,3),
                       label=c("Before","MoonUp","After")
  ) + labs(y=NULL,x=NULL,title="Moon Presence"
  ) + theme_minimal()
  
  Presdat = ggplot(masterLun,aes(x=as.numeric(as.factor(MoonPres)))
  )+geom_histogram(aes(y=..count../dim(masterLun)[1]),
                   fill='#66B2FF',
                   binwidth = 0.5,
                   alpha=0.5
  ) + scale_x_continuous(breaks=c(1,2,3),
                       label=c("Before","MoonUp","After")
  ) + labs(y=NULL,x=NULL
  ) + theme_minimal()
  
  png(file=paste(lunDir,'/',CTname,'/',site,"_LunHists.png",sep=""),width = 1200, height = 600, units = "px")
  grid.arrange(Ph,Alt,Pres,Phdat,Altdat,Presdat,ncol=3,nrow=2,top=paste(CTname,'at',site))
  while (dev.cur()>1) {dev.off()}
  
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


