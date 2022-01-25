library(tidyverse)
library(lubridate)
library(pracma)
library(installr)
library(suncalc)

tsDir = 'I:/TimeSeries_ScaledByEffortError'
seasDir = 'I:/Seasonal_Plots'
dielDir = 'I:/Diel Plots'
illumDir = 'I:/IlluminationFiles'
normTimeDir = 'I:/NormTimes'
int = "5minBin"
start = '2016-05-01 00:00:00'
end = '2019-04-30 23:59:59'
goodIdx = c(1,2,4,8,11,13:19)
goodSite = list(c(8:10),c(8:10),c(),c(1:7),c(),c(),c(),c(1,4:11),c(),c(),c(1:11),
                c(),c(1:5),c(1:11),c(1:6),
                c(1:11),c(1:11),c(1:11),c(1:6))


lats = as.numeric(c(41.06165,40.22999,39.83295,
                    39.19192,38.37337,37.16452,
                    35.30183,33.66992,32.10527,
                    30.58295,30.27818))
lons = as.numeric(c(-66.35155,-67.97798,-69.98194,
                    -72.22735,-73.36985,-74.46585,
                    -74.87895,-75.9977,-77.09067,
                    -77.39002,-80.22085))


# Cycle through time series and organize covariates for each species at each site --------

dateStart = as.POSIXct(start,format="%Y-%m-%d %H:%M:%S",tz="GMT");
dateEnd = as.POSIXct(end,format="%Y-%m-%d %H:%M:%S",tz="GMT");
tsFileList = list.files(path=tsDir,pattern=paste('*_',int,'.csv',sep=""),
                        full.names=TRUE,recursive=FALSE,
                        include.dirs=FALSE,no..=TRUE)
normTimeFileList = list.files(path=normTimeDir,pattern='*.csv',
                              full.names=TRUE,recursive=FALSE,
                              include.dirs=FALSE,no..=TRUE)

for (i in goodIdx){  # for each species' file
  
  thisCT = data.frame(read.csv(tsFileList[i]))  # load file
  attribs = attributes(thisCT)  # find names of sites (all columns but first)
  sites = attribs$names[-1]
  
  CTname = str_remove(tsFileList[i],paste(tsDir,'/',sep="")) # get the species/CT name
  CTname = str_remove(CTname,paste('_',int,'.csv',sep=""))
  if (str_detect(CTname,"Atl")){
    CTname = "Gervais"
  }
  
  #dateVec = as.POSIXlt(thisCT$Date,tz="GMT",format="%d-%b-%Y")
  if (int=="Daily") {
    dateVec = as.POSIXlt(thisCT$Date,format="%d-%b-%Y",tz="GMT")
  } else if (int=="Hourly") {
    dateVec = as.POSIXlt(thisCT$Hour,format="%d-%b-%Y %H:%M:%S",tz="GMT")
  } else if (int=="5minBin") {
    dateVec = as.POSIXct(thisCT$Bin,format="%d-%b-%Y %H:%M:%S",tz="GMT")
  }
  
  for (j in unlist(goodSite[i])){ # for each site
    
    numClicks = as.numeric(thisCT[,j+1]) # get presence data for this site
    #pres = which(numClicks!=0)
    if (any(i==c(11,16:19))){
      pres = which(numClicks>=50) # consider dolphins "present" if more than 50 clicks in bin
    } else if (any(i==c(1,2,4,8,13:15))){
      pres = which(numClicks>=20) # consider BWs/Kogia/Pm "present" if more than 20 clicks in bin
    }
    
    if (!numel(pres)==0){ # if there is any presence
      
      # # plot presence to determine if it's worth modeling
      # plot(dateVec,numClicks,main=paste(CTname,'at',sites[j],sep=" "),
      #      ylab="# Clicks Per Bin")
      # 
      # userVote = c()
      # while (numel(userVote)==0){
      #   
      #   userVote = readline(prompt="Enter 1 to proceed with modeling, or 0 to skip to next site: ")
      #   
      #   if (!numel(userVote)==0 && userVote!=1 && userVote!=0){
      #     message('WARNING: Entry not allowed')
      #     userVote = c()
      #     
      #   }
      #   else if (numel(userVote)==0){
      #     message('WARNING: Entry not allowed')
      #   }
      #   
      # }
      
      #if (as.numeric(userVote)==1){ # only proceed if user says "yes" to modeling
      
      Presence = rep(0,length(numClicks))
      Presence[pres] = 1
      
      # Determine autocorrelation and create grouping variable for GEEGLMs
      corr = acf(Presence[1:50000],lag.max=1500,na.action=na.exclude,plot=FALSE) 
      lagID = which(abs(corr$acf)<0.2) # determine lag at which autocorrelation is <0.2
      #itsVal = IntegralTimeScaleCalc(numClicks[1:50000])
      numClust = length(numClicks)/(lagID[1]-1)
      if (numClust<length(numClicks)){
        clustID = rep(1:ceiling(numClust),each=lagID[1])
        clustID = clustID[1:numel(numClicks)]
      } else {
        clustID = 1:length(numClicks)
      }
      
      # remove days/hours with no effort, create Julian day and year variables
      noDat = which(is.na(numClicks)) 
      if (!numel(noDat) == 0) {
        numClicks = numClicks[-noDat]
        Presence = Presence[-noDat]
        reducedDateVec = dateVec[-noDat]
        reducedClustID = clustID[-noDat]
        Jday = as.numeric(format(reducedDateVec, "%j"))
        yearGroup = year(reducedDateVec)
      } else if (numel(noDat) == 0) {
        reducedDateVec = dateVec
        reducedClustID = clustID
        Jday = as.numeric(format(reducedDateVec, "%j"))
        yearGroup = year(reducedDateVec)
      }
      
      # code for study year
      yrCode = rep(NaN,length(reducedDateVec))
      yrCode[which(reducedDateVec < (dateStart + (60*60*24*365)))] = 1
      yrCode[which(reducedDateVec >= (dateStart + (60*60*24*365)))] = 2
      yrCode[which(reducedDateVec >= (dateStart + (60*60*24*365*2)))] = 3
      
      if (j == 7) {
        # create variable coding for change of site at HAT
        hatSite = rep(2, length(reducedDateVec))
        hatAdates = which(reducedDateVec <= as.POSIXlt("2017-02-06 00:00:00",
                                                       format="%Y-%m-%d %H:%M:%S",tz="GMT"))
        hatSite[hatAdates] = 1
      }
      
      # account for leap day in 2016, shift Julian days by value of 1 
      # (so we don't try modeling Julian day 366 based on 1 data point)
      leapIdx = which(yearGroup==2016)
      Jday[leapIdx] = Jday[leapIdx]-1
      
      # # round presence data back to integers so it's Poisson distributed again
      # numClicks = round(numClicks)
      
      # load normalized bin times
      whichNormTimeFile = str_detect(normTimeFileList,sites[j])
      normTimes = data.frame(read.csv(normTimeFileList[whichNormTimeFile]))  # load file
      normBinTimes = normTimes$Bin[-noDat]
      
      # get lunar illuminance data
      illum = getMoonIllumination(date=as.character(reducedDateVec),keep="fraction")
      illum = illum[,2]
      
      # get day phase data
      dayData = getSunlightTimes(date=seq.Date(as.Date(dateStart),as.Date(dateEnd),by=1),
                                 lat=lats[j],lon=lons[j],
                                 keep=c("nauticalDawn","sunrise","sunset","nauticalDusk"),
                                 tz="UTC")
      dayData = dayData[,4:7]
      phaseBins = as.POSIXlt(c(t(dayData)),format="%Y-%m-%d %H:%M:%S",tz="GMT")
      tooLate = which(phaseBins>=dateEnd)
      phaseBins = phaseBins[-tooLate] # cutting out some good bins??
      phaseBins = c(dateStart,phaseBins,dateEnd)
      phaseVec = rep(list("Night","Dawn","Day","Dusk"),length.out=length(phaseBins)-1)
      
      # find which day phase each bin falls into
      whichBin = histc(as.numeric(reducedDateVec),as.numeric(phaseBins))
      # create vector coding for day phase
      dayPhase = phaseVec[whichBin$bin]
      
      # stitch everything together in a data frame
      if (j!=7){
        master = cbind(Presence,as.character(reducedDateVec),reducedClustID,Jday,yearGroup,yrCode,normBinTimes,dayPhase,illum)
        colnames(master) = c("Presence","TimeStamp","GroupID","JulianDay","Year","StudyYear","NormTime","DayPhase","LunarIllum")
      } else if (j==7){
        master = cbind(Presence,as.character(reducedDateVec),reducedClustID,Jday,yearGroup,yrCode,normBinTimes,dayPhase,illum,hatSite)
        colnames(master) = c("Presence","TimeStamp","GroupID","JulianDay","Year","StudyYear","NormTime","DayPhase","LunarIllum","HATSite")
      }
      
      # save as a .csv
      saveName = paste(tsDir,'/',CTname,'_at_',sites[j],"_",int,'_Master.csv',sep="")
      write.csv(master,saveName,row.names=FALSE)
      #}
      
      ### Create separate data frames for modeling nighttime presence relative to lunar variables -----
      
      ## calculate proportion of each night with clicks
      nightInd = which(dayPhase=='Night')
      gaps = which(diff(nightInd,lag=1)>1)
      nightEndInd = nightInd[gaps]
      nightStInd = nightInd[gaps+1]
      nightEndInd = nightEndInd[2:length(nightEndInd)] # data starts during nighttime, incomplete data for that night
      nightEndInd = c(nightEndInd,nightInd[length(nightInd)])
      nightSt = reducedDateVec[nightStInd]
      nightEnd = reducedDateVec[nightEndInd]
      nightBins = as.POSIXct(c(t(cbind(nightSt,nightEnd))),tz="GMT",origin="1970-01-01")
      # find which 5-min bins fall into each night & day
      whichNight = histc(as.numeric(reducedDateVec),as.numeric(nightBins))
      # tally number of bins in each night (odd bins are night, even are day)
      nightBinCounts = whichNight$cnt[seq(1,length(whichNight$cnt),by=2)] 
      
      # count how many bins each night have presence
      nightDF = data.frame(Bin=unlist(whichNight$bin),Pres=Presence)
      nightPres = nightDF %>% 
        group_by(Bin) %>%
        summarise(numBins = sum(Pres))
      if (nightPres$Bin[1]==0){ # get rid of bins that didn't fall into a full night
        nightPres = nightPres[2:dim(nightPres)[1],]
      }
      nightPres = nightPres[seq(1,dim(nightPres)[1],b=2),] # odd bins are night, even are day
      # divide # presence bins by # bins in each night
      nightPresProp = nightPres$numBins/nightBinCounts
      
      # calculate average lunar illuminance each night
      illumDF = data.frame(Bin=unlist(whichNight$bin),Illum=illum)
      nightIllum = illumDF %>%
        group_by(Bin) %>%
        summarise(avgIllum = mean(Illum))
      if (nightIllum$Bin[1]==0){
        nightIllum = nightIllum[2:dim(nightIllum)[1],]
      }
      nightIllum = nightIllum[seq(1,dim(nightIllum)[1],by=2),]
      
      # get times when moon is up
      lunTime = getMoonTimes(date=seq.Date(as.Date(dateStart),as.Date(dateEnd+(60*60*24)),by=1),
                             lat=lats[j],lon=lons[j],
                             keep=c("rise","set"),
                             tz="UTC")
      lunR = lunTime[,4]; lunR = lunR[is.na(lunR)==0] # some NAs showing up, get rid of them
      lunS = lunTime[,5]; lunS = lunS[is.na(lunS)==0]
      moonBins = as.POSIXlt(sort(cbind(lunR,lunS)),tz="GMT",origin="1970-01-01")# moon rise/set bins
      firstRise = which(moonBins==lunR[1]) # want to start with moon rise
      moonBins = moonBins[firstRise:length(moonBins),]
      
      tooLate = which(moonBins>=dateEnd) # don't go beyond study end date
      if (!is.empty(tooLate)){
        moonBins = moonBins[-tooLate]
        moonBins = c(moonBins,endDate)}
      
      # find which 5-min bins are during moon periods
      whichMoon = histc(as.numeric(reducedDateVec),as.numeric(moonBins))
      whichMoon$bin[whichMoon$bin==0] = NA
      moonInd = which(whichMoon$bin%%2!=0) # odd bins are moon up periods
      
      # find which 5-min bins during moon up periods are also during nighttime
      whichNightMoon = histc(as.numeric(reducedDateVec[moonInd]),as.numeric(nightBins))
      
      # calculate proportion of each night moon is up
      moonProp = whichNightMoon$cnt[seq(1,length(whichNightMoon$cnt),by=2)]/nightBinCounts
      
      # Determine autocorrelation and create grouping variable for GEEGLM
      corr = acf(nightPresProp,lag.max=60,na.action=na.exclude,plot=FALSE)
      lagID = which(abs(corr$acf)<0.2) # determine lag at which autocorrelation is <0.2
      numClust = ceiling(length(nightPresProp)/(lagID[1]-1))
      if (numClust<length(nightPresProp)){
        clustID = rep(1:numClust,each=lagID[1])
        clustID = clustID[1:length(nightPresProp)]
      } else {
        clustID = 1:length(nightPresProp)
      }
      
      # combine lunar variables into master data frame
      masterLun = cbind(as.character(nightSt),as.character(nightEnd),nightPresProp,nightIllum$avgIllum,moonProp,clustID)
      colnames(masterLun) = c("NightStart","NightEnd","PropPres","AvgLunIllum","PropMoonUp","GroupID")
      
      # save as a .csv
      saveName = paste(tsDir,'/',CTname,'_at_',sites[j],"_",int,'_MasterLun.csv',sep="")
      write.csv(masterLun,saveName,row.names=FALSE)
      
    }
    
  }
  
}



