library(tidyverse)
library(lubridate)
library(pracma)
library(installr)
library(suncalc)
library(pracma)

tsDir = 'J:/Chpt_2/TimeSeries_ScaledByEffortError'
seasDir = 'J:/Chpt_2/Seasonal_Plots'
dielDir = 'J:/Chpt_2/Diel Plots'
illumDir = 'J:/Chpt_2/IlluminationFiles'
normTimeDir = 'J:/Chpt_2/NormTimes'
int = "5minBin"
start = '2016-05-01 00:00:00'
end = '2019-04-30 23:59:59'
goodIdx = c(1,2,4,8,11,13:19)
goodSite = list(c(8:10),c(8:10),c(),c(1:8,10),c(),c(),c(),c(7:11),c(),c(),c(1:11),
                c(),c(1:6),c(1:11),c(1:6),
                c(1:11),c(1:11),c(1:11),c(1:6,8,11))



lats = as.numeric(c(41.06165,40.22999,39.83295,
                    39.19192,38.37337,37.16452,
                    35.5841,33.66992,32.10527,
                    30.58295,30.27818))
lons = as.numeric(c(-66.35155,-67.97798,-69.98194,
                    -72.22735,-73.36985,-74.46585,
                    -74.7499,-75.9977,-77.09067,
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
# lunFileList = list.files(path=illumDir,pattern='*.csv',
#                          full.names=TRUE,recursive=FALSE,
#                          include.dirs=FALSE,no..=TRUE)

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
    
    if (numel(pres)>=100){ # if there are at least 100 bins with presence
      
      Presence = rep(0,length(numClicks))
      Presence[pres] = 1
      
      # Determine autocorrelation and create grouping variable for GEEGLMs
      corr = acf(Presence[1:150000],lag.max=1500,na.action=na.exclude,plot=FALSE) 
      lagID = which(abs(corr$acf)<0.2) # determine lag at which autocorrelation is <0.2
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
      
      # # account for leap day in 2016, shift Julian days by value of 1 
      # # (so we don't try modeling Julian day 366 based on 1 data point)
      # leapIdx = which(yearGroup==2016)
      # Jday[leapIdx] = Jday[leapIdx]-1
      
      # # round presence data back to integers so it's Poisson distributed again
      # numClicks = round(numClicks)
      
      # load normalized bin times
      whichNormTimeFile = str_detect(normTimeFileList,sites[j])
      normTimes = data.frame(read.csv(normTimeFileList[whichNormTimeFile]))  # load file
      normBinTimes = normTimes$Bin[-noDat]
      
      # # get lunar illuminance data
      # illum = getMoonIllumination(date=as.character(reducedDateVec),keep="fraction")
      # illum = illum[,2]
      
      # get day phase data
      dayData = getSunlightTimes(date=seq.Date(as.Date(dateStart),as.Date(dateEnd),by=1),
                                 lat=lats[j],lon=lons[j],
                                 keep=c("nauticalDawn","sunrise","sunset","nauticalDusk"),
                                 tz="UTC")
      dayData = dayData[,4:7]
      phaseBins = as.POSIXct(c(t(dayData)),format="%Y-%m-%d %H:%M:%S",tz="GMT")
      tooLate = which(phaseBins>=dateEnd)
      phaseBins = phaseBins[-tooLate]
      phaseBins = c(dateStart,phaseBins,dateEnd)
      phaseVec = rep(list("Night","Dawn","Day","Dusk"),length.out=length(phaseBins)-1)
      
      # find which day phase each bin falls into
      whichBin = histc(as.numeric(reducedDateVec),as.numeric(phaseBins))
      # create vector coding for day phase
      dayPhase = phaseVec[whichBin$bin]
      
      # get times when moon is up
      lunTime = getMoonTimes(date=seq.Date(as.Date(dateStart),as.Date(dateEnd+(60*60*24)),by=1),
                             lat=lats[j],lon=lons[j],
                             keep=c("rise","set"),
                             tz="UTC")
      lunR = lunTime[,4]; lunR = lunR[!is.na(lunR)] # some NAs showing up, get rid of them
      lunS = lunTime[,5]; lunS = lunS[!is.na(lunS)]
      
      # determine if each time bin is before, during, or after moon presence
      BinDays = floor_date(reducedDateVec,unit="days")
      moonRiseDays = floor_date(lunR,unit="days")
      moonSetDays = floor_date(lunS,unit="days")
      allDays = unique(BinDays)
      moonPres = matrix(nrow=length(reducedDateVec),ncol=1)
      
      for (k in allDays){
        thisRise = lunR[which(moonRiseDays==k)] # time of this day's moon rise
        thisSet = lunS[which(moonSetDays==k)] # time of this day's moon set
        
        thisDay = which(BinDays==k) # which bins fall on this day
        
        if (length(thisRise)>0 & length(thisSet)>0){ # if the moon both rises and sets in this day
          # which bins on this day are before moon rise
          preMoon = which(reducedDateVec[thisDay]<thisRise)
          moonPres[thisDay[preMoon]] = "Pre"
          # which bins on this day are during moon up
          moonUp = which(reducedDateVec[thisDay]>=thisRise & reducedDateVec[thisDay]<=thisSet)
          moonPres[thisDay[moonUp]] = "MoonUp"
          # which bins on this day are after moon set
          postMoon = which(reducedDateVec[thisDay]>thisSet)
          moonPres[thisDay[postMoon]] = "Post"
        } else if (length(thisRise)>0 & length(thisSet)==0){ # if the moon rises in this day, but doesn't set (sets after midnight)
          # which bins on this day are before moon rise
          preMoon = which(reducedDateVec[thisDay]<thisRise)
          moonPres[thisDay[preMoon]] = "Pre"
          # remainder of bins are during moon up
          moonUp = setdiff(1:length(thisDay),preMoon)
          moonPres[thisDay[moonUp]] = "MoonUp"
        }  else if (length(thisRise)==0 & length(thisSet)>0){ # if the moon doesn't rise in this day, but does set (rose prior to midnight)
          # which bins on this day are after moon set
          postMoon = which(reducedDateVec[thisDay]>thisSet)
          moonPres[thisDay[postMoon]] = "Post"
          # remainder of bins are during moon up
          moonUp = setdiff(1:length(thisDay),postMoon)
          moonPres[thisDay[moonUp]] = "MoonUp"
        } else if (length(thisRise)==0 & length(thisSet)==0){ # if the moon doesn't rise or set in this day? 
          moonPres[thisDay] = "UhOh"
        }
      }
      
      # get phase of moon; data are for all bins, not just when moon is up
      allMoonPhase = getMoonIllumination(date=as.character(reducedDateVec),keep="phase")
      moonPhase = allMoonPhase$phase
      
      # create factor for night + moonUp
      LunFact = rep(0,length(reducedDateVec)) 
      LunFact[normBinTimes>0 & moonPres=="MoonUp"] = 1
      
      # stitch everything together in a data frame
      if (j!=7){
        master = cbind(Presence,
                       as.character(reducedDateVec),
                       reducedClustID,
                       Jday,yearGroup,
                       yrCode,
                       normBinTimes,
                       dayPhase,
                       moonPhase,
                       moonPres,
                       LunFact)
        colnames(master) = c("Presence",
                             "TimeStamp",
                             "GroupID",
                             "JulianDay",
                             "Year",
                             "StudyYear",
                             "NormTime",
                             "DayPhase",
                             "MoonPhase",
                             "MoonPres",
                             "LunFact")
      } else if (j==7){
        master = cbind(Presence,
                       as.character(reducedDateVec),
                       reducedClustID,
                       Jday,
                       yearGroup,
                       yrCode,
                       normBinTimes,
                       dayPhase,
                       moonPhase,
                       moonPres,
                       LunFact,
                       hatSite)
        colnames(master) = c("Presence",
                             "TimeStamp",
                             "GroupID",
                             "JulianDay",
                             "Year",
                             "StudyYear",
                             "NormTime",
                             "DayPhase",
                             "MoonPhase",
                             "MoonPres",
                             "LunFact",
                             "HATSite")
      }
      
      # save as a .csv
      saveName = paste(tsDir,'/',CTname,'_at_',sites[j],"_",int,'_MasterTempLun.csv',sep="")
      write.csv(master,saveName,row.names=FALSE)
      
      ### Create separate data frames for modeling nighttime presence relative to lunar variables -----
      
      # # get nighttime bins (dusk/night/dawn)
      # nb = which(dayPhase=="Dusk" | dayPhase=="Night" | dayPhase=="Dawn")
      # nightBins = reducedDateVec[nb]
      # nightPresence = Presence[nb]
      # nightPhase = dayPhase[nb]
      # 
      # # get times when moon is up
      # lunTime = getMoonTimes(date=seq.Date(as.Date(dateStart),as.Date(dateEnd+(60*60*24)),by=1),
      #                        lat=lats[j],lon=lons[j],
      #                        keep=c("rise","set"),
      #                        tz="UTC")
      # lunR = lunTime[,4]; lunR = lunR[!is.na(lunR)] # some NAs showing up, get rid of them
      # lunS = lunTime[,5]; lunS = lunS[!is.na(lunS)]
      # 
      # # determine if each time bin is before, during, or after moon presence
      # nightBinDays = floor_date(nightBins,unit="days")
      # moonRiseDays = floor_date(lunR,unit="days")
      # moonSetDays = floor_date(lunS,unit="days")
      # allDays = unique(nightBinDays)
      # moonPres = matrix(nrow=length(nightBins),ncol=1)
      # 
      # for (k in allDays){
      #   thisRise = lunR[which(moonRiseDays==k)] # time of this day's moon rise
      #   thisSet = lunS[which(moonSetDays==k)] # time of this day's moon set
      #   
      #   thisDay = which(nightBinDays==k) # which bins fall on this day
      #   
      #   if (length(thisRise)>0 & length(thisSet)>0){ # if the moon both rises and sets in this night
      #     # which bins on this day are before moon rise
      #     preMoon = which(nightBins[thisDay]<thisRise)
      #     moonPres[thisDay[preMoon]] = "Pre"
      #     # which bins on this day are during moon up
      #     moonUp = which(nightBins[thisDay]>=thisRise & nightBins[thisDay]<thisSet)
      #     moonPres[thisDay[moonUp]] = "MoonUp"
      #     # which bins on this day are after moon set
      #     postMoon = which(nightBins[thisDay]>thisSet)
      #     moonPres[thisDay[postMoon]] = "Post"
      #   } else if (length(thisRise)>0 & length(thisSet)==0){ # if the moon rises in this night, but doesn't set (sets after sunrise)
      #     # which bins on this day are before moon rise
      #     preMoon = which(nightBins[thisDay]<thisRise)
      #     moonPres[thisDay[preMoon]] = "Pre"
      #     # remainder of bins are during moon up
      #     moonUp = setdiff(1:length(thisDay),preMoon)
      #     moonPres[thisDay[moonUp]] = "MoonUp"
      #   }  else if (length(thisRise)==0 & length(thisSet)>0){ # if the moon doesn't rise in this night, but does set (rose prior to sunset)
      #     # which bins on this day are after moon set
      #     postMoon = which(nightBins[thisDay]>thisSet)
      #     moonPres[thisDay[postMoon]] = "Post"
      #     # remainder of bins are during moon up
      #     moonUp = setdiff(1:length(thisDay),postMoon)
      #     moonPres[thisDay[moonUp]] = "MoonUp"
      #   } else if (length(thisRise)==0 & length(thisSet)==0){ # if the moon doesn't rise or set in this night? 
      #     moonPres[thisDay] = "UhOh"
      #   }
      # }
      # 
      # # get phase of moon; data are for all night bins, not just when moon is up
      # allMoonPhase = getMoonIllumination(date=as.character(nightBins),keep="phase")
      # moonPhase = allMoonPhase$phase
      #
      # # get altitude of moon; data are for all night bins, not just when moon is up
      # allMoonAlt = getMoonPosition(date=as.character(nightBins),lat=lats[j],lon=lons[j],keep="altitude")
      # moonEl = allMoonAlt$altitude
      # 
      # # load apparent magnitude and elevation data (downloaded from Tethys)
      # # Note: data are only for when the moon is up, # bins in these files will not match # nightBins
      # lunFile = str_which(lunFileList,sites[j])
      # magElDat = data.frame(read.csv(lunFileList[lunFile]))
      # keepInd = which(as.POSIXct(magElDat$Date,format="%d-%b-%Y %H:%M:%S",tz="GMT") %in% nightBins)
      # putWhere = match(as.POSIXct(magElDat$Date,format="%d-%b-%Y %H:%M:%S",tz="GMT"),nightBins)
      # keepInd = keepInd[!is.na(keepInd)]
      # putWhere = putWhere[!is.na(putWhere)]
      # moonMag = matrix(nrow=length(nightBins),ncol=1)
      # moonMag[putWhere] = magElDat$Magnitude[keepInd]
      # # moonEl = matrix(nrow=length(nightBins),ncol=1)
      # # moonEl[putWhere] = magElDat$Elevation[keepInd]
      # 
      # # some bins missing in magEl data, interpolate
      # magTimeDiff = diff(as.POSIXct(magElDat$Date,format="%d-%b-%Y %H:%M:%S",tz="GMT"))
      # skippedBins = which(magTimeDiff==dminutes(x=10)) # time gaps
      # keepSkippedBins = intersect(keepInd,skippedBins) # time gaps in night bins that we care about
      # missMag = apply(cbind(magElDat$Magnitude[keepSkippedBins],magElDat$Magnitude[keepSkippedBins+1]),MARGIN=1,mean)
      # # missEl = apply(cbind(magElDat$Elevation[keepSkippedBins],magElDat$Elevation[keepSkippedBins+1]),MARGIN=1,mean)
      # missBins = as.POSIXct(magElDat$Date[keepSkippedBins],format="%d-%b-%Y %H:%M:%S",tz="GMT")+dminutes(x=5)
      # putWhere2 = match(missBins,nightBins) # indices where data is missing
      # moonMag[putWhere2] = missMag
      # # moonEl[putWhere2] = missEl
      # 
      # # # make sure no magnitude or elevation data when moon is not up
      # # moonMag[moonPres=="Pre"|moonPres=="Post"] = NA
      # # moonEl[moonPres=="Pre"|moonPres=="Post"] = NA
      # 
      # # Determine autocorrelation and create grouping variable for GEEGLM
      # corr = acf(nightPresence,lag.max=50000,na.action=na.exclude,plot=FALSE)
      # lagID = which(abs(corr$acf)<0.2) # determine lag at which autocorrelation is <0.2
      # numClust = ceiling(length(nightPresence)/(lagID[1]-1))
      # if (numClust<length(nightPresence)){
      #   clustID = rep(1:numClust,each=lagID[1])
      #   clustID = clustID[1:length(nightPresence)]
      # } else {
      #   clustID = 1:length(nightPresence)
      # }
      # 
      # # combine lunar variables into master data frame
      # masterLun = cbind(nightPresence,as.character(nightBins),nightPhase,moonPres,moonPhase,moonMag,moonEl,clustID)
      # colnames(masterLun) = c("NightPres","NightBinTimes","NightPhase","MoonPres","MoonPhase","MoonMag","MoonAltitude","GroupID")
      # 
      # # save as a .csv
      # saveName = paste(tsDir,'/',CTname,'_at_',sites[j],"_",int,'_MasterLun.csv',sep="")
      # write.csv(masterLun,saveName,row.names=FALSE)
      
    }
    
  }
  
}



##### Old code for lunar data frames

# ## calculate proportion of each night with clicks
# nightInd = which(dayPhase=='Night')
# gaps = which(diff(nightInd,lag=1)>1)
# nightEndInd = nightInd[gaps]
# nightStInd = nightInd[gaps+1]
# nightEndInd = nightEndInd[2:length(nightEndInd)] # data starts during nighttime, incomplete data for that night
# nightEndInd = c(nightEndInd,nightInd[length(nightInd)])
# nightSt = reducedDateVec[nightStInd]
# nightEnd = reducedDateVec[nightEndInd]
# nightBins = as.POSIXct(c(t(cbind(nightSt,nightEnd))),tz="GMT",origin="1970-01-01")
# # find which 5-min bins fall into each night & day
# whichNight = histc(as.numeric(reducedDateVec),as.numeric(nightBins))
# # tally number of bins in each night (odd bins are night, even are day)
# nightBinCounts = whichNight$cnt[seq(1,length(whichNight$cnt),by=2)] 
# 
# # count how many bins each night have presence
# nightDF = data.frame(Bin=unlist(whichNight$bin),Pres=Presence)
# nightPres = nightDF %>% 
#   group_by(Bin) %>%
#   summarise(numBins = sum(Pres))
# if (nightPres$Bin[1]==0){ # get rid of bins that didn't fall into a full night
#   nightPres = nightPres[2:dim(nightPres)[1],]
# }
# nightPres = nightPres[seq(1,dim(nightPres)[1],b=2),] # odd bins are night, even are day
# # divide # presence bins by # bins in each night
# nightPresProp = nightPres$numBins/nightBinCounts
# 
# # calculate average lunar illuminance each night
# illumDF = data.frame(Bin=unlist(whichNight$bin),Illum=illum)
# nightIllum = illumDF %>%
#   group_by(Bin) %>%
#   summarise(avgIllum = mean(Illum))
# if (nightIllum$Bin[1]==0){
#   nightIllum = nightIllum[2:dim(nightIllum)[1],]
# }
# nightIllum = nightIllum[seq(1,dim(nightIllum)[1],by=2),]
# 
# # get times when moon is up
# lunTime = getMoonTimes(date=seq.Date(as.Date(dateStart),as.Date(dateEnd+(60*60*24)),by=1),
#                        lat=lats[j],lon=lons[j],
#                        keep=c("rise","set"),
#                        tz="UTC")
# lunR = lunTime[,4]; lunR = lunR[is.na(lunR)==0] # some NAs showing up, get rid of them
# lunS = lunTime[,5]; lunS = lunS[is.na(lunS)==0]
# moonBins = as.POSIXlt(sort(cbind(lunR,lunS)),tz="GMT",origin="1970-01-01")# moon rise/set bins
# firstRise = which(moonBins==lunR[1]) # want to start with moon rise
# moonBins = moonBins[firstRise:length(moonBins),]
# 
# tooLate = which(moonBins>=dateEnd) # don't go beyond study end date
# if (!is.empty(tooLate)){
#   moonBins = moonBins[-tooLate]
#   moonBins = c(moonBins,endDate)}
# 
# # find which 5-min bins are during moon periods
# whichMoon = histc(as.numeric(reducedDateVec),as.numeric(moonBins))
# whichMoon$bin[whichMoon$bin==0] = NA
# moonInd = which(whichMoon$bin%%2!=0) # odd bins are moon up periods
# 
# # find which 5-min bins during moon up periods are also during nighttime
# whichNightMoon = histc(as.numeric(reducedDateVec[moonInd]),as.numeric(nightBins))
# 
# # calculate proportion of each night moon is up
# moonProp = whichNightMoon$cnt[seq(1,length(whichNightMoon$cnt),by=2)]/nightBinCounts
# 
# # Determine autocorrelation and create grouping variable for GEEGLM
# corr = acf(nightPresProp,lag.max=90,na.action=na.exclude,plot=FALSE)
# lagID = which(abs(corr$acf)<0.2) # determine lag at which autocorrelation is <0.2
# numClust = ceiling(length(nightPresProp)/(lagID[1]-1))
# if (numClust<length(nightPresProp)){
#   clustID = rep(1:numClust,each=lagID[1])
#   clustID = clustID[1:length(nightPresProp)]
# } else {
#   clustID = 1:length(nightPresProp)
# }
# 
# # combine lunar variables into master data frame
# masterLun = cbind(as.character(nightSt),as.character(nightEnd),nightPresProp,nightIllum$avgIllum,moonProp,clustID)
# colnames(masterLun) = c("NightStart","NightEnd","PropPres","AvgLunIllum","PropMoonUp","GroupID")
# 
# # save as a .csv
# saveName = paste(tsDir,'/',CTname,'_at_',sites[j],"_",int,'_MasterLun.csv',sep="")
# write.csv(masterLun,saveName,row.names=FALSE)
