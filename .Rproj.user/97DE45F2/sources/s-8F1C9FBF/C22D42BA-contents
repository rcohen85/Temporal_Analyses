library(tidyverse)
library(lubridate)
library(pracma)

tsDir = 'I:/TimeSeries_ScaledByEffortError'
seasDir = 'I:/Seasonal_Plots'
dielDir = 'I:/Diel Plots'
illumDir = 'I:/IlluminationFiles'
normTimeDir = 'I:/NormTimes'
int = "5minBin"
start = '2016-05-01 00:00:00'
end = '2019-04-30 23:59:59'
# goodIdx = c(1,2,4,8,11,13:19)
goodIdx = c(14:19)


lats = as.numeric(c(41.06165,40.22999,39.83295,
                    39.19192,38.37337,37.16452,
                    35.30183,33.66992,32.10527,
                    30.58295,30.27818))
lons = as.numeric(c(-66.35155,-67.97798,-69.98194,
                    -72.22735,-73.36985,-74.46585,
                    -74.87895,-75.9977,-77.09067,
                    -77.39002,-80.22085))


# Cycle through time series and organize covariates for each species at each site --------

dateStart = as.POSIXlt(start,format="%Y-%m-%d %H:%M:%S",tz="GMT");
dateEnd = as.POSIXlt(end,format="%Y-%m-%d %H:%M:%S",tz="GMT");
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
  
  #dateVec = as.POSIXlt(thisCT$Date,tz="GMT",format="%d-%b-%Y")
  if (int=="Daily") {
    dateVec = as.POSIXlt(thisCT$Date,format="%d-%b-%Y",tz="GMT")
  } else if (int=="Hourly") {
    dateVec = as.POSIXlt(thisCT$Hour,format="%d-%b-%Y %H:%M:%S",tz="GMT")
  } else if (int=="5minBin") {
    dateVec = as.POSIXlt(thisCT$Bin,format="%d-%b-%Y %H:%M:%S",tz="GMT")
  }
  
  
  # if it doesn't already exist, create directory to save figures
  if (!dir.exists(paste(seasDir,'/',CTname,sep=""))){
    dir.create(paste(seasDir,'/',CTname,sep=""))
  }
  
  for (j in 1:numel(sites)){ # for each site
    
    numClicks = as.numeric(thisCT[,j+1]) # get presence data for this site
    pres = which(numClicks!=0)
    
    if (!numel(pres)==0){ # if there is any presence
      
      # plot presence to determine if it's worth modeling
      plot(dateVec,numClicks,main=paste(CTname,'at',sites[j],sep=" "),
           ylab="# 5-min Bins w Presence")
      
      userVote = c()
      while (numel(userVote)==0){
        
        userVote = readline(prompt="Enter 1 to proceed with modeling, or 0 to skip to next site: ")
        
        if (!numel(userVote)==0 && userVote!=1 && userVote!=0){
          message('WARNING: Entry not allowed')
          userVote = c()
          
        }
        else if (numel(userVote)==0){
          message('WARNING: Entry not allowed')
        }
        
      }
      
      if (as.numeric(userVote)==1){ # only proceed if user says "yes" to modeling
        
        Presence = rep(0,length(numClicks))
        Presence[pres] = 1
      
        # Determine autocorrelation and create grouping variable for GEEGLMs
        corr = acf(numClicks,lag.max=1500,na.action=na.exclude,plot=FALSE) 
        lagID = which(abs(corr$acf)<0.2) # determine lag at which autocorrelation is <0.2
        itsVal = IntegralTimeScaleCalc(numClicks)
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
        
        # round presence data back to integers so it's Poisson distributed again
        numClicks = round(numClicks)
        
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
        
        # find which day phase each presence bin falls into
        whichBin = histc(as.numeric(reducedDateVec),as.numeric(phaseBins))
        # create vector coding for day phase
        dayPhase = phaseVec[whichBin$bin]
        
        # stitch everything together in a data frame
        if (j!=7){
          master = cbind(Presence,numClicks,Jday,yearGroup,normBinTimes,dayPhase,illum)
          colnames(master) = c("Presence","NumClicks","JulianDay","Year","NormTime","DayPhase","LunarIllum")
        } else if (j==7){
          master = cbind(Presence,numClicks,Jday,yearGroup,normBinTimes,dayPhase,illum,hatSite)
          colnames(master) = c("Presence","NumClicks","JulianDay","Year","NormTime","DayPhase","LunarIllum","HATSite")
        }
        
        # save as a .csv
        saveName = paste(tsDir,'/',CTname,'_at_',sites[j],"_",int,'_Master.csv',sep="")
        write.csv(master,saveName,row.names=FALSE)
      }
    }
    
  }
  
}

