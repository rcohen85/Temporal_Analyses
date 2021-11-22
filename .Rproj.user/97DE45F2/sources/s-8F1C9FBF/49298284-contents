

source("getTwilight.r")

inDir = 'I:/TimeSeries_ScaledByEffortError'
dielDir = 'I:/Diel Plots'
int = "5minBin"

datFileList = list.files(path=datDir,pattern=paste('*_',int,'.csv',sep=""),
                         full.names=TRUE,recursive=FALSE,
                         include.dirs=FALSE,no..=TRUE)
normTimeFileList = list.files(path=normTimeDir,pattern='*.csv',
                              full.names=TRUE,recursive=FALSE,
                              include.dirs=FALSE,no..=TRUE)
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

site = c("WAT_HZ","WAT_OC","WAT_NC","WAT_BC","WAT_WC","NFC",
  "HAT","WAT_GS","WAT_BP","WAT_BS","JAX")

for (i in goodIdx){  # for each species' file
  
  thisCT = data.frame(read.csv(datfileList[i]))  # load file
  attribs = attributes(thisCT)  # find names of sites (all columns but first)
  sites = attribs$names[-1]
  
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
      
      # load normalized bin times
      whichNormTimeFile = str_detect(normTimeFileList,sites[j])
      normTimes = data.frame(read.csv(normTimeFileList[whichNormTimeFile]))  # load file
      normBinTimes = normTimes$Bin[-noDat]
      
      # load lunar illuminance data
  
      # get day phase data
       TWdata = getSunlightTimes(date=dateRange,
                            lat=lat[i],lon=lon[i],
                            keep=c("nauticalDawn","nauticalDusk"),
                            tz="UTC")
       
       # determine which day phase each bin fall into
  
  # GEEGLM: # clicks per bin ~ s(norm time of day) + as.factor(day phase) + S(lunar illuminance)
  
    }
  }
}


