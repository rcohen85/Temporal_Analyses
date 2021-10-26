library(pracma)
library(stringr)

inDir = 'I:/TimeSeries'
distDir = 'I:/DataDists'
corrDir = 'I:/Autocorrelation'
int = 'Daily'

fileList = list.files(path=inDir,pattern=paste('*',int,'.csv',sep=""),
                      full.names=TRUE,recursive=FALSE,
                      include.dirs=FALSE,no..=TRUE)


# Plot presence data distributions and autocorrelation

for (i in 1:numel(fileList)){  # for each species' file
  
  thisCT = data.frame(read.csv(fileList[i]))  # load file
  attribs = attributes(thisCT)  # find names of sites (all columns but first)
  sites = attribs$names[2:numel(attribs$names)]
  
  CTname = str_remove(fileList[i],paste(inDir,'/',sep="")) # get the species/CT name
  CTname = str_remove(CTname,paste("_",int,".csv",sep=""))
  
  # if they don't already exist, create directories to save figures
  if (!dir.exists(paste(distDir,'/',CTname,sep=""))){
    dir.create(paste(distDir,'/',CTname,sep=""))
  }
  if (!dir.exists(paste(corrDir,'/',CTname,sep=""))){
    dir.create(paste(corrDir,'/',CTname,sep=""))
  }
  
  for (j in 1:numel(sites)){ # for each site
     # get presence data at this site and replace NaNs with NAs
    thisSite = as.numeric(thisCT[,j+1])
    noDat = which(is.na(thisSite))
    thisSite[noDat] = NA
    
    pres = which(thisSite!=0)
    
    if (!numel(pres)==0){ # if there is any presence
      
      # plot and save a histogram of the presence data
      saveName = paste(distDir,'/',CTname,'/',sites[j],"_",int,"_DataDist.png",sep="")
      png(saveName)
      hist(thisSite,main=paste(int,CTname,'at',sites[j]),xlab=c("Counts"))
      dev.off()
      
      # plot and save a plot of the autocorrelation
      saveName = paste(corrDir,'/',CTname,'/',sites[j],"_",int,"_Autocorr.png",sep="")
      png(saveName)
      acf(thisSite,lag.max=100,na.action=na.pass,
          main=paste(int,CTname,'at',sites[j]))
      dev.off()
      
    }  
    
  }
  
}

