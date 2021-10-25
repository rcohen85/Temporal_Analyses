library(pracma)
library(stringr)
library(lubridate)
library(ggplot2)
library(conover.test)
library(dunn.test)

# KW test & boxplots to compare months ------------------------------------

inDir = 'I:/TimeSeries'
seasDir = 'I:/Seasonal_Plots'
int = "Daily"

fileList = list.files(path=inDir,pattern=paste('*_',int,'.csv',sep=""),
                      full.names=TRUE,recursive=FALSE,
                      include.dirs=FALSE,no..=TRUE)


for (i in 1:numel(fileList)){  # for each species file
  
  thisCT = data.frame(read.csv(fileList[i]))  # load file
  attribs = attributes(thisCT)  # find names of sites (all columns but first)
  sites = attribs$names[2:numel(attribs$names)]
  
  CTname = str_remove(fileList[i],paste(inDir,'/',sep="")) # get the species/CT name
  CTname = str_remove(CTname,paste('_',int,'.csv',sep=""))
  
  dateVec = as.Date(thisCT$Date,format="%d-%b-%Y")
  
  # if it doesn't already exist, create directory to save figures
  if (!dir.exists(paste(seasDir,'/',CTname,sep=""))){
    dir.create(paste(seasDir,'/',CTname,sep=""))
  }
  
  for (j in 1:numel(sites)){ # for each site
    
    thisSite = as.numeric(thisCT[,j+1]) # get presence data for this site
    noDat = which(is.na(thisSite)) # remove days with no presence data
    thisSite = thisSite[-noDat]
    reducedDateVec = dateVec[-noDat]
    monthGroup = month(reducedDateVec) # find which month each day of data falls in
    
    pres = which(thisSite!=0)
    
    if (!numel(pres)==0){ # if there is any presence
      
      # test for significant differences between months
      kruskal.test(thisSite,monthGroup,na.action=na.pass) 
      pairCompCI = conover.test(thisSite,monthGroup,method='bonferroni',label=TRUE,
                                wrap=TRUE,table=TRUE,alpha=0.05)
      # pairCompDunn = dunn.test(thisSite,monthGroup,method='bonferroni',wrap=TRUE,
      #                      table=TRUE,alpha=0.05)
      CITable = matrix(NA,11,11)
      CITable[upper.tri(CITable,diag=TRUE)] = t(pairCompCI$P.adjusted)
      CITable = t(CITable)
      rownames(CITable) = c("10","11","12","2","3","4","5","6","7","8","9")
      colnames(CITable) = c("1","10","11","12","2","3","4","5","6","7","8")
      saveName = paste(seasDir,'/',CTname,'/',sites[j],"_Monthly_CI_Comparisons.csv",sep="")
      write.csv(CITable,file=saveName,row.names=TRUE)
      
      # Identify which pairwise comparisons were significant
      # sig = which(pairCompCI$P.adjusted<0.025)
      # sigPairs = pairCompCI$comparisons[sig]
      # sigNumbers = regmatches(sigPairs, gregexpr("[[:digit:]]+", sigPairs))
      # sigNumbers = as.numeric(unlist(sigNumbers))
      
      # Plot data by month
      saveName = paste(seasDir,'/',CTname,'/',sites[j],"_MonthlyBoxplot.png",sep="")
      png(saveName)
      boxplot(thisSite~monthGroup,main=paste(paste(CTname,'at',sites[j])),
              xlab=c("Month"),ylab=c("Daily Counts"))
      dev.off()
      
    }

  }
  
}



# KW test & boxplots to compare phases of the day -------------------------


