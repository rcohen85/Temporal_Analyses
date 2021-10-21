library(pracma)
library(stringr)
library(lubridate)
library(ggplot2)
library(conover.test)
library(dunn.test)

# KW test & boxplots to compare months ------------------------------------

inDir = 'I:/TimeSeries'
dielDir = 'I:/Diel Plots'

fileList = list.files(path=inDir,pattern='*.csv',
                      full.names=TRUE,recursive=FALSE,
                      include.dirs=FALSE,no..=TRUE)


for (i in 1:numel(fileList)){  # for each species file
  
  thisCT = data.frame(read.csv(fileList[i]))  # load file
  attribs = attributes(thisCT)  # find names of sites (all columns but first)
  sites = attribs$names[2:numel(attribs$names)]
  
  CTname = str_remove(fileList[i],paste(inDir,'/',sep="")) # get the species/CT name
  CTname = str_remove(CTname,".csv")
  
  dateVec = as.Date(thisCT$Date,format="%d-%b-%Y")
  
  # if it doesn't already exist, create directory to save figures
  if (!dir.exists(paste(dielDir,'/',CTname,sep=""))){
    dir.create(paste(dielDir,'/',CTname,sep=""))
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
      # pairComp = pairwise.wilcox.test(thisSite,monthGroup,p.adjust.method="BH") # don't use, doesn't use rankings or pool variance in same way as KW test
      pairCompCI = conover.test(thisSite,monthGroup,method='bonferroni',wrap=TRUE,
                              table=TRUE,alpha=0.05)
      pairCompDunn = dunn.test(thisSite,monthGroup,method='bonferroni',wrap=TRUE,
                           table=TRUE,alpha=0.05)
      
      # Identify which pairwise comparisons were significant
      sig = which(pairCompCI$P.adjusted<0.05)
      sigPairs = pairCompCI$comparisons[sig]
      sigPairs = str_replace(sigPairs," - "," ")
      sigPairs = as.numeric(sigPairs)
      
      # Plot data by month
      par(mfrow=c(2,1))
      
      boxplot(thisSite~monthGroup,main=paste(paste(CTname,'at',sites[j])),
              xlab=c("Month"),ylab=c("Daily Counts"))
      # p = ggplot(data.frame(thisSite, as.factor(monthGroup)),aes(x=monthGroup, y=thisSite))
      # p + geom_violin()
      
      
    }

  }
  
}



# KW test & boxplots to compare phases of the day -------------------------


