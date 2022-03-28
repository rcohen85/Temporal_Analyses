# Calculate how much presence each species has at each site, as an absolute count 
# of 5-min bins, and as a proportion of effort

library(pracma)

modelDFDir = 'J:/Chpt_2/TimeSeries_ScaledByEffortError'
outDir = 'J:/Chpt_2/ModelOutput'
int = "5minBin"

dfList = list.files(path=modelDFDir,pattern=paste('*',int,'.csv',sep=""),
                    full.names=TRUE,recursive=FALSE,
                    include.dirs=FALSE,no..=TRUE)
# lunList = list.files(path=modelDFDir,pattern=paste('*',int,'_MasterLun.csv',sep=""),
#                      full.names=TRUE,recursive=FALSE,
#                      include.dirs=FALSE,no..=TRUE)
sites = c("HZ","OC","NC","BC","WC","NFC","HAT","GS","BP","BS","JAX")

totPres = data.frame(matrix(ncol=0,nrow=11))
rownames(totPres) = sites

k=1 # counter to rename columns
for (i in c(1,2,4,8,11,13:(length(dfList)-2))){
  
  thisSpecies = data.frame(read.csv(dfList[i]))
  fileName = str_remove(dfList[i],paste(modelDFDir,'/',sep="")) # get the species/CT name
  CTname = str_remove(fileName,paste("_",int,".csv",sep=""))

  if (str_detect(CTname,"Gervais") & str_detect(CTname,"Atl")){
    CTname = "Gervais"
  }
  
  eval(parse(text=paste("totPres$",CTname,"Count=as.numeric(rep(NA,11,1))",sep="")))
  eval(parse(text=paste("totPres$",CTname,"Prop=as.numeric(rep(NA,11,1))",sep="")))
  
  for (j in 2:dim(thisSpecies)[2]){
    
    if (any(i==c(11,16:19))){
    count = sum((thisSpecies[,j]>50),na.rm=TRUE)
    } else {
      count = sum((thisSpecies[,j]>20),na.rm=TRUE)
    }
    prop = round(count/length(thisSpecies[,j]),digits=4)
    
    eval(parse(text=paste("totPres$",CTname,"Count[j-1]=",count,sep="")))
    eval(parse(text=paste("totPres$",CTname,"Prop[j-1]=",prop,sep="")))
  }
  
  k = k+2
  
  if (i==1){
    rownames(totPres) = colnames(thisSpecies)[2:dim(thisSpecies)[2]]
  }
}

write.csv(totPres,file=paste(outDir,'/PresenceQuantification.csv',sep=""),row.names=TRUE)
