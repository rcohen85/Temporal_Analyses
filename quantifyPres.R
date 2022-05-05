# Calculate how much presence each species has at each site, as an absolute count 
# of 5-min bins, and as a proportion of effort
library(stringr)
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

for (i in c(1,2,4,8,11,13:15,17:19)){
  
  thisSpecies = data.frame(read.csv(dfList[i]))
  fileName = str_remove(dfList[i],paste(modelDFDir,'/',sep="")) # get the species/CT name
  CTname = str_remove(fileName,paste("_",int,".csv",sep=""))

  if (str_detect(CTname,"Gervais") & str_detect(CTname,"Atl")){
    CTname = "Gervais"
  }
  
  if (i==1){
    totPres$EffortBins = as.numeric(rep(NA,11,1))
  }
  
  
  eval(parse(text=paste("totPres$",CTname,"Count=as.numeric(rep(NA,11,1))",sep="")))
  eval(parse(text=paste("totPres$",CTname,"Prop=as.numeric(rep(NA,11,1))",sep="")))
  
  for (j in 2:dim(thisSpecies)[2]){
    
    if (any(i==c(11,16:19))){
    count = sum((thisSpecies[,j]>50),na.rm=TRUE)
    } else {
      count = sum((thisSpecies[,j]>20),na.rm=TRUE)
    }
    effBins = sum(!is.na(thisSpecies[,j]))
    prop = round(count/effBins,digits=4)
    
    totPres$EffortBins[j-1] = effBins
    eval(parse(text=paste("totPres$",CTname,"Count[j-1]=",count,sep="")))
    eval(parse(text=paste("totPres$",CTname,"Prop[j-1]=",prop,sep="")))

  }
  
  if (i==1){
    rownames(totPres) = colnames(thisSpecies)[2:dim(thisSpecies)[2]]
  }

}

write.csv(totPres,file=paste(outDir,'/PresenceQuantification.csv',sep=""),row.names=TRUE)


# Plot presence for each species as a bar chart

countInd = which(str_detect(colnames(totPres),"Prop"))
plotList = list()
# xmax = max(apply(totPres[,countInd],MARGIN=2,max))*1.05

for (i in 1:length(countInd)){
  
  plotDF = data.frame(sites=sites,pres=totPres[,countInd[i]]*1000)
  plotDF$pres[plotDF$pres<1] = NA
  plotDF$pres = log10(plotDF$pres)
  
  spec = str_remove(colnames(totPres)[countInd[i]],"Prop")
  if (spec=="Risso"){ylabs=rev(sites)}else{ylabs=NULL}
  
  plotList[[spec]] = ggplot(plotDF
                            )+geom_col(aes(y=sites,x=pres),fill="#66B2FF",alpha=0.6
                            )+scale_y_discrete(limits=rev(sites),labels=ylabs
                            )+coord_cartesian(xlim=c(0,2.5)
                            )+scale_x_continuous(breaks=c(0,1,2),
                                                 labels=c(1,10,100)
                            )+annotation_logticks(sides="b",
                                                  scaled=TRUE,
                                                  alpha=0.4,
                                                  short=unit(0.5,"mm"),
                                                  mid=unit(1.2,"mm"),
                                                  long=unit(2,"mm")
                            )+labs(x=NULL,y=NULL,title=spec
                            )+theme_minimal(
                            )+theme(plot.title=element_text(size=12),
                                    panel.grid.minor.x=element_blank())
}

plotList =plotList[c("Risso","UD36","UD28","UD26","Blainville","Cuvier","Gervais","Sowerby","True","Kogia","SpermWhale")]

png(file=paste(outDir,'/PropPresComp.png',sep=""),width = 1600, height = 400, units = "px",res=125)
grid.arrange(grobs=plotList,ncol=11,nrow=1,bottom="Bins per Thousand with Presence",left="Site")
while (dev.cur()>1) {dev.off()}


