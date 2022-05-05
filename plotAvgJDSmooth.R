library(pracma)
library(ggplot2)

DFDir = 'J:/Chpt_2/TimeSeries_ScaledByEffortError'
modDir = 'J:/Chpt_2/ModelOutput'
dfList = list.files(path=DFDir,pattern=paste('*5minBin_MasterTempLun.csv',sep=""),
                    full.names=TRUE,recursive=FALSE,
                    include.dirs=FALSE,no..=TRUE)

species = list.dirs(modDir,recursive=FALSE)
species = species[-which(!is.na(str_match(species,"plots")))]
goodMods = list(c("HZ","BC","WC","HAT"), # only use models which met binned_resids criterion
             c("GS","BP","BS"),
             c('HZ','OC','NC','BC','WC','NFC','JAX'),
             c('HZ','OC','BC','WC'),
             c('HZ','OC','NC','WC','NFC','HAT','GS','BP','BS'),
             c('NC','NFC'),
             c('HZ','OC','NC','BC','WC','NFC','HAT','GS'),
             c('HZ','OC','NC','BC','WC','NFC','HAT','GS','BP','BS','JAX'),
             c('HZ','OC','NC','BC','WC','NFC'))

m=1
for ( i in c(2:3,5:8,10:numel(species))){
  
  JDFits = numeric()
  JDmodFiles = list.files(path=species[i],pattern="*Model_TempLun3.Rdata",
                          full.names=TRUE,recursive=FALSE,include.dirs=FALSE,no..=TRUE)
  CTname = str_remove(species[i],paste(modDir,'/',sep=""))
  sites = list()
  for (j in 1:numel(JDmodFiles)){
    
    site = str_remove(JDmodFiles[j],paste(modDir,"/",CTname,"/",sep=""))
    sites = c(sites,str_remove(site,"_5minBin_Model_TempLun3.Rdata"))
    
    if (any(str_detect(unlist(goodMods[m]),unlist(sites[j])))){
    load(JDmodFiles[j]) # load model
    
    # if model converged and JD was significant
    if (tempMod$geese$error==0 & PV$'p-value'[PV$'Variable'=="JDs"]<0.05){ 
    
      # find associated master dataframe
      thisSpec = which(str_detect(dfList,CTname))
      atSite = which(str_detect(dfList,unlist(sites[j])))
      thisModInd = intersect(thisSpec,atSite)
      thisSite = data.frame(read.csv(dfList[thisModInd])) # load master data frame
      
      JDInd = numeric()
      JDNTInd = numeric()
      JDInd = which(str_detect(names(tempMod$coefficients),"JDs")&!str_detect(names(tempMod$coefficients),"NTs"))
      JDNTInd = which(str_detect(names(tempMod$coefficients),"JDs")&str_detect(names(tempMod$coefficients),"NTs"))
      
      BootstrapParameters<-rmvnorm(10000, coef(tempMod), summary(tempMod)$cov.unscaled)
      JDayBootstrapCoefs<- BootstrapParameters[,JDInd]
      JDNTBootstrapCoefs<- BootstrapParameters[,JDNTInd]
      
      ## Recreate basis functions for smooths
      JDayForPlotting<- seq(min(thisSite$JulianDay), max(thisSite$JulianDay), length=1000)
      JDBasis<- mSpline(JDayForPlotting, 
                        knots=quantile(thisSite$JulianDay, probs=c(0.275,0.5,0.725)),
                        Boundary.knots=c(1,365),
                        periodic=T) 
      NTForPlotting<- seq(min(thisSite$NormTime), max(thisSite$NormTime), length=1000)
      NTBasis<- mSpline(NTForPlotting,  
                        knots=quantile(thisSite$NormTime, probs=c(0.275,0.5,0.725)),
                        Boundary.knots=c(-1,1),
                        periodic=T)
      
      if (!isempty(JDInd) & isempty(JDNTInd)){
        
        Fit<- JDBasis%*%(coef(tempMod)[JDInd]) # multiply basis functions by model coefficients to get values of spline at each X
        RealFit<- Fit+coef(tempMod)[1] # add intercept
        JDFits = cbind(JDFits,RealFit)
        
      } else if (!isempty(JDInd) & !isempty(JDNTInd)) {
        
        FitJD<- JDBasis%*%(coef(tempMod)[JDInd]) # multiply basis functions by model coefficients to get values of spline at each X
        IntPlotDF = numeric()
        IntCI = numeric()
        # Recreate interaction basis
        rowMat = numeric()
        IntBasis = numeric()
        for (l in 1:dim(JDBasis)[1]){ # for each row in JD basis
          for (k in 1:3){ # multiply by all values of NT
            rowMat = cbind(rowMat,NTBasis[,1]*JDBasis[l,k],NTBasis[,2]*JDBasis[l,k],NTBasis[,3]*JDBasis[l,k])
          }
          IntBasis = rbind(IntBasis,apply(rowMat,MARGIN=2,mean)) # average this JD across NT
          rowMat = numeric()
        }
        
        # Calculate JD smooth averaged across NT values
        IntAdjust = IntBasis%*%coef(tempMod)[JDNTInd]
        RealFit<- FitJD+coef(tempMod)[1]+IntAdjust # add intercept, interaction term
        JDFits = cbind(JDFits,RealFit)
        
      }
      
    }
  }
  }
  
  # Average smooths across sites for which they were significant
  regionalSmooth = apply(JDFits,MARGIN=1,mean)
  
  # Plot
  plotDF = data.frame(JD=JDayForPlotting,Fit=regionalSmooth)
  RS = ggplot(
  ) + geom_line(data=plotDF,aes(x=JD,y=Fit,),size=1,color="51D4A5"
  # ) + geom_ribbon(data=IntCI,aes(x=JulianDay,ymin=Ymin, ymax=Ymax),
  #                 fill="51D4A5",
  #                 alpha=0.2,
  #                 stat ="identity"
  ) + labs(x = "Julian Day",
           y = "s(JulianDay)"
  ) + ggtitle(paste('Average ',CTname,' Seasonal Pattern',sep="")
  ) + coord_cartesian(xlim=c(1,365)
  ) + scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),
                         label=c("J","F","M","A","M","J","J","A","S","O","N","D")
  ) + theme(legend.position="right",
            axis.line = element_line(size=0.2),
            panel.background = element_blank()
  )
  
  saveName = paste(modDir,'/',CTname,'_MeanJDPlot.png',sep="")
  ggsave(saveName,device="png", width=2, scale=4, height=1, units="in",dpi=600)
  while (dev.cur()>1) {dev.off()}
  
  m=m+1
}