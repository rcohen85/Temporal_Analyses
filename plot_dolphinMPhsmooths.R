library(ggplot2)
library(SimDesign)
library(splines2)
library(pracma)
library(stringr)
library(gridExtra)

modelDFDir = 'J:/Chpt_2/TimeSeries_ScaledByEffortError'
modDir = 'J:/Chpt_2/ModelOutput'
saveDir = 'J:/Chpt_2/Figures'
int = "5minBin"

specs = c("Risso","UD36","UD28","UD26")
specTitles = c('Gg1',"Gg2","Dd","Gm")
sites = list(c("HZ","OC","NC","BC","WC","NFC","HAT","GS","BP","BS","JAX"),
             c("HZ","OC","NC","BC","WC","NFC"),
             c("HZ","OC","NC","BC","WC","NFC","HAT","GS","BP","BS","JAX"),
             c("HZ","OC","NC","BC","WC","NFC","HAT",'GS'))

dfList = list.files(path=modelDFDir,pattern=paste('*',int,'_MasterTempLun.csv',sep=""),
                    full.names=TRUE,recursive=FALSE,
                    include.dirs=FALSE,no..=TRUE)

species = list.dirs(modDir,recursive=FALSE)
species = species[-which(!is.na(str_match(species,"plots")))]
gridList = list()

for (i in 1:length(specs)){
  
  eval(parse(text=paste(specs[i],'PlotList = list()',sep="")))
  
  modDirSpec = which(str_detect(species,specs[i]))
  JDmodFiles = list.files(path=species[modDirSpec],pattern="*Model_TempLun3.Rdata",
                          full.names=TRUE,recursive=FALSE,include.dirs=FALSE,no..=TRUE)
  for (j in 1:length(sites[[i]])){
    
    JDInd = numeric()
    NTInd = numeric()
    MPhInd = numeric()
    NTMPhInd = numeric()
    
    # find master DF for this species & site
    whichDFSpec = which(str_detect(dfList,specs[i]))
    whichDFSite = which(str_detect(dfList,sites[[i]][j]))
    DFind = intersect(whichDFSpec,whichDFSite)
    
    # find model for this species & site
    whichModSite = which(str_detect(JDmodFiles,sites[[i]][j]))
    
    if (!isempty(DFind) & !isempty(whichModSite)){
      
      # load master DF & model
      thisSite = data.frame(read.csv(dfList[DFind]))
      load(JDmodFiles[whichModSite])
      
      JDInd = which(str_detect(names(tempMod$coefficients),"JDs")&!str_detect(names(tempMod$coefficients),"NTs"))
      NTInd = which(str_detect(names(tempMod$coefficients),"NTs")&!str_detect(names(tempMod$coefficients),"JDs")&!str_detect(names(tempMod$coefficients),"MPhs"))
      MPhInd = which(str_detect(names(tempMod$coefficients),"MPhs")&!str_detect(names(tempMod$coefficients),"NTs")&!str_detect(names(tempMod$coefficients),"LF"))
      NTMPhInd = which(str_detect(names(tempMod$coefficients),"MPhs")&str_detect(names(tempMod$coefficients),"NTs"))
      
      BootstrapParameters1<-rmvnorm(10000, coef(tempMod), summary(tempMod)$cov.unscaled)
      JDayBootstrapCoefs<- BootstrapParameters1[,JDInd]
      NTBootstrapCoefs<- BootstrapParameters1[,NTInd]
      MPhBootstrapCoefs<- BootstrapParameters1[,MPhInd]
      NTMPhBootstrapCoefs<- BootstrapParameters1[,NTMPhInd]
      
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
      MPhForPlotting<- seq(min(thisSite$MoonPhase), max(thisSite$MoonPhase), length=1000)
      MPhBasis<- mSpline(MPhForPlotting,
                         knots=quantile(thisSite$MoonPhase, probs=c(0.275,0.5,0.725)),
                         Boundary.knots=c(0,1),
                         periodic=T)
      
      # Calculate each smooth term at its mean value
      if (!isempty(JDInd)){JDmean = as.vector(apply(JDBasis,MARGIN=2,mean)%*%coef(tempMod)[JDInd])
      } else {JDmean = 0}
      if (!isempty(NTInd)){NTmean = as.vector(apply(NTBasis,MARGIN=2,mean)%*%coef(tempMod)[NTInd])
      } else {NTmean = 0}
      if (!isempty(MPhInd)) {MPhmean = as.vector(apply(MPhBasis,MARGIN=2,mean)%*%coef(tempMod)[MPhInd])
      } else {MPhmean=0}
      
      if (!isempty(MPhInd) & isempty(NTMPhInd)){ # if no interaction, plot moon phase on its own
        Fit<- MPhBasis%*%(coef(tempMod)[MPhInd]) # multiply basis functions by model coefficients to get values of spline at each X
        RealFit<- Fit+coef(tempMod)[1]+JDmean+NTmean # add intercept and other terms at their mean values
        BootstrapFits<- (MPhBasis%*%t(MPhBootstrapCoefs))+coef(tempMod)[1]+JDmean+NTmean # get spread of spline values at each X based on distributions of each coefficient
        quant.func<- function(x){quantile(x, probs=c(0.0275,0.975))}
        cis<-apply(BootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
        CIs = cbind(MPhForPlotting,t(cis))
        
        plotDF = data.frame(MPhForPlotting,RealFit)
        colnames(plotDF) = c("MoonPhase","Fit")
        
        CIDF = as.data.frame(CIs)
        colnames(CIDF) = c("MoonPhase","Ymin","Ymax")
        
        thisFig = ggplot(
        ) + geom_line(data=plotDF,aes(x=MoonPhase,y=Fit,),size=0.5,color="#AD1DC9"
        ) + geom_ribbon(data=CIDF,aes(x=MoonPhase,ymin=Ymin, ymax=Ymax),
                        fill="#AD1DC9",
                        alpha=0.2,
                        stat ="identity",
        ) + labs(x = NULL,
                 y = NULL
        ) + scale_y_continuous(breaks=c(round(min(CIDF$Ymin),digits=1),round(max(CIDF$Ymax),digits=1)),
                               labels=c(round(min(CIDF$Ymin),digits=1),round(max(CIDF$Ymax),digits=1))
        ) + coord_cartesian(ylim=c(round(min(CIDF$Ymin),digits=1),round(max(CIDF$Ymax),digits=1))                       
        ) + scale_x_continuous(breaks=c(0,0.5,1),
                               label=c("New","Full","New")
        ) + theme(axis.line = element_line(size=0.2),
                  panel.background = element_blank(),
                  plot.margin=margin(c(0,0,0,5),unit="pt")
        )
        
        eval(parse(text=paste(specs[i],'PlotList$',sites[[i]][j],'=thisFig',sep="")))
        
      } else if (!isempty(NTMPhInd)){ # if interaction, plot moon phase during daytime & nighttime
        
        FitMPh<- MPhBasis%*%(coef(tempMod)[MPhInd]) # multiply basis functions by model coefficients to get values of spline at each X
        FitNT<- NTBasis%*%(coef(tempMod)[NTInd]) # multiply basis functions by model coefficients to get values of spline at each X
        
        ## MPh at during daytime and nighttime
        IntPlotDF = numeric()
        IntCI = numeric()
        dayRows = which(NTForPlotting<0)
        nightRows = which(NTForPlotting>0)
        testRows = list(dayRows,nightRows)
        # testRows = list(1,500)
        
        for (l in 1:length(testRows)){
          IntBasis = numeric()
          # Recreate interaction basis
          for (m in 1:dim(MPhBasis)[1]){ # for each row in MPh basis
            rowMat = numeric()
            for (k in 1:3){ # multiply by all daytime OR nighttime values of NT
              rowMat = cbind(rowMat,NTBasis[unlist(testRows[l]),1]*MPhBasis[m,k],NTBasis[unlist(testRows[l]),2]*MPhBasis[m,k],NTBasis[unlist(testRows[l]),3]*MPhBasis[m,k])
            # IntBasis = cbind(IntBasis,NTBasis[unlist(testRows[l]),1]*MPhBasis[,k],NTBasis[unlist(testRows[l]),2]*MPhBasis[,k],NTBasis[unlist(testRows[l]),3]*MPhBasis[,k])
              }
            IntBasis = rbind(IntBasis,apply(rowMat,MARGIN=2,mean)) # average this MPh across NT
            
          }
          
          # Calculate MPh smooth in daytime/nighttime
          IntAdjust = IntBasis%*%coef(tempMod)[NTMPhInd]
          RealFitInt<- FitMPh+coef(tempMod)[1]+IntAdjust+JDmean # adjust offset
          IntBootstrapFits<- (MPhBasis%*%t(MPhBootstrapCoefs))+coef(tempMod)[1]+(IntBasis%*%t(NTMPhBootstrapCoefs))+JDmean # get spread of spline values at each X based on distributions of each coefficient
          quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
          cisInt<-apply(IntBootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
          IntCI = rbind(IntCI,cbind(MPhForPlotting,t(cisInt),rep(l,times=1000))) # get lower and upper CI bounds
          
          IntPlotDF = rbind(IntPlotDF,cbind(MPhForPlotting,RealFitInt,rep(l,times=1000)))
          
        }
        
        IntPlotDF = as.data.frame(IntPlotDF)
        colnames(IntPlotDF) = c("MoonPhase","Data","Fit")
        IntCI = as.data.frame(IntCI)
        colnames(IntCI) = c("MoonPhase","Ymin","Ymax","Fit2")
        IntPlotDF$Fit = as.factor(IntPlotDF$Fit)
        IntCI$Fit2 = as.factor(IntCI$Fit2)
        
        thisFig = ggplot(
        ) + geom_line(data=IntPlotDF,aes(x=MoonPhase,y=Data,group=Fit,color=Fit),size=0.5
        ) + scale_color_manual("Fit",values=c("#EC2934","#07B4CC"),
                               labels=c("Day","Night"),
                               name=NULL
        ) + geom_ribbon(data=IntCI,aes(x=MoonPhase,ymin=Ymin, ymax=Ymax,group=Fit2,fill=Fit2),
                        alpha=0.2,
                        stat ="identity"
        ) + scale_fill_manual("Fit2",values=c("#EC2934","#07B4CC"),
                              guide=NULL
        ) + labs(x = NULL,
                 y = NULL
        ) + scale_y_continuous(breaks=c(round(min(IntCI$Ymin),digits=1),round(max(IntCI$Ymax),digits=1)),
                               labels=c(round(min(IntCI$Ymin),digits=1),round(max(IntCI$Ymax),digits=1))
        ) + coord_cartesian(ylim=c(round(min(IntCI$Ymin),digits=1),round(max(IntCI$Ymax),digits=1))
        ) + scale_x_continuous(breaks=c(0,0.5,1),
                               label=c("New","Full","New")
        ) + theme(legend.position="none",
                  axis.line = element_line(size=0.2),
                  panel.background = element_blank(),
                  plot.margin=margin(c(0,0,0,5),unit="pt")
        )
        
        
        eval(parse(text=paste(specs[i],'PlotList$',sites[[i]][j],'=thisFig',sep="")))
        
      } else if (isempty(MPhInd) & isempty(NTMPhInd)){
        eval(parse(text=paste(specs[i],'PlotList$',sites[[i]][j],'=ggplot()',sep="")))
      }
      
    } else {eval(parse(text=paste(specs[i],'PlotList$',sites[[i]][j],' = ggplot()',sep="")))}
  }
  
  if (specs[i]=="Risso"){
    eval(parse(text=paste(specs[i],'PlotList$OC = ggplot()',sep="")))
    eval(parse(text=paste(specs[i],'PlotList$HAT = ggplot()',sep="")))
    eval(parse(text=paste(specs[i],'PlotList$GS = ggplot()',sep="")))
    eval(parse(text=paste(specs[i],'PlotList$BP = ggplot()',sep="")))
    eval(parse(text=paste(specs[i],'PlotList$BS = ggplot()',sep="")))
  }
  
  eval(parse(text=paste(specs[i],'Grid = grid.arrange(grobs=',specs[i],'PlotList,ncol=1,nrow=length(sites[[i]]),top=specTitles[i])',sep="")))
  
  gridList[[specs[i]]] = eval(parse(text=paste(specs[i],'Grid',sep="")))
  
}

png(file=paste(saveDir,'/DolphinMPhComp.png',sep=""),width=625,height=660)
grid.arrange(grobs=gridList,ncol=4,nrow=11,layout_matrix=rbind(c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(NA,NA,3,4),
                                                               c(NA,NA,3,4),
                                                               c(NA,NA,3,NA),
                                                               c(NA,NA,3,NA),
                                                               c(1,NA,3,NA)))
while (dev.cur()>1) {dev.off()}
pdf(file=paste(saveDir,'/DolphinMPhComp.pdf',sep=""),width=6.25,height=6.6)
grid.arrange(grobs=gridList,ncol=4,nrow=11,layout_matrix=rbind(c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(NA,NA,3,4),
                                                               c(NA,NA,3,4),
                                                               c(NA,NA,3,NA),
                                                               c(NA,NA,3,NA),
                                                               c(1,NA,3,NA)))
while (dev.cur()>1) {dev.off()}

