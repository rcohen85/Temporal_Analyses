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

specs = c("Sowerby","Cuvier","True","Sperm Whale","Blainville","Gervais","Kogia")
specTitles = c('Mb',"Zc","Mm","Pm","Md","Me","Kg")
sites = list(c("HZ","OC","NC","BC","WC"),
          c("HZ","OC","NC","BC","WC","NFC","HAT"),
          c("HZ","OC","NC","BC","WC","NFC"),
          c("HZ","OC","NC","BC","WC","NFC","HAT","GS","BP","BS"),
          c("GS","BP","BS"),
          c("GS","BP","BS"),
          c("GS","BP","BS"))

dfList = list.files(path=modelDFDir,pattern=paste('*',int,'_MasterTempLun.csv',sep=""),
                    full.names=TRUE,recursive=FALSE,
                    include.dirs=FALSE,no..=TRUE)

species = list.dirs(modDir,recursive=FALSE)
species = species[-which(!is.na(str_match(species,"plots")))]
gridList = list()

for (i in 1:length(specs)){

  modDirSpec = which(str_detect(species,specs[i]))
  JDmodFiles = list.files(path=species[modDirSpec],pattern="*Model_TempLun3.Rdata",
                          full.names=TRUE,recursive=FALSE,include.dirs=FALSE,no..=TRUE)
  
  if (specs[i]=="Sperm Whale"){
    # specs[i] = str_remove(specs[i]," ")
    eval(parse(text=paste(str_remove(specs[i]," "),'PlotList = list()',sep="")))
  } else {
    eval(parse(text=paste(specs[i],'PlotList = list()',sep="")))
  }
  
  
  for (j in 1:length(sites[[i]])){
    
    JDInd = numeric()
    NTInd = numeric()
    MPhInd = numeric()
    JDNTInd = numeric()
    
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
      JDNTInd = which(str_detect(names(tempMod$coefficients),"JDs")&str_detect(names(tempMod$coefficients),"NTs"))
      
      BootstrapParameters1<-rmvnorm(10000, coef(tempMod), summary(tempMod)$cov.unscaled)
      JDayBootstrapCoefs<- BootstrapParameters1[,JDInd]
      NTBootstrapCoefs<- BootstrapParameters1[,NTInd]
      MPhBootstrapCoefs<- BootstrapParameters1[,MPhInd]
      JDNTBootstrapCoefs<- BootstrapParameters1[,JDNTInd]
      
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
                         Boundary.knots=c(-1,1),
                         periodic=T)
      
      # Calculate each smooth term at its mean value
      if (!isempty(JDInd)){JDmean = as.vector(apply(JDBasis,MARGIN=2,mean)%*%coef(tempMod)[JDInd])
      } else {JDmean = 0}
      if (!isempty(NTInd)){NTmean = as.vector(apply(NTBasis,MARGIN=2,mean)%*%coef(tempMod)[NTInd])
      } else {NTmean = 0}
      if (!isempty(MPhInd)) {MPhmean = as.vector(apply(MPhBasis,MARGIN=2,mean)%*%coef(tempMod)[MPhInd])
      } else {MPhmean=0}
      
      if (!isempty(JDInd) & isempty(JDNTInd)){ # If JD is significant w no interaction
        Fit<- JDBasis%*%(coef(tempMod)[JDInd]) # multiply basis functions by model coefficients to get values of spline at each X
        RealFit<- Fit+coef(tempMod)[1]+NTmean+MPhmean # add intercept and other terms at their mean values
        BootstrapFits<- (JDBasis%*%t(JDayBootstrapCoefs))+coef(tempMod)[1]+NTmean+MPhmean # get spread of spline values at each X based on distributions of each coefficient
        quant.func<- function(x){quantile(x, probs=c(0.0275,0.975))}
        cis<-apply(BootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
        CIs = cbind(JDayForPlotting,t(cis))
        
        plotDF = data.frame(JDayForPlotting,RealFit)
        colnames(plotDF) = c("JDay","Fit")
        
        CIDF = as.data.frame(CIs)
        colnames(CIDF) = c("JDay","Ymin","Ymax")
        
        thisFig = ggplot(
        ) + geom_line(data=plotDF,aes(x=JDay,y=Fit,),size=0.5,color="#16A7CA"
        ) + geom_ribbon(data=CIDF,aes(x=JDay,ymin=Ymin, ymax=Ymax),
                        fill="#16A7CA",
                        alpha=0.2,
                        stat ="identity",
        ) + labs(x = NULL,
                 y = NULL
        ) + scale_y_continuous(breaks=c(round(min(CIDF$Ymin),digits=1),round(max(CIDF$Ymax),digits=1))
        ) + scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),
                               label=c("J","","","A","","","J","","","O","","")
        ) + coord_cartesian(ylim=c(round(min(CIDF$Ymin),digits=1),round(max(CIDF$Ymax),digits=1))
        ) + theme(axis.line = element_line(size=0.2),
                  panel.background = element_blank(),
                  plot.margin=margin(c(0,5,0,0),unit="pt")
        )
        
        if (specs[i]=="Sperm Whale"){
          eval(parse(text=paste(str_remove(specs[i]," "),'PlotList$',sites[[i]][j],'=thisFig',sep="")))
        } else {
          eval(parse(text=paste(specs[i],'PlotList$',sites[[i]][j],'=thisFig',sep="")))
        }
        
      } else if (!isempty(JDNTInd)){ # If there is an interaction btwn JD & NT
        
        FitJD<- JDBasis%*%(coef(tempMod)[JDInd]) # multiply basis functions by model coefficients to get values of spline at each X
        FitNT<- NTBasis%*%(coef(tempMod)[NTInd]) # multiply basis functions by model coefficients to get values of spline at each X
        
        ## Calculate JD smooth averaged across NT values
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
        
        IntAdjust = IntBasis%*%coef(tempMod)[JDNTInd]
        RealFitInt<- FitJD+coef(tempMod)[1]+IntAdjust+MPhmean # add intercept, interaction term, and other terms at their mean values
        IntBootstrapFits<- (JDBasis%*%t(JDayBootstrapCoefs))+coef(tempMod)[1]+(IntBasis%*%t(JDNTBootstrapCoefs))+MPhmean # get spread of spline values at each X based on distributions of each coefficient
        quant.func<- function(x){quantile(x, probs=c(0.025,0.975))}
        cisInt<-apply(IntBootstrapFits, 1, quant.func) # confidence interval of smooth function estimate
        IntCI = cbind(JDayForPlotting,t(cisInt)) # get lower and upper CI bounds
        
        IntPlotDF = cbind(JDayForPlotting,RealFitInt)
        
        IntPlotDF = as.data.frame(IntPlotDF)
        colnames(IntPlotDF) = c("JDay","Data")
        IntCI = as.data.frame(IntCI)
        colnames(IntCI) = c("JDay","Ymin","Ymax")
        
        thisFig = ggplot(
        ) + geom_line(data=IntPlotDF,aes(x=JDay,y=Data,),size=0.5,color="51D4A5"
        ) + geom_ribbon(data=IntCI,aes(x=JDay,ymin=Ymin, ymax=Ymax),
                        fill="51D4A5",
                        alpha=0.2,
                        stat ="identity",
        ) + labs(x = NULL,
                 y = NULL
        ) + coord_cartesian(xlim=c(1,365),
                            ylim=c(round(min(IntCI$Ymin),digits=1),round(max(IntCI$Ymax),digits=1))
        ) + scale_y_continuous(breaks=c(round(min(IntCI$Ymin),digits=1),round(max(IntCI$Ymax),digits=1))
        ) + scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),
                               label=c("J","","","A","","","J","","","O","","")
        ) + theme(legend.position="right",
                  axis.line = element_line(size=0.2),
                  panel.background = element_blank(),
                  plot.margin=margin(c(0,5,0,0),unit="pt")
        )
        
        if (specs[i]=="Sperm Whale"){
          eval(parse(text=paste(str_remove(specs[i]," "),'PlotList$',sites[[i]][j],'=thisFig',sep="")))
        } else {
        eval(parse(text=paste(specs[i],'PlotList$',sites[[i]][j],'=thisFig',sep="")))
        }
      }
      
    } else {eval(parse(text=paste(specs[i],'PlotList$',sites[[i]][j],' = ggplot()',sep="")))}
  }
  
  if (specs[i]=="Sperm Whale"){
    eval(parse(text=paste(str_remove(specs[i]," "),'Grid = grid.arrange(grobs=',str_remove(specs[i]," "),'PlotList,ncol=1,nrow=length(sites[[i]]),top=specTitles[i])',sep="")))
  } else {
  eval(parse(text=paste(specs[i],'Grid = grid.arrange(grobs=',specs[i],'PlotList,ncol=1,nrow=length(sites[[i]]),top=specTitles[i])',sep="")))
  }
  
  gridList[[str_remove(specs[i]," ")]] = eval(parse(text=paste(str_remove(specs[i]," "),'Grid',sep="")))
  
}

png(file=paste(saveDir,'/DeepDiverJDComp.png',sep=""),width=625,height=600)
grid.arrange(grobs=gridList,ncol=4,nrow=10,layout_matrix=rbind(c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(NA,2,3,4),
                                                               c(NA,2,NA,4),
                                                               c(5,6,7,4),
                                                               c(5,6,7,4),
                                                               c(5,6,7,4)))
while (dev.cur()>1) {dev.off()}
pdf(file=paste(saveDir,'/DeepDiverJDComp.pdf',sep=""),width=6.25,height=6)
grid.arrange(grobs=gridList,ncol=4,nrow=10,layout_matrix=rbind(c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(1,2,3,4),
                                                               c(NA,2,3,4),
                                                               c(NA,2,NA,4),
                                                               c(5,6,7,4),
                                                               c(5,6,7,4),
                                                               c(5,6,7,4)))
while (dev.cur()>1) {dev.off()}

