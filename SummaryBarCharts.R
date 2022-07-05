library(ggplot2)
library(ggpattern)
library(stringr)
library(gridExtra)

outDir = 'J:/Chpt_2/Figures'
totPres = data.frame(read.csv('J:/Chpt_2/PresenceQuantification.csv'))
load('J:/Chpt_2/ModelOutput/PresStats.Rdata')
names(PresStats)[7] = 'SpermWhale'
sites = c("HZ","OC","NC","BC","WC","NFC","HAT","GS","BP","BS","JAX")


# Plot bars of presence colored by season of peak ---------------------------
my_colors=rev(c("#2E22EA","#DA3DFB","#FA4855","#FCDA7B","#D8EC28","#4BBA0F","#2EB991","#2C24E9"))


hues_df = data.frame(Jday = 1:365,
    label="",
    colors = colorRampPalette(my_colors)(365))
hues_df$label[which(hues_df$Jday%in%c(1,92,183,274))] = c("Win","Spr","Sum","Fall")

colorWheel = ggplot(hues_df
)+geom_rect(aes(ymin=2,
                ymax=4,
                xmin=Jday-0.5,
                xmax=Jday+0.5,
                color=NULL,
                fill=colors),
            alpha=0.6
)+coord_polar(direction=1,
              start=0
)+scale_color_identity(
)+scale_fill_identity(
)+guides(fill='none',color='none'
)+theme(legend.position="none"
)+theme_void(
)+ylim(c(1,6)
)+geom_text(aes(x=Jday,
                y=5.5,
                label=label),
            size=4) 


countInd = which(str_detect(colnames(totPres),"Prop"))
plotList = list()
specName = c("Me","Md","Zc","Kg","Gg1","Mb","Pm","Mm","Gm","Dd","Gg2")
nonSignif = "#CCCCCC"

for (i in 1:length(countInd)){
  spec = str_remove(colnames(totPres)[countInd[i]],"Prop")
  
  plotDF = data.frame(sites=sites,pres=totPres[,countInd[i]]*1000)
  plotDF$pres[plotDF$pres<1] = 1
  plotDF$pres = log10(plotDF$pres)
  plotDF$Pattern = 'none'
  
  # get line/fill colors for each site according to season of peak presence
  plotDF$LinCol = hues_df$colors[round(as.numeric(PresStats[[spec]]$JDPeak))]
  plotDF$Fill = hues_df$colors[round(as.numeric(PresStats[[spec]]$JDPeak))]
  # set line/fill color to grey for sites where JD wasn't significant
  badInd = c(which(as.numeric(PresStats[[spec]]$JDSignif)>=0.05),which(is.na(PresStats[[spec]]$JDSignif)))
  plotDF$LinCol[badInd] = nonSignif
  plotDF$Fill[badInd] = nonSignif
  # set line/fill color to black/white w stripes for sites where model performance was poor
  badInd = which(PresStats[[spec]]$ModPerf<0.6)
  plotDF$LinCol[badInd] = "#CCCCCC"
  plotDF$Fill[badInd] = "#FFFFFF"
  plotDF$Pattern[badInd] = 'stripe'
  # set line/fill color to black/white for sites w no model
  badInd = which(is.na(PresStats[[spec]]$ModPerf))
  plotDF$LinCol[badInd] = "#CCCCCC"
  plotDF$Fill[badInd] = "#FFFFFF"
  
  # if (spec=="SpermWhale"){ # catch for non-converged model at BC
  #   plotDF$LinCol[4] = "#CCCCCC"
  #   plotDF$Fill[4] = "#FFFFFF"
  #   plotDF$Pattern[4] = 'stripe'
  # }
  if (spec=="Gervais"){ # catch for non-modeled site (HAT) for Me
    plotDF$LinCol[7] = "#CCCCCC"
    plotDF$Fill[7] = "#FFFFFF"
    plotDF$Pattern[7] = 'stripe'
  }
  
  # if (spec=="UD28"){ylabs=rev(sites)}else{ylabs=NULL}
  ylabs=NULL
  
  plotList[[spec]] = ggplot(plotDF,aes(y=sites,x=pres,fill=sites,color=sites)
  )+geom_col_pattern(alpha=0.6,
                     pattern_colour="#CCCCCC",
                     pattern_fill="#FFFFFF",
                     pattern_angle=45,
                     pattern=plotDF$Pattern
  )+scale_fill_manual(breaks=sites,values=c(plotDF$Fill),guide=NULL
  )+scale_color_manual(breaks=sites,values=plotDF$LinCol,guide=NULL
  )+scale_y_discrete(limits=rev(sites),labels=ylabs
  )+coord_cartesian(xlim=c(0,2.5)
  )+scale_x_continuous(breaks=c(0,1,2),
                       labels=c(0.1,1,10)
                       # labels=c(expression("10"^0),expression("10"^1),expression("10"^2))
  )+annotation_logticks(base=10,
                        sides="b",
                        scaled=TRUE,
                        alpha=0.4,
                        short=unit(0,"points"),
                        mid=unit(0,"points"),
                        long=unit(2,"mm")
  )+labs(x=NULL,y=NULL,
         title=specName[i]
  )+theme_minimal(
  )+theme(legend.position="none",
            plot.title=element_text(size=16),
          panel.grid.major=element_line(color="#CCCCCC"),
          panel.grid.minor.x=element_blank(),
          # axis.text.x=element_text(size=14,angle=315,vjust=0),
          axis.text.x=element_text(size=14),
          axis.text.y=element_text(size=14))
}

plotList =plotList[c("UD28","Risso","UD36","UD26","Blainville","Gervais","Cuvier","Sowerby","True","Kogia","SpermWhale")]
plotList$CW = colorWheel

png(file=paste(outDir,'/SeasonalPres.png',sep=""),width = 1600, height = 450, units = "px",res=125)
grid.arrange(grobs=plotList,ncol=12,nrow=2,layout_matrix=rbind(c(seq(11),NA),c(seq(12))),
             bottom="Bins per Thousand with Presence",left="Site")
while (dev.cur()>1) {dev.off()}
pdf(file=paste(outDir,'/SeasonalPres.pdf',sep=""),width = 7, height = 2)
grid.arrange(grobs=plotList,ncol=12,nrow=2,layout_matrix=rbind(c(seq(11),NA),c(seq(12))),
             bottom="Bins per Thousand with Presence",left="Site")
while (dev.cur()>1) {dev.off()}


# Plot bars of presence colored by diel prop -------------------------
my_colors=c("#0F52EE","#FFFFFF","#EEED0F")

hues_df = data.frame(Prop = (0:100)/100,
                     label="",
                     colors = colorRampPalette(my_colors)(101))
hues_df$label[which(hues_df$Prop%in%c(0,1))] = c(" Nocturnal","Diurnal")

colorBar = ggplot(hues_df
)+geom_rect(aes(ymin=2.5,
                ymax=3.5,
                xmin=Prop-0.01,
                xmax=Prop+0.01,
                color=NULL,
                fill=colors),
            alpha=0.6
)+scale_color_identity(
)+scale_fill_identity(
)+guides(fill='none',color='none'
)+theme(legend.position="none"
)+theme_void(
)+ylim(c(1,5)
# )+geom_text(aes(x=Prop,
#                 y=4,
#                 label=label),
#             size=4
)+coord_flip() 


countInd = which(str_detect(colnames(totPres),"Prop"))
plotList = list()
specName = c("Me","Md","Zc","Kg","Gg1","Mb","Pm","Mm","Gm","Dd","Gg2")
nonSignif = "#CCCCCC"

for (i in 1:length(countInd)){
  spec = str_remove(colnames(totPres)[countInd[i]],"Prop")
  
  plotDF = data.frame(sites=sites,pres=totPres[,countInd[i]]*1000)
  plotDF$pres[plotDF$pres<1] = 1
  plotDF$pres = log10(plotDF$pres)
  plotDF$LinCol = "#CCCCCC"
  plotDF$Fill = NA
  plotDF$Pattern = 'none'
  # set fill color based on value of DielProp
  whichCols = match(round(as.numeric(PresStats[[spec]]$DielProp),2),hues_df$Prop)
  whichSites = which(!is.na(PresStats[[spec]]$DielProp))
  # plotDF$LinCol[whichSites] = hues_df$colors[whichCols[!is.na(whichCols)]]
  plotDF$Fill[whichSites] = hues_df$colors[whichCols[!is.na(whichCols)]]
  # set fill color to grey for sites where NT wasn't significant, either on its own or in an interaction
  goodInd = c(which((as.numeric(PresStats[[spec]]$NTSignif)<=0.05|as.numeric(PresStats[[spec]]$NTJDSignif)<=0.05|as.numeric(PresStats[[spec]]$NTMPhSignif)<=0.05)))
  badInd = setdiff(seq(11),goodInd)
  # plotDF$LinCol[badInd] = nonSignif
  plotDF$Fill[badInd] = nonSignif
  # set fill color to white w stripes for sites where model performance was poor
  badInd = which(PresStats[[spec]]$ModPerf<0.6)
  # plotDF$LinCol[badInd] = "#CCCCCC"
  plotDF$Fill[badInd] = "#FFFFFF"
  plotDF$Pattern[badInd] = 'stripe'
  # set fill color to white w crosshatch for sites w no model
  badInd = which(is.na(PresStats[[spec]]$ModPerf))
  # plotDF$LinCol[badInd] = "#CCCCCC"
  plotDF$Fill[badInd] = "#FFFFFF"
  plotDF$Pattern[badInd] = 'crosshatch'
  
  # if (spec=="SpermWhale"){ # catch for non-converged model at BC
  #   # plotDF$LinCol[4] = "#CCCCCC"
  #   plotDF$Fill[4] = "#FFFFFF"
  #   plotDF$Pattern[4] = 'stripe'
  # }
  if (spec=="Gervais"){ # catch for non-modeled site (HAT) for Me
    # plotDF$LinCol[7] = "#CCCCCC"
    plotDF$Fill[7] = "#FFFFFF"
    # plotDF$Pattern[7] = 'crosshatch'
    plotDF$Pattern[7] = 'stripe'
    
  }
  
  # if (spec=="UD28"){ylabs=rev(sites)}else{ylabs=NULL}
  ylabs=NULL
  
  plotList[[spec]] = ggplot(plotDF,aes(y=sites,x=pres,fill=sites,color=sites)
  )+geom_col_pattern(pattern_colour="#CCCCCC",
                     pattern_fill="#FFFFFF",
                     pattern_angle=45,
                     pattern=plotDF$Pattern
  )+scale_fill_manual(breaks=sites,values=c(plotDF$Fill),guide=NULL
  # )+scale_color_manual(breaks=sites,values=plotDF$LinCol,guide=NULL
  )+scale_color_manual(breaks=sites,values=c(rep("#CCCCCC",11)),guide=NULL
  )+scale_y_discrete(limits=rev(sites),labels=ylabs
  )+coord_cartesian(xlim=c(0,2.5)
  )+scale_x_continuous(breaks=c(0,1,2),
                       labels=c(1,10,100)
                       # labels=c(expression("10"^0),expression("10"^1),expression("10"^2))
  )+annotation_logticks(sides="b",
                        scaled=TRUE,
                        alpha=0.4,
                        short=unit(0,"points"),
                        mid=unit(0,"points"),
                        long=unit(2,"mm")
  )+labs(x=NULL,y=NULL,
         title=specName[i]
  )+theme_minimal(
  )+theme(legend.position="none",
          plot.title=element_text(size=16),
          panel.grid.major=element_line(color="#CCCCCC"),
          panel.grid.minor.x=element_blank(),
          # axis.text.x=element_text(size=14,angle=315,vjust=0),
          axis.text.x=element_text(size=14),
          axis.text.y=element_text(size=14))
}

plotList =plotList[c("UD28","Risso","UD36","UD26","Blainville","Gervais","Cuvier","Sowerby","True","Kogia","SpermWhale")]
plotList$CW = colorBar

png(file=paste(outDir,'/DielPres.png',sep=""),width = 1600, height = 450, units = "px",res=125)
grid.arrange(grobs=plotList,ncol=12,nrow=4,layout_matrix=rbind(c(seq(11),NA),c(seq(12)),c(seq(12)),c(seq(11),NA)),
             bottom="Bins per Thousand with Presence",left="Site")
while (dev.cur()>1) {dev.off()}
pdf(file=paste(outDir,'/DielPres.pdf',sep=""),width = 7, height = 2)
grid.arrange(grobs=plotList,ncol=12,nrow=4,layout_matrix=rbind(c(seq(11),NA),c(seq(12)),c(seq(12)),c(seq(11),NA)),
             bottom="Bins per Thousand with Presence",left="Site")
while (dev.cur()>1) {dev.off()}
