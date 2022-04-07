---
  title: "Figures call characteristics"
author: "Annebelle Kok"
date: "13 januari 2019"
output: html_document
---
  # NB: save figures with dimensions: 1200 x 868!
  # #3D plot of duration, depth and frequency
  # ```{r setup, include=FALSE}
setwd("~/Promotie/Call propagation and noise/Resultaten")
res2<-read.delim("PSDdata_peakf_15022019_nonegSNR.txt")
resSNR10<-subset(res2,res2$snr>=10)
resSNR10$SNRrec<-resSNR10$TOLRecHz-resSNR10$SPLtotRecHz
summary(resSNR10$SNRrec)
resSNR1010<-resSNR10[-86,] #there is one call with SNRrec<0 that has to be removed. 
summary(resSNR1010$SNRrec)
resfinal<-resSNR1010
resfinal$Dcorr.uncorr<-resfinal$SPLcorrHz-resfinal$SPLuncorrHz
resfinal.flow<-subset(resfinal,resfinal$Dcorr.uncorr>=6)
resfinal.138a<-subset(resfinal,resfinal$Tag=="gm138a")
resfinal.high<-subset(resfinal,resfinal$Frequency>6400)
resfinal.flow2<-rbind(resfinal.flow,resfinal.138a,resfinal.high)
resfinal.flow.na<-resfinal.flow2[,c(1,2,3,4,5,6,12,14,19,20,25,26,30,33)]
subset(resfinal.flow.na,resfinal.flow.na$Depth<0)
resfinal.flow.na2<-resfinal.flow.na[complete.cases(resfinal.flow.na),]
require(ggplot2)
lowF<-subset(resfinal.flow.na2,resfinal.flow.na2$Frequency<2500)
mediumF<-subset(resfinal.flow.na2,resfinal.flow.na2$Frequency>2500&resfinal.flow.na2$Frequency<7500)
#mediumF2<-subset(resfinal.flow.na2,resfinal.flow.na2$Frequency>5000&resfinal.flow.na2$Frequency<7500)
highF<-subset(resfinal.flow.na2,resfinal.flow.na2$Frequency>=7500)
intercF<-subset(resfinal.flow.na2,resfinal.flow.na2$Frequency==1000)
intercFD<-subset(intercF,intercF$Depth<0.05)
# ```
#Theme for plot lay-out
# ```{r}
mijnthema<-theme(panel.spacing=unit(2,"lines"),legend.position="none",plot.title=element_text(size=32, face="bold",hjust=0.5,vjust=1),axis.text.x = element_text(size=26,colour="black"), axis.text.y=element_text(size=26,colour="black"),axis.title.x=element_text(size=30,colour="black"), axis.title.y=element_text(size=30,colour="black"),  axis.line = element_line(colour="black",size=1), panel.background=element_rect(fill= "transparent",colour =NA),legend.background= element_rect(fill="transparent",colour=NA),axis.line.x = element_line(colour="black"),axis.line.y=element_line(colour="black"))
# ```
#Fig 5a: 3D plot of duration, frequency and SPL at producer (TOL)
# ```{r}
g<-ggplot(data=resfinal.flow.na2,aes(x=TOLHz,y=Duration))
g+geom_point(data=lowF,fill="black",shape=21,size=4)+
  geom_point(data=mediumF,fill="darkgrey",shape=21,size=4)+
  geom_point(data=highF,fill="red",shape=21,size=4)+
  stat_smooth(data=lowF,method=lm,colour="black",fullrange=F,linetype=2)+
  stat_smooth(data=mediumF,method=lm,colour="darkgrey",fullrange=F)+
  stat_smooth(data=highF,method=lm,colour="red",fullrange=F,linetype=3)+
  annotation_logticks(base = 10,sides="l",scaled = TRUE)+
  mijnthema+
  labs(x=expression("Call PSD at producer (dB re 1 ?"*Pa^2/"Hz)"),y=expression("Duration (s)"))+
  scale_fill_manual(values=c("#FFFFFF","#000000"))+
  #scale_x_log10()+
  scale_y_log10()
# ```
# Fig 5b: depth vs duration and peak frequency
# ```{r}
g<-ggplot(data=resfinal.flow.na2,aes(x=Depth,y=Duration))
g+geom_point(data=lowF,fill="black",shape=21,size=4)+
  geom_point(data=mediumF,fill="darkgrey",shape=21,size=4)+
  geom_point(data=highF,fill="red",shape=21,size=4)+
  stat_smooth(data=lowF,method=lm,colour="black",fullrange=F,linetype=2)+
  stat_smooth(data=mediumF,method=lm,colour="darkgrey",fullrange=F)+
  stat_smooth(data=highF,method=lm,colour="red",fullrange=F,linetype=3)+
  annotation_logticks(base = 10,sides="bl",scaled = TRUE)+
  mijnthema+
  labs(x=expression("Depth of production (m)"),y=expression("Duration (s)"))+
  scale_x_log10(breaks = c(0.1,10,500),labels = c(0.1,10,500))+
  scale_y_log10()
# ```
#Fig 5c: 3D plot of SPL at producer (TOL), depth and frequency
# ```{r}
g<-ggplot(data=resfinal.flow.na2,aes(x=Depth,y=TOLHz))
g+geom_point(data=lowF,fill="black",shape=21,size=4)+
  geom_point(data=mediumF,fill="darkgrey",shape=21,size=4)+
  geom_point(data=highF,fill="red",shape=21,size=4)+
  stat_smooth(data=lowF,method=lm,colour="black",fullrange=F,linetype=2)+
  stat_smooth(data=mediumF,method=lm,colour="darkgrey",fullrange=F)+
  stat_smooth(data=highF,method=lm,colour="red",fullrange=F,linetype=3)+
  annotation_logticks(base = 10,sides="b",scaled = TRUE)+
  mijnthema+
  labs(x=expression("Depth of production (m)"),y=expression("Call PSD at producer (dB re 1 ?"*Pa^2/"Hz)"))+
  scale_fill_manual(values=c("#FFFFFF","#000000"))+
  #scale_x_log10()+
  scale_x_log10(breaks = c(0.1,10,500),labels = c(0.1,10,500))
# ```
#Fig 5d: plot of ambient noise at producer vs call amplitude and peak frequency
# ```{r}
g<-ggplot(data=resfinal.flow.na2,aes(x=SPLtotHz,y=TOLHz))
g+geom_point(data=lowF,fill="black",shape=21,size=4)+
  geom_point(data=mediumF,fill="darkgrey",shape=21,size=4)+
  geom_point(data=highF,fill="red",shape=21,size=4)+
  stat_smooth(method=lm,fullrange=F)+
  mijnthema+
  labs(x=expression("Ambient noise PSD at producer (dB re 1 ?"*Pa^2/"Hz)"),y=expression("Call PSD at producer (dB re 1 ?"*Pa^2/"Hz)"))
# ```
#3D plot of duration, depth, and SPL at producer (TOL)
# ```{r}
#grouping SPLsig into categories
#lowtol<-subset(resfinal.na3,resfinal.na3$TOLHz<90)
#midtol<-subset(resfinal.na2,resfinal.na2$TOL>=90&resfinal.na2$TOL<110)
#hightol<-subset(resfinal.na2,resfinal.na2$TOL>=110)
lowdur<-subset(resfinal.flow.na2,resfinal.flow.na2$Duration<0.2)
middur<-subset(resfinal.flow.na2,resfinal.flow.na2$Duration>=0.2&resfinal.flow.na2$Duration<1)
highdur<-subset(resfinal.flow.na2,resfinal.flow.na2$Duration>=1)
g<-ggplot(data=resfinal.flow.na2,aes(x=TOLHz,y=Depth))
g+geom_point(shape=21,size=2)+
  stat_smooth(data=lowdur,method=lm,colour="black",fullrange=F)+
  stat_smooth(data=middur,method=lm,colour="grey",fullrange=F)+
  stat_smooth(data=highdur,method=lm,colour="red",fullrange=F)+
  annotation_logticks(base = 10,sides="l",scaled = TRUE)+
  theme(panel.spacing=unit(2,"lines"),legend.position="none",axis.text.x = element_text(size=26), axis.text.y=element_text(size=26),axis.title.x=element_text(size=30), axis.title.y=element_text(size=30), panel.background=element_rect(fill= "transparent",colour =NA),legend.background= element_rect(fill="transparent",colour=NA),axis.line.x = element_line(colour="black"),axis.line.y=element_line(colour="black"))+
  xlab("\nCall SPSD at producer (dB re 1 ?Pa2/Hz)")+
  ylab("Depth of production (m)\n")+
  scale_fill_manual(values=c("#FFFFFF","#000000"))+
  #scale_x_log10()+
  scale_y_log10(breaks = c(0.1,10,500),labels = c(0.1,10,500))
# ```
# ```{r}
longcall<-subset(highSNR.na,highSNR.na$Duration>=1)
shortcall<-subset(highSNR.na,highSNR.na$Duration<=1)
highSNR.na2<-highSNR.na[c(1:352,354:421),]
g<-ggplot(data=highSNR.na2,aes(x=Duration,y=DepthP))
g+geom_point(aes(fill=SPLsigcorr),shape=21,size=2)+
  stat_smooth(method=lm,colour="black",fullrange=F)+
  theme(panel.spacing=unit(2,"lines"),legend.position="none",axis.text.x = element_text(size=26), axis.text.y=element_text(size=26),axis.title.x=element_text(size=30), axis.title.y=element_text(size=30), panel.background=element_rect(fill= "transparent",colour =NA),legend.background= element_rect(fill="transparent",colour=NA),axis.line.x = element_line(colour="black"),axis.line.y=element_line(colour="black"))+
  xlab("\nDuration (s)")+
  ylab("Depth of production (m)\n")+
  scale_y_reverse()
# ```
# ```{r}
g<-ggplot(data=resfinal.na2,aes(x=SPLtot,y=Frequency))
g+geom_point()+
  stat_smooth(method=lm, fullrange=F,colour="black")+
  theme(panel.spacing=unit(2,"lines"),legend.position="none",axis.text.x = element_text(size=26), axis.text.y=element_text(size=26),axis.title.x=element_text(size=30), axis.title.y=element_text(size=30), panel.background=element_rect(fill= "transparent",colour =NA),legend.background= element_rect(fill="transparent",colour=NA),axis.line.x = element_line(colour="black"),axis.line.y=element_line(colour="black"))+
  xlab("\nNoise level at producer (dB re 1 ?Pa)")+
  ylab("SPL of call at producer (dB re 1 ?Pa)\n")
# ```
#Figure S1: Signal excess at receiver vs peak frequency at producer.
# ```{r}
# Add critical noise levels to sound pressure spectral density of ambient noise
CR<-c(18,18,18,18,18,18,20,23,23,22,26,26,28,28,30,32)
resfinal.flow2$DT<-resfinal.flow2$SPLtotRecHz
a<-levels(as.factor(resfinal.flow2$Frequency))
Freqfac<-as.factor(resfinal.flow2$Frequency)
resfinal.flow2<-cbind(resfinal.flow2,Freqfac)
for(n in 1:length(a)){
  b<-which(resfinal.flow2$Freqfac==a[n])
  resfinal.flow2$DT[b]<-resfinal.flow2$SPLtotRecHz[b]+CR[n]
}
resfinal.flow2$snrRecDT<-resfinal.flow2$TOLRec-resfinal.flow2$DT
g<-ggplot(data=resfinal.flow2,aes(x=Frequency,y=snrRecDT))
g+geom_point(size=4)+
  mijnthema+
  xlab("\npeak frequency of call at producer (Hz)")+
  ylab("Signal excess at receiver (dB)\n")
# ```