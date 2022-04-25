library(stringr)
library(performance)
source("binned_residuals_RC.R")

modDir = 'J:/Chpt_2/ModelsForMorgan2'
species = list.dirs(modDir,recursive=FALSE)

for (i in 1:length(species)){
  
  modFiles = list.files(path=species[i],pattern="*Model_TempLun3.Rdata",
                          full.names=TRUE,recursive=FALSE,include.dirs=FALSE,no..=TRUE)
  CTname = str_remove(species[i],paste(modDir,'/',sep=""))
  
  if (CTname=="Plots"){
    next
  } else {
  sites = list()
  for (j in 1:length(modFiles)){
    
    site = str_remove(str_remove(modFiles[j],paste(modDir,"/",CTname,"/",sep="")),"_5minBin_Model_TempLun3.Rdata")
    load(modFiles[j]) # load model
    
    png(file=paste(modDir,'/Plots/',CTname,'_',site,"_BinResid.png",sep=""),width = 400, height = 300, units = "px")
    print(binned_residuals_RC(tempMod))
    while (dev.cur()>1) {dev.off()}
  }
  }
}