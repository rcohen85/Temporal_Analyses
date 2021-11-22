IntegralTimeScaleCalc = function(ts){
  
  avt = mean(ts,na.rm=TRUE)
  dcc = ts-avt
  
  q = acf(dcc,na.action=na.pass,plot=FALSE);
   
  i = ceiling(length(q$acf)/2);
  its = 0;
  
  while (q$acf[i]>0){
  its = its+q$acf[i];
  i = i+1;
  }

  return(its)
  
}