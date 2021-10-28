IntegralTimeScaleCalc = function(ts){
  
  avt = mean(ts)
  dcc = ts-avt
  
  [y,lag] = xcorr(dcc,'coeff');
  

  plot(lag,y)
   
  i = ceiling(length(y)/2);
  its = 0;
  
  while (y[i]>0){
  its = its+y[i];
  i = i+1;
  }

  return(its)
  
}