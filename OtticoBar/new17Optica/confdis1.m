load old
lami=lambda*1000;
h=figure, 
%set(h,'pos',[841   102   434   811])
set(h,'pos',[  358   177   667   483])
subplot(211), plot(lami+fou*fafr,aou,lami+zev*fafr,zeros(size(zev)),'wo') ,
title([' Dlam = ',num2str(zev*fafr),'          Lam= ',num2str(zev*fafr+lambda*1000)])

grid
subplot(212), semilogy(lami+fou*fafr,gou,lami+zev*fafr,g0v,'wo')
grid

  load new
  figure(h),
subplot(211), hold on, plot(lami+fou*fafr,aou,'.--',lami+zev*fafr,zeros(size(zev)),'wo') ,

subplot(212), hold on, semilogy(lami+fou*fafr,gou,'.--',lami+zev*fafr,g0v,'wo')  
 
  load new0
  figure(h),
subplot(211), hold on, plot(lami+fou*fafr,aou,'s',lami+zev*fafr,zeros(size(zev)),'wo') ,

subplot(212), hold on, semilogy(lami+fou*fafr,gou,'s',lami+zev*fafr,g0v,'wo')   

  load speed0
  figure(h),
subplot(211), hold on, plot(lami+fou*fafr,aou,'+',lami+zev*fafr,zeros(size(zev)),'wo') ,

subplot(212), hold on, semilogy(lami+fou*fafr,gou,'+',lami+zev*fafr,g0v,'wo')  
  
