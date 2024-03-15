clear all
close all

Fi='Anders1060-6nmQw318In-Bar10nm_17Dic19T'


TV=[300 325 350];

for T=TV

eval(['load Anders1060-6nmQw318In-Bar10nm_17Dic19T',num2str(T)]')


 

 h=figure;
 set(h,'pos',[120         415        1196         525])
 %figure(h)
 subplot(121)
 En=squeeze(Es(:,:,1));
 Ma=max(En,[],2);
 En=diag(1./Ma)*En;
 %hold on
 plot(lambda,En)
grid
 xlabel('Wavelenght nm')
 ylabel(' Normalized Es')
 legend(num2str(1e-12*N'),'location','best')

% subplot(133)
% plot(N,Rsp)
%  xlabel('Carrier Density') 
%  ylabel('R_{sp} 1/(cm^3 s)')

 subplot(122)
 %hold on 
  plot(lambda,squeeze(G(:,:,1)))
   xlabel('Wavelenght nm')
grid
    ylabel('Gain 1/cm ')
title(['Temp =', num2str(Tvet)])

pausak

end


return

load Anders1060-6nmQw318In-Bar10nm_17Dic19T350

 h=figure;
 set(h,'pos',[120         415        1196         525])
 %figure(h)
 subplot(121)
 En=squeeze(Es(:,:,1));
 Ma=max(En,[],2);
 En=diag(1./Ma)*En;
 grid
 %hold on
 plot(lambda,En)
 xlabel('Wavelenght nm')
 ylabel(' Normalized Es')
legend(num2str(1e-12*N'),'location','best')
grid
% subplot(133)
% plot(N,Rsp)
%  xlabel('Carrier Density') 
%  ylabel('R_{sp} 1/(cm^3 s)')

 subplot(122)
 %hold on 
  plot(lambda,squeeze(G(:,:,1)))
   xlabel('Wavelenght nm')
grid
    ylabel('Gain 1/cm ')
title(['Temp =', num2str(Tvet)])