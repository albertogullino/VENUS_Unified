clear all
close all

barrier=input(' 4 [0] or 6 [1] nm barrier? ')
if barrier==0
 load LUT5QW-Bar4nm_Jun19
else
 load LUT5QW-Bar6nm_Jun19
end

 h=figure;
 set(h,'pos',[120         415        1796         525])
 subplot(131)
 plot(lambda,squeeze(Es(:,:,1)))
 xlabel('Wavelenght nm')
 ylabel('Es 1/s')

 subplot(133)
 plot(N,Rsp)
  xlabel('Carrier Density') 
  ylabel('R_{sp} 1/(cm^3 s)')
  
 subplot(132)
  plot(lambda,squeeze(G(:,:,1)))
   xlabel('Wavelenght nm')

    ylabel('Gain 1/cm ')
  
