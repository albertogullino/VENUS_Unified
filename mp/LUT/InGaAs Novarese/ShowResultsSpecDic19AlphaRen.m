clear all
close all




TV=[300 325 350];

%TV=[350];

for T=TV

eval(['load Anders1060-6nmQw318In-Bar10nm_17Dic19SpRen',num2str(T)]')


 
Densityv=N;
lambdavet=lambda;

h=figure
set(h,'pos',[  420         340        1316         605])
subplot(121)


  plot(lambdavet,G)
  xlabel('Wavelength (nm)')
  grid


    ylabel('Gain 1/cm ')
if length(Tvet)>1    
legend(num2str(Tvet'),'location','best')
else
legend(num2str(N'*1e-12),'location','best')
end





LaMat=ones(length(Densityv),1)*(.02*pi./(lambdavet*1e-9));    
gg=(G-Gi)./LaMat;  
nn=squeeze(Dep);
nn=squeeze(Dep-Depi);

al=nn./gg;


subplot(122)
plot(lambdavet,al)
ylim([-2 .2])
grid



%ylabel('(1+\alpha^2)')
ylabel('\alpha')
xlabel('Wavelength (nm)')  
title(['Temp = ', num2str(Tvet)])
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