clear all
close all
load Anders1060-6nmQw27In-Bar10nm_5Set19_300
%load Anders1060-6nmQw27In-Bar10nm_5Set19_350

Densityv=N;
lambdavet=lambda;

h=figure
set(h,'pos',[  420         340        1316         605])
subplot(131)

if length(size(G))==3
  plot(Densityv,squeeze(G(:,1,:)),'.-')
   if length(lambdavet)>length(Densityv)
    xlabel('Wavelenght nm')
   else
    xlabel('Carrier Dens')
   end
else   
  plot(Densityv,G)
  xlabel('Carrier Dens')

end
    ylabel('Gain 1/cm ')
if length(Tvet)>1    
legend(num2str(Tvet'),'location','best')
else
legend(num2str(lambdavet'),'location','best')
end

subplot(132)
plot(Densityv,EFspV)
ylabel('EF split (eV)')
    xlabel('Carrier Dens')

title(['Temp = ', num2str(Tvet)])

LaMat=ones(length(Densityv),1)*(.02*pi./(lambdavet*1e-9));    
gg=G./LaMat;  
nn=squeeze(Dep);

subplot(133)
plot(Densityv(1:end-1),diff(gg)./diff(nn))
ylim([-10 0])



ylabel('Henry factor')
xlabel('Carrier Dens')  
