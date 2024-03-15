clear all
close all





TV=[300 325 350];

for T=TV

eval(['load Anders1060-6nmQw318In-Bar10nm_17Dic19Alpha',num2str(T)]')

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
plot(Densityv(1:end-1),(diff(nn)./diff(gg)).^1)
%ylim([ 0 5])



%ylabel('(1+\alpha^2)')
ylabel('\alpha')
xlabel('Carrier Dens')  

pausak

end

return

load Anders1060-6nmQw318In-Bar10nm_17Dic19Alpha350
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
plot(Densityv(1:end-1),(diff(nn)./diff(gg)).^1)
%ylim([0 5])


%ylabel('(1+\alpha^2)')
ylabel('\alpha')
xlabel('Carrier Dens')  