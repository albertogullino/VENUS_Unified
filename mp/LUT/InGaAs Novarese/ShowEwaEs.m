clear all
close all





TV=[300 325 350];

for T=TV

eval(['load Anders1060-6nmQw318In-Bar10nm_17Dic19renF',num2str(T)]')

Densityv=N;
lambdavet=lambda;

h=figure
set(h,'pos',[  206         341        1631         605])
subplot(141)

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

subplot(142)
plot(Densityv,EFspV)
ylabel('EF split (eV)')
    xlabel('Carrier Dens')
title(['Temp = ', num2str(Tvet)])

hcel=6.62e-34*3e17/1.6e-19;    

hcl=ones(length(Densityv),1)*(hcel./(lambdavet));    

DE=repmat(EFspV,1,length(lambdavet));
DEL=(hcl-DE)/.025;

subplot(143)
nsp=1./(1-exp(DEL));
plot(Densityv,nsp.*G)
ylabel('EF split (eV)')
    xlabel('Carrier Dens')

fe=1e-12;
figure, plot(Densityv*fe,G.*nsp*0.8e10)
figure, plot(Densityv*fe,Es)
figure, plot(Densityv*fe,Es,'linewidth',2)
chold
plot(Densityv*fe,G.*nsp*0.8e10,'--')
ylim([0 8e13])
grid
xlabel(' N  (1e12/cm^2)')
ylabel('R_{sp} (1/s)')
pausak
figure, plot(Densityv*fe,Es./G,'linewidth',2)
xlabel(' N  (1e12/cm^2)')
ylabel('Effective n_{sp}')
ylabel('Effective n_{sp} x v_g')




LaMat=ones(length(Densityv),1)*(.02*pi./(lambdavet*1e-9));    
gg=G./LaMat;  
nn=squeeze(Dep);

subplot(144)
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