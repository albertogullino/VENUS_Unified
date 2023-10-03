clear
close all
clc

addpath E:\Dati\mvcsel\new10Optica
addpath E:\Dati\mvcsel\new10Optica\AlbertoBW
lambda=3300;
DC=0.5;
nin=4.67;
nout=1;
n1=3.43;
n2=1;
thickness=200;
NModes=11;
DC=.5;

vPeriod=logspace(0,3,101)*.7;
vPeriod=[500 700 900];
vTh=[20 40 80 160 320];
vTh=[380];
vDC=linspace(.2,.8,61);

%vPeriod=250;
%for indDC=1:length(vDC)
%  DC=vDC(indDC);
 thickness=vTh
%for indDC=1:length(vDC)
%  DC=vDC(indDC);
 for indPeriod=1:length(vPeriod)
    Period=vPeriod(indPeriod);
for indDC=1:length(vDC)
    DC=vDC(indDC);

    [neff,nBW,flagt] = f_EffectiveIndex(Period,lambda,DC,nin,nout,n1,n2,thickness,NModes);
    neffTE(indPeriod,indDC)=neff(1);
    neffTM(indPeriod,indDC)=neff(2);
    nbwTE(indPeriod,indDC)=nBW(1);
    nbwTM(indPeriod,indDC)=nBW(2);    
    flag(indPeriod,indDC)=flagt;
 end
end

vPlot=vDC;
ilo=0

if ilo==1
h=figure, 
set(h,'pos',[ 235         419        1373         504])
subplot(131)
semilogx(vPlot,neffTE,vPlot,nbwTE,'g.','linewidth',2), 
hold on
semilogx(vPlot,neffTM,'--',vPlot,nbwTM,'r.','linewidth',2)
title('TE'), 
%hold on
subplot(132)
semilogx(vPlot,neffTM,'--',vPlot,nbwTM,'r.','linewidth',2)
title('TM'), 
subplot(133)
semilogx(vPlot,neffTE.^2-neffTM.^2,'linewidth',2)
else
h=figure, 
set(h,'pos',[ 235         419        1373         504])
subplot(131)
plot(vPlot,neffTE,vPlot,nbwTE,'g.','linewidth',2), 
hold on
plot(vPlot,neffTM,'--',vPlot,nbwTM,'r.','linewidth',2)
title('TE TM'), 
plot(vPlot, ones(size(vPlot))*sqrt(4.7),'w--')
%hold on
subplot(132)
plot(vPlot,neffTM,'--',vPlot,nbwTM,'r.','linewidth',2)
title('TM'), 
subplot(133)
plot(vPlot,neffTE.^2-neffTM.^2,'linewidth',2)

end
