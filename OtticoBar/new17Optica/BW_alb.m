clear
close all
clc

addpath E:\Dati\mvcsel\new10Optica
addpath E:\Dati\mvcsel\new10Optica\AlbertoBW
lambda=3300;
DC=0.5;
nin=1.41;
nout=1;
n1=3.43;
n2=1;
thickness=250;
NModes=11;
DC=.4;

vPeriod=logspace(2,3,21);
vTh=[20 40 80 160 320];
vTh=[10 20];
%vTh=[250];
%vPeriod=250;
for indDC=1:length(vTh)
   thickness=vTh(indDC);
 for indPeriod=1:length(vPeriod)
    Period=vPeriod(indPeriod);
    [neff,nBW,flagt] = f_EffectiveIndex(Period,lambda,DC,nin,nout,n1,n2,thickness,NModes);
    neffTE(indPeriod,indDC)=neff(1);
    neffTM(indPeriod,indDC)=neff(2);
    nbwTE(indPeriod,indDC)=nBW(1);
    nbwTM(indPeriod,indDC)=nBW(2);    
    flag(indPeriod,indDC)=flagt;
 end
end

h=figure, 
set(h,'pos',[ 235         419        1373         504])
subplot(131)
semilogx(vPeriod,neffTE,vPeriod,nbwTE,'g.','linewidth',2), 
title('TE'), 
%hold on
subplot(132)
semilogx(vPeriod,neffTM,'--',vPeriod,nbwTM,'r.','linewidth',2)
title('TM'), 
subplot(133)
semilogx(vPeriod,neffTE.^2-neffTM.^2,'linewidth',2)
