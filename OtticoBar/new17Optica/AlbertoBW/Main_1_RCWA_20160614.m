
clear
close all
clc

lambda=3300;
DutyCycle=0.5;
nin=5.03;
nout=1;
n1=3.43;
n2=1;
thickness=50;
NModes=11;

vPeriod=logspace(0,3,101);
%vPeriod=250;

for indPeriod=1:length(vPeriod)
    Period=vPeriod(indPeriod);
    [neff,nBW,flagt] = f_EffectiveIndex(Period,lambda,DutyCycle,nin,nout,n1,n2,thickness,NModes);
    neffTE(indPeriod)=neff(1);
    neffTM(indPeriod)=neff(2);
    flag(indPeriod)=flagt;
end

figure,semilogx(vPeriod,neffTE,vPeriod,ones(size(vPeriod))*nBW(1))

figure,semilogx(vPeriod,neffTM,vPeriod,ones(size(vPeriod))*nBW(2))