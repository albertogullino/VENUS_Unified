
clear all
close all

Lv=.168;
nv=2.3;

Lb=[.298 Lv];
nb=[1.3 nv];
kt=0;

la=1.6;
N=5;
ring=1;
rr=1;
rout=1;

Nv=3:15;
for k=1:length(Nv)
Nb=Nv(k);
[Ged,Gmd]=gaperdm(kt,0,la,Lv,nv,Lb,nb,Nb,rout,rr,0,[],[],ring);
Gev(k)=Ged;
Gmv(k)=Gmd;
end
figure, plot(Nv,-log10(1-abs(Gev).^2),Nv,-log10(1-abs(Gmv).^2),'o'), pausak

kv=linspace(0,1.5,101);

clear Gev Gmv
for k=1:length(kv)
kt=kv(k);
[Ged,Gmd]=gaperdm(kt,0,la,Lv,nv,Lb,nb,N,rout,rr,0,[],[],ring);
Gev(k)=Ged;
Gmv(k)=Gmd;
end

figure, plot(kv,abs(Gev),kv,abs(Gmv))