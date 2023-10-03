close all
clear all
dF=2*pi/3;
dfi=[0 1 2]*dF+.1;
dff=[0 1 2]*dF+.3;
mv=[3:3:30];
clear iS iC
ico=0;
for km=mv
ico=ico+1
iS(ico)=sum(sin(km*dff)-sin(km*dfi))/km;
iC(ico)=sum(cos(km*dff)-cos(km*dfi))/km;
end

figure, plot(mv,iS,mv,iC), pausak

muv=0:20;

 veNa=muv*NaN;
ico=0;
for km=muv
ico=ico+1
 du=muv'+km;
 duf=du/3-fix(du/3);
 fi=find(duf==0);
 veNal=veNa;
 veNal(fi)=du(fi);
 ks(:,ico)=veNal; 
  du=abs(muv'-km);
  duf=du/3-fix(du/3);
  fi=find(duf==0);
  veNal=veNa;
  veNal(fi)=du(fi);
 kd(:,ico)=veNal
end
figure, plot(muv,ks,'s'),
figure, plot(muv,kd,'o')