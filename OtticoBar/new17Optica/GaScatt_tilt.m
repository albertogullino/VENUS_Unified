function [E,S]=GaScat_Tilt(P);
%'entro gae', pausak


  [V,autov]=eig(P);
  daut=diag(autov);
%  ver=P*V-V*autov;

% Mdef = V * diag(exp(daut)) / V;
% Mpad=expm(P); 
% map(log10(abs(Mdef-Mpad))), pausak

  
diaut=daut;

%'qui', keyboard

mar=max(abs(real(diaut)));
sogliaLog=10;
% ripulisco autovalori
for ks=1:length(diaut)
 a=abs(real(diaut(ks)));
 b=abs(imag(diaut(ks)));
 if abs(log10(a/b))>sogliaLog
  if (log10(a/b))>sogliaLog
   diaut(ks)=real(diaut(ks));
  else 
   diaut(ks)=j*imag(diaut(ks));
  end
 end
end
%diaut
%'autovalori ripuliti', pausak

%'sopo pul', pausak
ve=V;
disort=[];
ves=[];
[fis]=find(imag(diaut)<=0 & real(diaut)<=0);
iautov=1;

if mar<.1
 [du,fis0]=sort(imag(diaut));
 lme=length(du)/2;
 pum=1:lme;
 fis=[fis0(pum); flipud(fis0(pum+lme))];
 iautov=0;
end


 dim=diaut(fis);
 vim=ve(:,fis);
%' qui fis', keyboard

if mean(abs(real(diaut)))>.1
 iautov=2;
end
iautov=1;
if iautov==2
 
 firm=find(real(diaut)>0);
 [du,firmu]=sort(abs(real(diaut(firm))));
 fit1=(firm(firmu));
 fitotr=[fit1];
 
  firm=find(real(diaut)<0);
  [du,firmu]=sort(abs(real(diaut(firm))));
  fit1=(firm(firmu));
 fitotp=[fit1];
 
 fis=[fitotr; fitotp];
 disort=diaut(fis);
 ves=V(:,fis);
elseif iautov==1
 fir=find(abs(real(diaut))>abs(imag(diaut)));
 fii=find(abs(real(diaut))<abs(imag(diaut)));
 diar=diaut(fir);
 firm=find(real(diar)<0);
 [du,firmu]=sort(abs(real(diar(firm))));
 fit1=fir(firm(firmu));
 diai=diaut(fii);
 fiim=find(imag(diai)<0);
 [du,fiimu]=sort((imag(diai(fiim))));
 fit2=fii(fiim(fiimu));
 fitotr=[fit2; fit1];
 
  firm=find(real(diar)>0);
  [du,firmu]=sort(abs(real(diar(firm))));
  fit1=fir(firm(firmu));
  fiim=find(imag(diai)>0);
  [du,fiimu]=sort(-(imag(diai(fiim))));
  fit2=fii(fiim(fiimu));
 fitotp=[fit2; fit1];
 
 fis=[fitotr; fitotp];
 disort=diaut(fis);
 ves=V(:,fis);
else
 disort=dim;
 ves=vim;
end


% [du,fis]=sort(imag(diaut));
% disort=diaut(fis);
% ves=ve(:,fis);

lve=length(V)/2;
l1=1:lve;
l2=l1+lve;
 
ve=ves;
areg0 =disort(l1);
%areg=real(areg0)-j*abs(imag(areg0));
areg=areg0;

apro=disort(l2);

%'dopo sort', keyboard 
%'areg', keyboard
vei=inv(ve);
lve=length(ve)/2;
l1=1:lve;
l2=l1+lve;

Ta=ve;
Ta11=Ta(l1,l1);
Ta12=Ta(l1,l2);
Ta21=Ta(l2,l1);
Ta22=Ta(l2,l2);
iT=inv(Ta11);
sa22=-iT*Ta12;
sa21=iT;
sa12=Ta22-Ta21*iT*Ta12;
sa11=Ta21*iT;

Ta=vei;
Ta11=Ta(l1,l1);
Ta12=Ta(l1,l2);
Ta21=Ta(l2,l1);
Ta22=Ta(l2,l2);
iT=inv(Ta11);
sb22=-iT*Ta12;
sb21=iT;
sb12=Ta22-Ta21*iT*Ta12;
sb11=Ta21*iT;


Ev=(exp(areg));
%Ev=(exp(apro));
fie=find(abs(Ev)<1e-8);
Ev(fie)=0;
E=diag(Ev);
uno=diag(ones(lve,1));

sab11=E*sb11*E;
sab12=E*sb12;
sab21=sb21*E;
sab22=sb22;

De2=inv(uno-sa22*sab11);
De1=inv(uno-sab11*sa22);
s11u=sa11+sa12*De1*sab11*sa21;
s12u=sa12*De1*sab12;
s21u=sab21*De2*sa21;
s22u=sab22+sab21*De2*sa22*sab12;


S=[s11u s12u; s21u s22u];
