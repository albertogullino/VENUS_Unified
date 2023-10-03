%' entro tscatuu ', keyboard
%' entro tscatuu ', keyboard
%save sca1
' entro tscatuu ', 

  [V,autov]=eig(P);
  daut=diag(autov);
diaut=daut;



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
diaut
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

if iautov==1
 dim=diaut(fis);
 vim=ve(:,fis);
 [du,fis]=sort(-abs(dim));
 disort=[disort; dim(fis)];
 ves=[ves vim(:,fis)];
 ves1=ves;
 fis1=fis;
 if length(fis)<length(V)
 [fis]=find(imag(diaut)>=0 & real(diaut)>=0);
  dim=diaut(fis);
  vim=ve(:,fis);
  [du,fis]=sort(-abs(dim));
  disort=[disort; dim(fis)];
  ves=[ves vim(:,fis)];
 end 
else
 disort=dim;
 ves=vim;
end
%'dopo sort', keyboard 
% [du,fis]=sort(imag(diaut));
% disort=diaut(fis);
% ves=ve(:,fis);

lve=length(V)/2;
l1=1:lve;
l2=l1+lve;
 
ve=ves;
areg =disort(l1);

'areg', keyboard
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

is=inv(s12u);
Tu(l1,l1)=s21u-s22u*is*s11u;
Tu(l1,l2)=s22u*is;
Tu(l2,l1)=-is*s11u;
Tu(l2,l2)=is;
xTras=Tu;
x=Tu;
