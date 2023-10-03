function prif(Tt,iauto,dv,nv,fmlstot,rr,lambda,freq,KK,Gas,iLP,istop,iem,fape)
%close all
%ikba=0;
%if ikba==1
% load tr
%else
% load tr1
%end
iLPi=0;
Tc=Tt;
iproga=1
ii=find(iauto(:,1)==2)-1;
ispe=2:ii;
ispe1=ispe-1;
ds=dv(ispe)/1000;
ns=nv(ispe,1);
ifm1=fmlstot(ispe1,1);
ifm2=fmlstot(ispe1,2);
npmir=max(ifm2);
fid=find(ifm2==npmir)';
fiu=fid(1)-1;
puo=fliplr(1:fiu);
pum=fliplr(fid);
pui=fliplr(fid(end)+1:length(ispe));
ring=rr; %Refractive index of the infinite medium before the Bragg-Mirror
Lv=ds(pui); %Thickness of layer before the periodic section
rv=ns(pui); %Refractive index of the layer before the periodic section
Luv=ds(puo); %Thickness of layer after the periodic section
ruv=ns(puo); %Refractive index of the layer after the periodic section
rf=nv(1); %Refractive index of the infinite Medium at the end of the Bragg-Mirror
Lb=ds(pum);
rb=ns(pum);
ring=rr;
n=npmir;
icont=0;
if icont==1
' '
' UP '
[Luv ruv]
pausak
' '
' Mirr '
n
[Lb rb]
pausak
' '
' LOW '
[Lv rv]
keyboard
end

[Gae,Gam]=gaperd(KK,freq,lambda,Lv,rv,Lb,rb,n,rf,rr,iLPi,Luv,ruv,ring,fape);

[Gae0,Gam0]=gaemms(KK,freq,lambda,Lv,rv,Lb,rb,n,rf,rr,iLPi,Luv,ruv,ring);

isa=0;

if isa==1
KK1    =KK  ;
freq1  =freq ;
lambda1=lambda;
Lv1    =Lv;
rv1    =rv ;
Lb1    =Lb  ;
rb1    =rb;
n1     =n ;
rf1    =rf;
rr1    =rr;
iLP1   =iLPi;
Luv1   =Luv;
ruv1   =ruv;
ring1  =ring;

save pgamma1 KK1 freq1 lambda1 Lv1 rv1 Lb1 rb1 n1 rf1 rr1 iLP1 Luv1 ruv1 ring1

end


   if iproga==1
%     ' qui per gamma',  keyboard

     Ga1=diag(Gas);
     s=size(Tc);
     l1=s(1)/2;
     l2=s(1)/2+1;
     l3=s(1);
     T11=Tc(1:l1,1:l1);
     T12=Tc(1:l1,l2:l3);
     T21=Tc(l2:l3,1:l1);
     T22=Tc(l2:l3,l2:l3);
     It=inv(T22);
     S11=-It*T21;
     S12=It;
     S21=T11-T12*It*T21;
     S22=T12*It;

     GAU=(T12+T11*Ga1)*inv(T22+T21*Ga1);
%     GAU1=S11+S12*inv(eye(length(Ga1))-Ga1*S22)*Ga1*S21;
     GAU1=S22+S21*inv(eye(length(Ga1))-Ga1*S11)*Ga1*S12;
%     GAU=(T21+T22*Ga1)/(T11+T12*Ga1);
%     GAU1=(T21+T22*Ga1)*inv(T11+T12*Ga1);
%     mapab(GAU), pausak
     if iLP==1
      puD=1:length(KK);
     else
      puD=1:length(KK);
      puD=[puD puD+length(GAU)/2];
     end

     D=(diag(GAU));
     D=D(puD);
     D1=(diag(GAU1));
     D1=D1(puD);
    if iLP==0
     Dvte=D(1:end/2);
     Dvtm=D(1+end/2:end);
     Dte1=D1(1:end/2);
     Dtm1=D1(1+end/2:end);
     Dte=reshape(Dvte,length(KK),1);
     Dtm=reshape(Dvtm,length(KK),1);
     figure, plot(KK,abs(Dte),KK,abs(Dtm)),
     hold on, plot(KK,abs(Dte1),'o',KK,abs(Dtm1),'o'),
     hold on, plot(KK,abs(Gae0),'s',KK,abs(Gam0),'s'),
     hold on, plot(KK,abs(Gae),'.-',KK,abs(Gam),'.-'),
    else
     figure, plot(KK,abs(D)),
     hold on, plot(KK,abs(D1),'o'),
     if iem==1
      plot(KK,abs(Gae0),'s',KK,abs(Gae),'.-'),
     else
%      plot(KK,abs(Gam),'.'),
      plot(KK,abs(Gam0),'s',KK,abs(Gam),'.-'),
     end
    end
    title(' punto-linea: TMM con perdite, quad. senza, linea e tondi: CMT ')



   end
   if istop==1
    ' in prif '
    pausak
    close
   end
%' controllo'
%keyboard
%keyboard
