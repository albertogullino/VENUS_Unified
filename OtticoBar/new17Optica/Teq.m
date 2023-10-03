function [Tg1,Tg2,G1,G2]=Teq(KK,lambda,thicki,period,DC,r1,r2,r_in, r_out,rr,mbv,ifp)

lami=lambda;
thick=thicki;
dv0=period;
d1=dv0*DC;
d2=dv0*(1-DC);
%r_in=rr;
%r_out=rr;

NModi=11;

fim=find(mbv>0);
if length(mbv)==1
 mv=[0 2*mbv];
else
 mv=0:2*max(mbv);
end

k0=2*pi/lambda;
kr=k0*rr;

%kt=KK*kr;

tev=asin(KK*rr);


tetav=tev/pi*180;
phiv=linspace(0,180,31);

ifit=1;
if ifit==0
phiv=phiv(2:end);
end

%' qui Teq ing' , keyboard
ical=1;

if ical==1

for ite=1:length(tetav)
tetai=tetav(ite);
for ife=1:length(phiv)
phii=phiv(ife);


  [T11,T22,T12,T21,s11]=orta_skewTOTr(phii,tetai,r_in,r_out,r1,r2,d1,d2,thick,lambda,NModi,0);

  Gte0(ife,ite,1)=T11(1,1);
  Gtm0(ife,ite,1)=T11(2,2);
  Gtem0(ife,ite,1)=T11(1,2);
  Gtme0(ife,ite,1)=T11(2,1);
  Gte0(ife,ite,2)=T22(1,1);
  Gtm0(ife,ite,2)=T22(2,2);
  Gtem0(ife,ite,2)=T22(1,2); 
  Gtme0(ife,ite,2)=T22(2,1);
  Gte0(ife,ite,3)=T12(1,1);
  Gtm0(ife,ite,3)=T12(2,2);
  Gtem0(ife,ite,3)=T12(1,2);
  Gtme0(ife,ite,3)=T12(2,1);
  Gte0(ife,ite,4)=T21(1,1);
  Gtm0(ife,ite,4)=T21(2,2);
  Gtem0(ife,ite,4)=T21(1,2); 
  Gtme0(ife,ite,4)=T21(2,1);
  Gte0(ife,ite,5)=s11(1,1);
  Gtm0(ife,ite,5)=s11(2,2);
  Gtem0(ife,ite,5)=s11(1,2); 
  Gtme0(ife,ite,5)=s11(2,1);

end 
%'fine phi in Teq', keyboard
end 
save cal
else
load cal
end


%'vedo Teq', keyboard


if ifit==1

phiv0=phiv;
phiv=linspace(0,180,201)';

 for k=1:5
  for kt=1:length(tetav)
  Gte1(:,kt,k)=spline(phiv0',Gte0(:,kt,k),phiv);
  Gtm1(:,kt,k)=spline(phiv0',Gtm0(:,kt,k),phiv);
  Gtem1(:,kt,k)=spline(phiv0',Gtem0(:,kt,k),phiv);
  Gtme1(:,kt,k)=spline(phiv0',Gtme0(:,kt,k),phiv);
  end
 end
 phiv=phiv';
else
  Gte1=Gte0;
  Gtm1=Gtm0;
  Gtem1=Gtem0;
  Gtme1=Gtme0;
end

phir=phiv*pi/180;
dfi=diff(phir(1:2))/(pi);

for k=1:length(mv)
 m=mv(k);
 esc=cos(phir*m);
 ess=sin(phir*m);
  if ifit==1
   esc([1 end])=esc([1 end])/2;
   ess([1 end])=ess([1 end])/2;
  end 
 for kt=1:5      % ij
  for ki=1:4      %TE-TM  
   if ki==1
    G=Gte1(:,:,kt);
   elseif ki==2
    G=Gtm1(:,:,kt);
   elseif ki==3
    G=Gtem1(:,:,kt);
   elseif ki==4
    G=Gtme1(:,:,kt);    
   end
   RIc=esc*G*dfi;
   RIs=ess*G*dfi;
   Tc(ki,:)=RIc;
   Ts(ki,:)=RIs;
  end 
  
  Im1c(k,:,kt)=Tc(1,:);
  Im2c(k,:,kt)=Tc(2,:);
  Im3c(k,:,kt)=Tc(3,:);
  Im4c(k,:,kt)=Tc(4,:);
  Im1s(k,:,kt)=Ts(1,:);
  Im2s(k,:,kt)=Ts(2,:);
  Im3s(k,:,kt)=Ts(3,:);
  Im4s(k,:,kt)=Ts(4,:);  
 end %kt 
end  %mv 
mv0=mv; 

isk=1
if isk==0
Im1=Im1s;
Im2=Im2s;
Im3=Im3s;
Im4=Im4s;

mv=KK;
figure, semilogy(mv,abs(Im1(:,:,1))), hold on, 
semilogy(mv,abs(Im1(:,:,2)),'.'), semilogy(mv,abs(Im1(:,:,3)),'o'), 
semilogy(mv,abs(Im1(:,:,4)),'s'), 
title(' EE ')
pausak

figure, semilogy(mv,abs(Im2(:,:,1))), hold on, 
semilogy(mv,abs(Im2(:,:,2)),'.'), semilogy(mv,abs(Im2(:,:,3)),'o'), 
semilogy(mv,abs(Im2(:,:,4)),'s'), 
title(' MM ')
pausak

figure, semilogy(mv,abs(Im3(:,:,1))), hold on, 
semilogy(mv,abs(Im3(:,:,2)),'.'), semilogy(mv,abs(Im3(:,:,3)),'o'), 
semilogy(mv,abs(Im3(:,:,4)),'s'), 
title(' EM ')
pausak

figure, semilogy(mv,abs(Im1(:,:,1))), hold on, 
semilogy(mv,abs(Im2(:,:,1)),'.'), semilogy(mv,abs(Im3(:,:,1)),'o'), 
title(' misti ')
pausak

Im1=Im1c;
Im2=Im2c;
Im3=Im3c;
Im4=Im4c;

mv=KK;
figure, semilogy(mv,abs(Im1(:,:,1))), hold on, 
semilogy(mv,abs(Im1(:,:,2)),'.'), semilogy(mv,abs(Im1(:,:,3)),'o'), 
semilogy(mv,abs(Im1(:,:,4)),'s'), 
title(' EE ')
pausak

figure, semilogy(mv,abs(Im2(:,:,1))), hold on, 
semilogy(mv,abs(Im2(:,:,2)),'.'), semilogy(mv,abs(Im2(:,:,3)),'o'), 
semilogy(mv,abs(Im2(:,:,4)),'s'), 
title(' MM ')
pausak

figure, semilogy(mv,abs(Im3(:,:,1))), hold on, 
semilogy(mv,abs(Im3(:,:,2)),'.'), semilogy(mv,abs(Im3(:,:,3)),'o'), 
semilogy(mv,abs(Im3(:,:,4)),'s'), 
title(' EM ')
pausak

figure, semilogy(mv,abs(Im1(:,:,1))), hold on, 
semilogy(mv,abs(Im2(:,:,1)),'.'), semilogy(mv,abs(Im3(:,:,1)),'o'), 
title(' misti ')
pausak


end  %isk

clear Tc Ts
for kti=1:5;    %ij
%for kti=1;    %f-b
 for kem=1:4   %TE-TM
  if kem==1
   Imc=Im1c(:,:,kti);
  elseif kem==2
   Imc=Im2c(:,:,kti);
  elseif kem==3
   Imc=Im3s(:,:,kti);
%   Imc=Im3c(:,:,kti);
  elseif kem==4
   Imc=Im4s(:,:,kti);
%   Imc=Im4c(:,:,kti);
  end 
  %'kem', pausak
   Im=Imc;
     if kem<=2
      eec=Im(1,:)+Im(2,:);
      ees=Im(1,:)-Im(2,:);
     else 
      eec=-Im(2,:);
      ees=Im(2,:);  
     end  
     if kem==1
      Tc11=ees;
      Ts11=eec;
     elseif kem==2 
      Tc22=eec;
      Ts22=ees;
     elseif kem==3 
      Tc12=ees;
      Ts12=eec;
     elseif kem==4 
      Tc21=ees;
      Ts21=eec;
     end
 end %kem
  Tc=[diag(Tc11) diag(Tc12); diag(Tc21) diag(Tc22)];
  Ts=[diag(Ts11) diag(Ts12); diag(Ts21) diag(Ts22)];
  
  if kti==1
   T11c=Tc;
   T11s=Ts;
  elseif kti==2
   T22c=Tc;
   T22s=Ts;  
  elseif kti==3
   T12c=Tc;
   T12s=Ts;    
  elseif kti==4 
   T21c=Tc;
   T21s=Ts;  
  elseif kti==5
   G11c=Tc;
   G11s=Ts;     
  end   
end 

Ttotc=[T11c T12c; T21c T22c];
Ttots=[T11s T12s; T21s T22s];

[Gei,Gmi]=gaperdm(KK,0,lambda,[],[],[],[],0,r_in,rr,0,[],[],rr);
Gti=[Gei; Gmi];

[Geu,Gmu]=gaperdm(KK,0,lambda,[],[],[],[],0,rr,rr,0,[],[],r_out);
Gtu=[Geu; Gmu];

o=ones(size(Gti));
oM=diag(o);
gM=diag(Gtu);
Vu=[oM gM; gM oM];
gM=-diag(Gti);
Vi=[oM gM; gM oM];

Tg1=Vu*Ttotc*Vi;
Tg2=Vu*Ttots*Vi;


G1=G11c;
G2=G11s;

% save saG G1 G2
%'fine TeqD', keyboard
if ifp>=0

figure, plot(abs(diag(Ttots))), hold on, plot(abs(diag(Ttotc)),'r')
hold on
plot(abs(diag(Oo1)),'.'), hold on, plot(abs(diag(Oo2)),'r.'), 
title(' Diag 0')
pausak

figure, plot(abs(diag(Ttots,20))), hold on, plot(abs(diag(Ttotc,20)),'r')
hold on
plot(abs(diag(Oo1,20)),'.'), hold on, plot(abs(diag(Oo2,20)),'r.'), 
title(' Diag 20')
pausak


figure, plot(abs(diag(Ttots,10))), hold on, plot(abs(diag(Ttotc,10)),'r')
hold on
plot(abs(diag(Oo1,10)),'.'), hold on, plot(abs(diag(Oo2,10)),'r.'), 
title(' Diag 10')

end
