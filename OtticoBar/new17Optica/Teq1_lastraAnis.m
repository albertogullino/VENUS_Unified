function [Tg1,Tg2,G1,G2]=Teq_lastraAnis(KK,lambda,thicki,period,DC,r1,r2,r_in, r_out,rr,mbv,ifp,iplotCa)

isca=0;

lami=lambda;
thick=thicki;
dv0=period;
d1=dv0*DC;
d2=dv0*(1-DC);
%r_in=rr;
%r_out=rr;

par_grat.per=period;
par_grat.r_in=r_in;
par_grat.r_out=r_out;
par_grat.r1=r1;
par_grat.r2=r2;
par_grat.DC=DC;

iBW=0;

if length(mbv)==1
 mv=[0 2*mbv];
else
 mv=0:2:2*max(mbv);
end

tev=asin(KK*rr/r_in);

%keyboard


tetav=tev/pi*180;
phiv=linspace(0,180,31);

ifit=1;
if ifit==0
phiv=phiv(2:end);
end


k0=2*pi/lambda;

if iBW==0



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

if iplotCa==2
 phii=55;
 tetai=40;
 iplotCa=1;
end 


  [T11,T12,T21,T22,s11,s12,s21,s22,Tin]=LastraAn_mia(tetai, phii,lambda,thick,par_grat,iplotCa);
  %[T11,T22,T12,T21,s11,du,s21]=orta_skewTOTr(phii,tetai,r_in,r_out,r1,r2,d1,d2,thick,lambda,NModi,0);
%' guardo lastra isca', keyboard

Ti(:,:,ife)=Tin;
Tl=[T11 T12; T21 T22];
Tlm(:,:,ife)=Tl;
Sl=[s11 s12; s21 s22];
Slm(:,:,ife)=Sl;

if isca==1
T11=s11;
T12=s12;
T21=s21;
T22=s22;

S=[s11 s12; s21 s22];
U=S*S';
verun=sum(sum(abs(U-eye(length(U)))));
vev(ife)=verun;
end
%' unit ', keyboard

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

  for kt=1:length(tetav)
   for ir=1:4
    for ic=1:4
     Slm1(ir,ic,:)=spline(phiv0',Slm(ir,ic,:),phiv);
     Tlm1(ir,ic,:)=spline(phiv0',Tlm(ir,ic,:),phiv);
     Ti1(ir,ic,:)=spline(phiv0',Ti(ir,ic,:),phiv);
    end
   end 
 end
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

if isca==1
figure, plot(phiv0,vev), title('unitarieta')
' qui plot var fphi', keyboard
end


 esc=cos(phir*mbv);
 ess=sin(phir*mbv);
  if ifit==1
   esc([1 end])=esc([1 end])/2;
   ess([1 end])=ess([1 end])/2;
  end 

dfi1=dfi*2;
 Ed1(1,1,:)=ess.^2.*dfi1;
 Ed1(1,2,:)=-ess.*esc.*dfi1;
 Ed1(2,1,:)=-ess.*esc.*dfi1;
 Ed1(2,2,:)=esc.^2.*dfi1;

 Ed2(1,1,:)=esc.^2.*dfi1;
 Ed2(1,2,:)=ess.*esc.*dfi1;
 Ed2(2,1,:)=ess.*esc.*dfi1;
 Ed2(2,2,:)=ess.^2.*dfi1;
 
for irB=[0 2] 
 for icB=[0 2] 
  for ir=1:2
    IR=ir+irB;
   for ic=1:2
    IC=ic+icB;
    Es1(IR,IC,:)=Ed1(ir,ic,:);
    Es2(IR,IC,:)=Ed2(ir,ic,:);
   end
  end
 end
end
 
 
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
   Tc11=[];
   Tc12=[];
   Tc21=[];
   Tc22=[];
   Ts11=[];
   Ts12=[];
   Ts21=[];
   Ts22=[];   

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

   for kmodir=mbv;
    Tcd=[];
    Tsd=[];
    Nr=1;
    if kmodir==0
     Nr=2;
    end 
     
    for kmodic=mbv;
    Nc=1;
    if kmodic==0
     Nc=2;
    end
     N=1/sqrt(Nr*Nc);
     idif=kmodir-kmodic;
     Jey=j^(idif);
%     pausak
     iseg=sign(idif);
     
     if idif<0
      idif=abs(idif);
     end
     idifv=find(mv==idif);
     isum=kmodir+kmodic;
     isumv=find(mv==isum);     
     if kem<=2
%   eec=Im(1,:)+Im(2,:);
%   ees=Im(1,:)-Im(2,:);
      eec=Im(idifv,:)+Im(isumv,:);
      ees=Im(idifv,:)-Im(isumv,:);   
     else 
%   eec=-Im(2,:);
%   ees=Im(2,:);  
      eec=(-iseg*Im(idifv,:)+Im(isumv,:));   
      ees=(iseg*Im(idifv,:)+Im(isumv,:));
     end  
     eec=diag(eec)*N*Jey;
     ees=diag(ees)*N*Jey;
%      Tcd=[Tcd ees];
%      Tsd=[Tsd eec];
     if kem==1
%      Tc11=ees;
%      Ts11=eec;
      Tcd=[Tcd ees];
      Tsd=[Tsd eec];
%      ' 11 ', keyboard
     elseif kem==2 
%      Tc22=eec;
%      Ts22=ees;
      Tcd=[Tcd eec];
      Tsd=[Tsd ees];
     elseif kem==3 
%      Tc12=ees;
%      Ts12=eec;
      Tcd=[Tcd -ees];
      Tsd=[Tsd eec];
%      'kem=3'
%      pausak
     elseif kem==4 
%      Tc21=ees;
%      Ts21=eec;
      Tcd=[Tcd -eec];
      Tsd=[Tsd ees];
%      'kem=4'
%      pausak
     end      
    end   % fine righe
     if kem==1
      Tc11=[Tc11; Tcd];
      Ts11=[Ts11; Tsd];
     elseif kem==2 
      Tc22=[Tc22; Tcd];
      Ts22=[Ts22; Tsd];
     elseif kem==3 
      Tc12=[Tc12; Tcd];
      Ts12=[Ts12; Tsd];
     elseif kem==4 
      Tc21=[Tc21; Tcd];
      Ts21=[Ts21; Tsd];
     end    
   end   %fine colonne
%   ' cont kem', keyboard
 end %kem
  Tc=[Tc11 Tc12; Tc21 Tc22];
  Ts=[Ts11 Ts12; Ts21 Ts22];
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

%'dopo', keyboard

Ttotc=[T11c T12c; T21c T22c];    
Ttots=[T11s T12s; T21s T22s];



 %' qui ', keyboard


[Gei,Gmi,Tei,Tmi]=gaperdm(KK,0,lambda,[],[],[],[],0,r_in,rr,0,[],[],rr);
Geia=repmat(Gei./Tei,length(mbv),1);
Gmia=repmat(Gmi./Tmi,length(mbv),1);
Gti=[Geia; Gmia];
Teia=repmat(1./Tei,length(mbv),1);
Tmia=repmat(1./Tmi,length(mbv),1);
Tti=[Teia; Tmia];
tnD=-diag(Gti);
tD=diag(Tti);
Vi=[tD tnD; tnD tD];


[Gei,Gmi,Tei,Tmi]=gaperdm(KK,0,lambda,[],[],[],[],0,rr,rr,0,[],[],r_out);
%Geua=repmat(Geu,length(mbv),1);
%Gmua=repmat(Gmu,length(mbv),1);
%Gtu=[Geua; Gmua];
Geia=repmat(Gei./Tei,length(mbv),1);
Gmia=repmat(Gmi./Tmi,length(mbv),1);
Gtu=[Geia; Gmia];
Teia=repmat(1./Tei,length(mbv),1);
Tmia=repmat(1./Tmi,length(mbv),1);
Ttu=[Teia; Tmia];

tnD=-diag(Gtu);
tD=diag(Ttu);
Vu=[tD tnD; tnD tD];

fi=find(abs(Ttotc)<1e-12);
Ttotc(fi)=0;

fi=find(abs(Ttots)<1e-12);
Ttots(fi)=0;


Tg1=Vu*Ttotc*Vi;
Tg2=Vu*Ttots*Vi;

%' Teq1_ver contr', keyboard

G1=G11c;
G2=G11s;

P=Es1.*Tlm1;
T1=sum(P,3);
P=Es1.*Ti1;
Tic=sum(P,3);
P=Es1.*Slm1;
Sic=sum(P,3);

pu1=1:2;
pu2=3:4;
T1O=T1;
            T11=T1O(pu1,pu1);
            T12=T1O(pu1,pu2);
            T21=T1O(pu2,pu1);
            T22=T1O(pu2,pu2);
            
	    s12=inv(T22);
	    s21=T11-s12*T12*T21;
            s22=s12*T12;
            s11=-s12*T21;   

SdaT=[s11 s12; s21 s22];            
vun=sum(sum(abs(SdaT'*SdaT)-eye(size(SdaT))))

vunSic=sum(sum(abs(Sic'*Sic)-eye(size(SdaT))))
% funzioni base
		ro=linspace(0,10,20)';
		phi=pi/4;
		
		mpi=mbv+1;
		mme=mbv-1;
		Kt=rr*k0*KK;
	        fnup=cos(mpi*phi);
	    	gnup=-sin(mpi*phi);
	        fnum=cos(mme*phi);
	     	gnum=-sin(mme*phi);
		Exe=besselj(mpi,ro*Kt)*fnup+besselj(mme,ro*Kt)*fnum;
		Eye=-besselj(mpi,ro*Kt)*gnup+besselj(mme,ro*Kt)*gnum;
     		Exm=besselj(mpi,ro*Kt)*fnup-besselj(mme,ro*Kt)*fnum;
     		Eym=-besselj(mpi,ro*Kt)*gnup-besselj(mme,ro*Kt)*gnum;

		cpiuver=zeros(2,1);
		cmenover=zeros(2,1);
		cpiuver(:,1)=[0;1];
		cmenover(:,1)=s11*cpiuver(:,1);
                %s11dir=Sic(pu1,pu1);
		%cmenover(:,1)=s11dir*cpiuver(:,1);
		ind=1;
	        statoin=[cpiuver(:,ind); cmenover(:,ind)];
	        
	    	statoout=Tic*statoin;            
                stv=[statoin statoout];

        for seg=[1 -1]
            for isez=[1 2]
                sta=stv(:,isez);
     		SiE=sta(1)+sta(3);
     		SiM=sta(2)+sta(4);
     		Ex(:,isez)=SiE*Exe+seg*SiM*Exm;
     		Ey(:,isez)=SiE*Eye+seg*SiM*Eym;
            end
            if seg==1
             Epx=Ex;
             Epy=Ey;
            else
             Emx=Ex;
             Emy=Ey;
            end
     		figure, plot(ro,abs(Ex)),
     		title(' Campi Ex sinistra- destra  '), 
     		'segno = ', seg
     		pausak
        end


' Teq1_Anis', keyboard

else

 er1=r1^2;
 er2=r2^2;
 tt=period;
 t1=period*DC;
 t2=period*(1-DC);
 ex=tt*er1*er2/(t2*er1+t1*er2);
 ey=(t1*er1+t2*er2)/tt;

 rint=sqrt((ex+ey)/2);
 rintz=rint;
 n_di=(ex-ey);
 Delta_gr=n_di/rr^2;
 dos=thick; 
 kcav=k0*rr;
 
 
 bb=sqrt(1-KK.^2);
 be=[bb; bb];
 mm=mbv;
 nube=1;
 pasnu=2;
 segem=-1;
 %segem=1;
 nubesi=mm;
 nubesu=mm;
 sm1=length(KK);
 

 
     ZEv=1./bb;
     ZMv=bb;
     Ideltad=([ZEv; ZMv])/2;
     Ideltazd=([ZEv*0; ZMv.*KK.^2./(1-KK.^2)])/2;
     Idelta=diag(Ideltad);
     Ideltaz=diag(Ideltazd);
     
% anisotropia
' vale solo per 1 modo: no accopp !!!!!!!! '
       Ipeze=ZEv/2;
       Ipezm=ZMv/2;
       
       Idelteac=[];
       Ideltmac=[];
       Ideltea=[ Ipeze; Idelteac];
       Ideltma=[ Ipezm; Ideltmac];
       ddu=1;
       Idelteah=[ Ipeze/ddu; Idelteac];
       Ideltmah=[ Ipezm/ddu; Ideltmac];
       Ideltadiag=zeros(size(Ideltea));

       nupd=nube/2-fix(nube/2);
       if mm<=1
        isub=0;
        inmod=1+(nubesu-nubesi)/pasnu;
        inmodi=1;
       else
        isub=1;
        inmod=1+(nubesu-nubesi)/pasnu;
%        inmodi=(fix((nube)/2)-numodiacc)*length(KK)+1;
        inmodi=(fix((nubesi)/2))*length(KK)+1;
       end
%       iipv=[(numodi-2*numodiacc-isub)*length(KK)+1:numodi*length(KK)];
%       keyboard

       iip=[inmodi:round((nubesu+1)/2)*length(KK)];

          Ideldiae=[ Ipeze/ddu; Ideltadiag];
          Ideldiam=[ Ipezm/ddu; Ideltadiag];
%          Ideldiae=[ Ipeze; Ideltadiag];
%          Ideldiam=[ Ipezm; Ideltadiag];
          Ideltea=[ Ipeze; Idelteac];
          Ideltma=[ Ipezm; Ideltmac];

   

        % even
          Ideltee=diag(Idelteah,sm1)+diag(Ideltea,-sm1)+diag(Ideldiae);
          Ideltmm=-(diag(Ideltmah,sm1)+diag(Ideltma,-sm1))+diag(Ideldiam);
          Ideltem=-diag(Idelteah,sm1)+diag(Ideltea,-sm1)-diag(Ideldiae);
          Ideltme=diag(Ideltmah,sm1)-diag(Ideltma,-sm1)-diag(Ideldiam);

          Ideltee=Ideltee(iip,iip);
          Ideltmm=Ideltmm(iip,iip);
          Ideltem=segem*Ideltem(iip,iip);
          Ideltme=segem*Ideltme(iip,iip);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Matrice tale che Kt=Kt+Delta(strato)*Ideltaap
%% Cosi' anche le altre per i casi odd e con m pari

%          Ideltaap=[Ideltee -Ideltem; -Ideltme -Ideltmm]/4;
          Ideltaap=[Ideltee -Ideltem; -Ideltme Ideltmm]/4;

          Ideltee=diag(Idelteah,sm1)+diag(Ideltea,-sm1)-diag(Ideldiae);
          Ideltmm=-(diag(Ideltmah,sm1)+diag(Ideltma,-sm1))-diag(Ideldiam);
          Ideltem=-diag(Idelteah,sm1)+diag(Ideltea,-sm1)+diag(Ideldiae);
          Ideltme=diag(Ideltmah,sm1)-diag(Ideltma,-sm1)+diag(Ideldiam);

          Ideltee=Ideltee(iip,iip);
          Ideltmm=Ideltmm(iip,iip);
          Ideltem=segem*Ideltem(iip,iip);
          Ideltme=segem*Ideltme(iip,iip);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         PIdeltaam=[Ideltee -Ideltem; -Ideltme Ideltmm]/4;
          Ideltaam=[Ideltee -Ideltem; -Ideltme Ideltmm]/4; 
 
%' fine Ideltamm', keyboard          
 
       depe=(rint^2-rr^2)/rr^2;
       depez=1-rr^2/rintz^2;
       ck=-j*kcav*dos;
       Beli=ck*be;
       Mod=([Beli; -Beli]);
       
           Kan_gr=Ideltaap;
           KOt=ck*(Idelta*depe+Delta_gr*Kan_gr);
           KOz=ck*(Ideltaz*depez);
           P=[KOt+KOz  KOt-KOz; -(KOt-KOz)  -(KOt+KOz)]; 
           Pe=diag(Mod)+P;
           TeBW=expm(Pe);
           Tg1=TeBW;
           
           G1=-Tg1(2,1)/Tg1(2,2);
           
%           ' pausa TE', keyboard
           
           Kan_gr=Ideltaam;
           KOt=ck*(Idelta*depe+Delta_gr*Kan_gr);
           KOz=ck*(Ideltaz*depez);
           P=[KOt+KOz  KOt-KOz; -(KOt-KOz)  -(KOt+KOz)]; 
           Pe=diag(Mod)+P;
           TmBW=expm(Pe);           
           Tg2=TmBW;

           G2=-Tg2(2,1)/Tg2(2,2);



end



