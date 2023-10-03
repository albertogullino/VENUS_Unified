function [Tg1,Tg2,G1,G2]=Teq1_modif(KK,lambda,thicki,period,DC,r1,r2,r_in, r_out,rr,mbv,ifp,segem)

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
iplotCa=0;
NModi=11;


tev=asin(KK*rr/r_in);



tetav=tev/pi*180;
phiv=linspace(0,180,31)';

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


for ite=1:length(tetav)
tetai=tetav(ite);

for ife=1:length(phiv)
phii=phiv(ife);

%if iplotCa==1
% phii=55;
% tetai=40;
%end 


%  [T11,T12,T21,T22,s11,s12,s21,s22,Tin]=LastraAn_mia(tetai, phii,lambda,thick,par_grat,iplotCa);
   [T,S]=orta_skewTOTu(phii,tetai,r_in,r_out,r1,r2,d1,d2,thick,lambda,NModi,0,iplotCa,segem);
%' guardo lastra isca', keyboard
Tlm(:,:,ife)=T;
%Slm(:,:,ife)=S;

end 

if ifit==1

phif=linspace(0,180,201)';
   for ir=1:4
    for ic=1:4
%     Slm1(ir,ic,:)=spline(phiv0',Slm(ir,ic,:),phiv);
      Tlm1(ir,ic,:)=spline(phiv,reshape(Tlm(ir,ic,:),length(phiv),1),phif);
%     Ti1(ir,ic,:)=spline(phiv0',Ti(ir,ic,:),phiv);
    end
   end 
else
%Slm1=Slm;
 Tlm1=Tlm;
% Ti1=Ti;
end

phir=phif*pi/180;
dfi=diff(phir(1:2))/(pi);
dfi1=dfi*2*ones(size(phir));
if ifit==1
   dfi1([1 end])=dfi1([1 end])/2;
end 

 puop=[];
 puor=[];
 pd0=[1];
 for ik=1:length(mbv)
  puop=[puop pd0+(ik-1)*4];
  puor=[puor pd0+(2+(ik-1)*4)];
 end 
 puord=[puop puop+1 puor puor+1];
 %' puor', keyboard

    T1d=[];
    T2d=[]; 
    
   for kmodir=mbv;
    T1r=[];
    T2r=[];
    Nr=1;
    if kmodir==0
     Nr=2;
    end 
	 cr=cos(phir*kmodir);
	 sr=sin(phir*kmodir);
     
    for kmodic=mbv;
    Nc=1;
    if kmodic==0
     Nc=2;
    end
     N=1/sqrt(Nr*Nc);
     idif=kmodir-kmodic;     
     Jey=j^(idif);


	 cc=cos(phir*kmodic);
	 sc=sin(phir*kmodic);	 

	 fd1(1,1,:)=sr.*sc.*dfi1;
	 fd1(1,2,:)=-sr.*cc.*dfi1;
	 fd1(2,1,:)=-cr.*sc.*dfi1;
	 fd1(2,2,:)=cr.*cc.*dfi1;

	 fd2(1,1,:)=cr.*cc.*dfi1;
	 fd2(1,2,:)=cr.*sc.*dfi1;
	 fd2(2,1,:)=sr.*cc.*dfi1;
	 fd2(2,2,:)=sr.*sc.*dfi1;	


% fattori angolari traformata cilindrica

		for irB=[0 2] 
		 for icB=[0 2] 
		  for ir=1:2
		    IR=ir+irB;
		   for ic=1:2
		    IC=ic+icB;
		    fa1(IR,IC,:)=fd1(ir,ic,:);
		    fa2(IR,IC,:)=fd2(ir,ic,:);
		   end
		  end
		 end
		end
% integrali

		P=fa1.*Tlm1*N*Jey;
		T1i=sum(P,3);
		P=fa2.*Tlm1*N*Jey;
		T2i=sum(P,3);
		
	        T1r=[T1r T1i];
	        T2r=[T2r T2i];
	        %' colomma', keyboard
     end     
        T1d=[T1d; T1r];
        T2d=[T2d; T2r];     
   end   % fine righe

 T1s{ite}=T1d(puord,puord);   
 T2s{ite}=T2d(puord,puord);   
 %'fine Kk', keyboard

end   % fine KK
G1=0
G2=0

lK=length(KK);
siz=length(mbv)*4;
M1=zeros(siz,siz);
M2=M1;
pu=0:lK:lK*length(mbv)*4-1;



for l=1:lK
 pui=pu+l;
 M1(pui,pui)=T1s{l};
 M2(pui,pui)=T2s{l};
end
Ttotc=M1;
Ttots=M2;


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

'fine nuovo', keyboard
return

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
                s11dir=Sic(pu1,pu1);
		cmenover(:,1)=s11dir*cpiuver(:,1);
		ind=1;
	        statoin=[cpiuver(:,ind);cmenover(:,ind)];
	        
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



