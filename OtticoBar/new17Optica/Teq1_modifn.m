function [Tg1,Tg2,G1,G2]=Teq1_modif(KK,lambda,thicki,par_grat,rr,mbv,ifp,segem)

lambda=9.954296420207723e-001;
isca=0;
NDIS_fi=51;

lami=lambda;
thick=thicki;
if isfield(par_grat,'th')==1
 thick=par_grat.th/1000;
end
dv0=par_grat.per;
DC=par_grat.DC;
d1=dv0*DC;
d2=dv0*(1-DC);
%r_in=rr;
%r_out=rr;

%par_grat.per=period;
%par_grat.r_in=r_in;
%par_grat.r_out=r_out;
%par_grat.r1=r1;
%par_grat.r2=r2;
%par_grat.DC=DC;

iBW=0;
iplotCa=0;
NModi=par_grat.NModi;
r_in=par_grat.r_in;
r_out=par_grat.r_out;
r1=par_grat.r1;
r2=par_grat.r2;


tev=asin(KK*rr/r_in);




tetav=tev/pi*180;
phiv=linspace(0,180,NDIS_fi)';

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
Slm(:,:,ife)=S;

end 

if ite==3
 'controllo Orta',
 keyboard
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


end
