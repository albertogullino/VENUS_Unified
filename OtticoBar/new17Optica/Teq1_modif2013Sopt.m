function [Tg1,Tg2,Tsav]=Teq1_modif2013S(KK,lambda,thicki,par_grat,rr,mbv,ifp,segem,icrit,isca,Tlm_in)

if exist('Tlm_in')==1
 icalcOrt=0;
 Tlms=Tlm_in;
else 
 icalcOrt=1;
end

if exist('icrit')==0
 icrit=0;
end

if exist('isca')==0
 isca=1;
end

%'qui icrit', keyboard
if icrit==1
iScat=2;     %0; media T;  1; media S e trasforma in T;  2, lavora tutto con S, Tg1 e' S
else
iScat=0;
end
%lambda=9.954296420207723e-001;
NDIS_fi=31;

phif=linspace(0,180,501)';

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

%'tev', keyboard




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

 if icalcOrt==1
	for ife=1:length(phiv)
	phii=phiv(ife);

%if iplotCa==1
% phii=55;
% tetai=40;
%end 


%  [T11,T12,T21,T22,s11,s12,s21,s22,Tin]=LastraAn_mia(tetai, phii,lambda,thick,par_grat,iplotCa);
   [T,S]=orta_skewTOTu(phii,tetai,r_in,r_out,r1,r2,d1,d2,thick,lambda,NModi,0,iplotCa,segem);
%' guardo lastra isca Optica', keyboard
if iScat==0
 Tlm(:,:,ife)=T;
else

% per usare la matrice scattering, devo riportarla all'impedenza di riferimento invece che 
% alle impedenze alle porte; uso il formalismo delle matrici di impedenza, vedi Graglia verde
  if isca==0
  ' Esse', keyboard
  z=-(S+eye(4))*inv(S-eye(4));
  k=KK(ite);
  kti=k*rr/real(r_in);
  fri=sqrt(1-kti^2);
  ktu=k*rr/real(r_out);
  fru=sqrt(1-ktu^2);
  sZver=diag(sqrt([[1/fri fri]/real(r_in) [1/fru fru]/real(r_out) ]));  % radice Z_vero (mezzi ingresso e uscita)
  fref=sqrt(1-k^2);
  sZrifI=diag(sqrt([[fref 1/fref] [fref 1/fref] ]*rr));  % inverso di radice Z_rif
  zetaref=sZrifI*sZver*z*sZver*sZrifI;
  Sr=(zetaref-eye(4))*inv(zetaref+eye(4));

  Tlm(:,:,ife)=Sr;
  else
  Tlm(:,:,ife)=S;
  end
 
 end %iScat
 
	end  %ife



if ifit==1


   for ir=1:4
    for ic=1:4
%     Slm1(ir,ic,:)=spline(phiv0',Slm(ir,ic,:),phiv);
      Tlm1(ir,ic,:)=spline(phiv,reshape(Tlm(ir,ic,:),length(phiv),1),phif);
%      Slm1(ir,ic,:)=spline(phiv,reshape(Slm(ir,ic,:),length(phiv),1),phif);
%     Ti1(ir,ic,:)=spline(phiv0',Ti(ir,ic,:),phiv);
    end
   end 
else
 Slm1=Slm;
 Tlm1=Tlm;
% Ti1=Ti;
end
else
 Tlm1=Tlms(:,:,:,ite);
end % icalcOrt

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

	 fd1(1,1,:)=-sr.*cc.*dfi1;
	 fd1(1,2,:)=-sr.*sc.*dfi1;
	 fd1(2,1,:)=cr.*cc.*dfi1;
	 fd1(2,2,:)=cr.*sc.*dfi1;

	 fd2(1,1,:)=-cr.*sc.*dfi1;
	 fd2(1,2,:)=cr.*cc.*dfi1;
	 fd2(2,1,:)=-sr.*sc.*dfi1;
	 fd2(2,2,:)=sr.*cc.*dfi1;	

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
	        if ite==length(KK)
%	         ' colomma', keyboard
	        end 
     end     
        T1d=[T1d; T1r];
        T2d=[T2d; T2r];     
%        ' cont qui', keyboard
   end   % fine righe

 T1s{ite}=T1d(puord,puord);   
 T2s{ite}=T2d(puord,puord);   
  %'fine kk ite', keyboard
 if icalcOrt==1
  Tsav(:,:,:,ite)=Tlm1;
 else 
  %'fine kk vedo', keyboard  
 end
end   % fine KK
G1=0;
G2=0;

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


if iScat<2

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

else  %iScat

	Tg1=Ttotc;
	Tg2=Ttots;
	if isca==1

	 k=repmat(KK,length(mbv),1);
	 kti=k.*rr/real(r_in);
	 fri=msqrt(1-kti.^2);
	 ktu=k.*rr/real(r_out);
	 fru=msqrt(1-ktu.^2);
	 sZver=diag(sqrt([[1./fri; fri]/real(r_in); [1./fru; fru]/real(r_out) ]));  % radice Z_vero (mezzi ingresso e uscita)
	 fref=msqrt(1-k.^2);
	 sZrifI=diag(sqrt([[fref; 1./fref]; [fref; 1./fref] ]*rr));  % inverso di radice Z_rif
	 Oo1=Tg1;
  	 lo=length(Oo1);
	 S=Oo1;
	 Id=eye(lo);
	 z=-(S+Id)*inv(S-Id);
	 zetaref=sZrifI*sZver*z*sZver*sZrifI;
	 Tg1=(zetaref-Id)*inv(zetaref+Id);

	 Oo1=Tg2;
	 S=Oo1;
	 z=-(S+Id)*inv(S-Id);
	 zetaref=sZrifI*sZver*z*sZver*sZrifI;
	 Tg2=(zetaref-Id)*inv(zetaref+Id);	 
	
	end

end

'fine Vecchio Optica', keyboard
return


end
