function [Ez,Hz,nz,zet,Gaqw,NQW_ef,az,uL]=Flam_field(L_in,n_i1,iat,fiQW,la0,rr,kt,ibast,par_grat,g0,Nx_in,a_i);



if isfield(par_grat,'itetm')==1
iem=par_grat.itetm;
else
iem=1;
end

L_i=L_in/1000;
n_i=n_i1(:,1);

if ~exist('a_i')
 a_i=zeros(size(L_i));
end 

if exist('rr')
 n=rr;
else
 n=pi*50;
 n=2.5;
end

if exist('kt')
 kk=abs(kt);
 if kt>0
  ite=1;
 else 
  ite=0;
 end 
else
 ite=0;
 kk=0;
end



La=L_i(iat);


n1=n_i(1);
n3=n_i(end);
%n=n1;
%rr=n;

[Gei,Gmi,Tei,Tmi]=ga_simp(n1,rr,kt,rr);
[Geu,Gmu,Teu,Tmu]=ga_simp(rr,n3,kt,rr);

if iem==1
 gMu=Geu;
 gMi=Gei;
 Tru=Teu;
 Tri=Tei;
else
 gMu=Gmu;
 gMi=Gmi;
 Tru=Tmu;
 Tri=Tmi;
end

Vu=[1 -gMu; -gMu 1]/Tru;
Vi=[1 -gMi; -gMi 1]/Tri;


la=la0;
k0=2*pi/la;
be=j*k0*n;
na=n_i(iat);
nato=n_i1(iat,:);
anato=a_i(iat,:);

 ns=n_i(2:iat-1);
 nsto=n_i1(2:iat-1,:);
 asto=a_i(2:iat-1,:);
 nd=n_i(iat+1:end-1);
 ndto=n_i1(iat+1:end-1,:);
 adto=a_i(iat+1:end-1,:);
 Ls=L_i(2:iat-1);
 Ld=L_i(iat+1:end-1);

fiqw=zeros(size(n_i));
fiqw(fiQW)=1;

 
 fiQWs=fiqw(2:iat-1);
 fiQWd=fiqw(iat+1:end-1);

Lsop=L_i(1);
Lsot=L_i(end);

iwi=0;


if iwi==1
while ns(1)==n1
 Lsop=Lsop+Ls(1);
 ns=ns(2:end);
 Ls=Ls(2:end);
 nsto=nsto(2:end,:);
 asto=asto(2:end,:); 
 fiQWs=fiQWs(2:end);
end



%Lsot=L_i(end);
while nd(end)==n3
 Lsot=Lsot+Ld(end);
 nd=nd(1:end-1);
 Ld=Ld(1:end-1);
 ndto=ndto(1:end-1,:);
 adto=adto(1:end-1,:);  
 fiQWd=fiQWd(1:end-1);
end
end


%'contr na', keyboard

GA=0.5*g0*1e-4/k0;
na=real(n_i(iat))+j*GA;

    bb=conj(sqrt(1-kk.^2));
    ZEv=(1./bb);
     if ite==0
      ZEv=(bb);
     end
%    Ideltad=([ZEv; ZMv])/2;    
%    ZMv=real(mr.*bb);

T=eye(2);
nx=ns;
nxto=nsto;
axto=asto;
Lx=Ls;
fiqwi=fiQWs;

 iBa=0;
if ibast>0
% 'ibast', keyboard
  iBa=ibast;
 if ibast<iat
  iBa=ibast-1;
 end
end

Nx=Nx_in;
f_in=[0; 1];
uLong=[];
zed=linspace(0,Lsop,Nx+1);
zed=zed(2:end);
k01=sqrt((n1*k0)^2-(rr*kt)^2);

fie_in=[exp(-j*k01*(zed-Lsop))*f_in(1); exp(j*k01*(zed-Lsop))*f_in(2)];

fie=Vi*f_in;

fie_0=fie;
%'fie-0', keyboard
%fie(:,end)=Vi*fie(:,end)/sqrt(rr/n1);

%fie(:,end)=Vi*fie(:,end)*sqrt(rr/n1);

%fie(:,end)=Vi*fie(:,end);

fieIn=fie(:,end);
ze=ones(1,Nx)*diff(zed(1:2));
nz=ones(Nx,1)*n_i1(1,:);
az=ones(Nx,1)*a_i(1,:);

ze_in=ze;
nz_in=nz;
az_in=az;

%'nz Flam_filed', keyboard

   uLong=[uLong zeros(1,Nx)];

   uLongT=[zeros(1,Nx)];

%'nz Flam_filed', keyboard

%' sopra', keyboard
for in=1:length(nx)
 Nx=Nx_in;

% ' Nx', keyboard
 Li=Lx(in);
 if Li/Nx/la>.02
  Nx=fix(Li/(la*.02));
 end
% in
% pausak
 if in==iBa 
 %'fel', keyboard 
  [Te,Tm,Nx,neq]=TgratingCar(la,Li,rr,kk,par_grat,Nx_in);  
%  [Te,Tm,Nx,neq]=Tgrating(la,Li,0,kk,par_grat,Nx_in);  
  if iem==1
   Ti=Te;
  else
   Ti=Tm;
  end

%  ze=[ze Li];
%  nz=[nz 0];
%  fie=[fie Ti*fie(:,end)];
  dx=Li/Nx;
  nin=neq*ones(size(nz(end,:)));
  if(size(nz(end,:),2)>1)
   ain=0*ones(1,size(nz(end,:),2)-1);
  else
   ain=0;
  end
  %'ain', keyboard
  global Geometry Parameters  Options
  iAlb=Options.iAlb;
  iAlb=0;
 if iAlb==0  
  [fie,ze,nz,az]=fieldz_siyi(dx,nin,Ti,fie,ze,nz,Nx,az,ain);  
  uLong=[uLong zeros(1,Nx)];
 else
  [fied,zedu,nzdu,azdu]=fieldz_siyi(dx,nin,Ti,fie,ze,nz,Nx,az,ain);   
  lambdanm=la0*1000;
     [S11temp,S21temp,S12temp,S22temp,JunctionInfo,LayerInfo,HalfSpaceInfo]=f_EvalSMatrixRCWA(lambdanm,Parameters,Geometry,Options,kt/1000+1e-14,0);
 
% 	Rx(1,1)=S11temp(1,1);
% 	Rx(2,2)=S11temp(Parameters.NModes_in+1,Parameters.NModes_in+1);
% 	Rx(1,2)=S11temp(1,Parameters.NModes_in+1);
% 	Rx(2,1)=S11temp(Parameters.NModes_in+1,1);
% 	
% 	Tx(1,1)=S21temp(1,1);
% 	Tx(2,2)=S21temp(Parameters.NModes_in+1,Parameters.NModes_in+1);
% 	Tx(1,2)=S21temp(1,Parameters.NModes_in+1);
% 	Tx(2,1)=S21temp(Parameters.NModes_in+1,1);    
% 	
% 	Ga=Rx;
% 	Tr=Tx;
%     Gte(ikx,iky)=Ga(1,1); 
%     Gtm(ikx,iky)=Ga(2,2); 	

         a1inc=zeros(2*Parameters.NModes_in,1);
      if par_grat.itetm==1
         indMode=[1];
       else
         indMode=[Parameters.NModes_in+1];      
       end
       a1inc(indMode)=1;
        z=Geometry.z;
     %    [VContent,IContent,indMode]=f_PlotHarmonicContent_Forced(z,a1inc,Parameters,Geometry,JunctionInfo,LayerInfo);
 	[VContent]=f_PlotHarmonicContent_Forced(z,a1inc,Parameters,Geometry,JunctionInfo,LayerInfo); 
         Vee=VContent(:,indMode);  
        zd=z'/1000;
%        zd=-zd+zd(end);
        zd=diff(zd(1:2))*ones(size(z(1:end-1)));
        zd=[zd 0];
        Ve=abs(flipud(Vee));
        Ve=Ve/Ve(1)*fie(2,end);
        Vep=Ve';
        VeM=[zeros(size(Vep)); Vep ];
        VeM(:,end)=fied(:,end);
        fie=[fie VeM];
        ze=[ze zd];
        uLong=[uLong zeros(1,length(zd))];        
        nz=[nz; neq*ones(size(zd'))];
%    'qui campi', keyboard
   end

 else

 ni=nx(in);
 del=(ni^2-n^2)/(2*n^2)*ZEv;
 M=[-(bb+del) -del; del (bb+del)];
  if ite==0  
   delz=(1-(n/ni)^2)/2*ZEv*(kt/bb)^2;
   Mz=[-delz delz; -delz delz];
   M=M+Mz;
  end
 dx=Li/Nx;
 Mv=dx*be*M;
 Ti=expm(Mv); 
 nin=nxto(in,:);
 ain=axto(in,:);

  if in==1
   Tirr=Ti;
%   Ti=Vi*Ti*sqrt(rr/n1);
  end

 [fie,ze,nz,az]=fieldz_siyi(dx,nin,Ti,fie,ze,nz,Nx,az,ain);
  if  fiqwi(in)==1
   uLong=[uLong ones(1,Nx)];
  else
   uLong=[uLong zeros(1,Nx)];
  end
% 'ze ', pausak
 end
% 'Ti ', pausak
end


 ni=na;
 Li=La;
 del=(ni^2-n^2)/(n^2)*ZEv;
 Ma=[-(del/2) -del/2; del/2 (del/2)]*be*Li;
 M=[-(bb+del/2) -del/2; del/2 (bb+del/2)];
  if ite==0  
   delz=(1-(n/ni)^2)/2*ZEv*(kt/bb)^2;
   Mz=[-delz delz; -delz delz];
   M=M+Mz;
  end 
 dx=Li/Nx;
 Mv=dx*be*M;
 Ti=expm(Mv); 
 
 nin=nato;
 ain=anato;
 [fie,ze,nz,az]=fieldz_siyi(dx,nin,Ti,fie,ze,nz,Nx,az,ain);  
   uLong0=uLong*0;
   uLong0=[uLong0 ones(1,Nx)];
   uLong=[uLong ones(1,Nx)];

nx=nd;
nxto=ndto;
axto=adto;
Lx=Ld;
fiqwi=fiQWd;

%' dopo qw', keyboard
 iBa=0;
if ibast>0
% 'ibast', keyboard
 if ibast>iat
  iBa=ibast;
 end
end

for in=1:length(nx)
 Li=Lx(in);

 if in==iBa 
 %'fel', keyboard
  [Te,Tm]=TgratingCar(la,Li,rr,kk,par_grat);  
  if iem==1
   Ti=Te;
  else
   Ti=Tm;
  end
  ze=[ze Li];
  nz=[nz; zeros(size(nz(end,:)))];
  az=[az; zeros(size(nz(end,:)))];


  
  fie=[fie Ti*fie(:,end)];
   uLong0=[uLong0 zeros(1,Nx)];
   uLong=[uLong zeros(1,Nx)];
 else

 ni=nx(in);
 del=(ni^2-n^2)/(2*n^2)*ZEv;
 M=[-(bb+del) -del; del (bb+del)];
  if ite==0  
   delz=(1-(n/ni)^2)/2*ZEv*(kt/bb)^2;
   Mz=[-delz delz; -delz delz];
   M=M+Mz;
  end
 dx=Li/Nx;
 Mv=dx*be*M;
 Ti=expm(Mv); 
 


 nin=nxto(in,:);
 ain=axto(in,:);
 [fie,ze,nz,az]=fieldz_siyi(dx,nin,Ti,fie,ze,nz,Nx,az,ain);
   uLong0=[uLong0 zeros(1,Nx)]; 
 if  fiqwi(in)==1
   uLong=[uLong ones(1,Nx)];
  else
   uLong=[uLong zeros(1,Nx)];
  end 
 end
% 'sotto in', in
% pausak
 %'Ti ', pausak
end

%'sotto', keyboard

zep=linspace(0,Lsot,Nx+1);
zeu=zep(2:end);
zea=ones(size(zeu))*diff(zeu(1:2));
ze=[ze zea];

nzu=ones(size(zea'))*n_i1(end,:);
azu=zeros(size(zea'))*a_i(end,:);
%nz=[nz ones(size(zeu))*n3];
nz=[nz; nzu];
az=[az; azu];

k03=sqrt((n3*k0)^2-(rr*kt)^2);
%'qui campo u', keyboard, keyboard
%fied=[exp(-j*k03*zeu)*sum(fie(:,end)); exp(j*k03*zeu)*0];
fuv=Vu*fie(:,end);   %passo fuori dal risonatore, da impedenza rr a n3
fu=fie(:,end);   %passo fuori dal risonatore, da impedenza rr a n3
%fied=[exp(-j*k03*zeu)*fu(1); exp(j*k03*zeu)*fu(2)]*sqrt(rr/n3);
fied=[exp(-j*k03*zeu)*fuv(1); exp(j*k03*zeu)*fuv(2)];

fies=fie;

fie_in=[exp(-j*k01*(zed-Lsop))*f_in(1); exp(j*k01*(zed-Lsop))*f_in(2)];

fiedV=[exp(-j*k03*zeu)*sum(fuv); exp(j*k03*zeu)*0];
fiedI=[exp(-j*k03*zeu)*diff(fuv); exp(j*k03*zeu)*0];
%fieV=[fie_in*sqrt(rr/n1) fie fiedV*sqrt(rr/n3)];
%fieI=[fie_in/sqrt(rr/n1)  fie -fiedI/sqrt(rr/n3);];
sqr=sqrt(rr/n1);
fieV=[fie_in fie/sqr fiedV*sqrt(rr/n3)/sqr];
fieI=[fie_in/sqr^2  fie/sqr -fiedI/sqrt(rr/n3)/sqr];
%fie=[fie fied];

   uLongT=[uLongT ones(1,length(uLong0)-length(uLongT))];

   uLong0=[uLong0 zeros(1,Nx)];
   uLong=[uLong zeros(1,Nx)];
   uLongT=[uLongT zeros(1,Nx) ];




ze=[0 ze]; 
uLongT=[0 uLongT]; 
uLong0=[0 uLong0]; 
uLong=[0 uLong]; 

zet=cumsum(ze);
%' ferma fiez', keyboard
%' ferma fiez', keyboard

Ez=sum(fieV);
%Ez=Ez/max(Ez);
Hz=-diff(fieI);
%Hz=Hz/max(Hz);
%' ferma fiez', keyboard

nz=[nz_in; nz];
az=[az_in; az];

Ezn=abs(Ez.*nz(:,1).'.*uLongT).^2.*ze;

Ezf=abs(Ez).^2;


   
   Izt=sum(Ezn);
   Izqtot=sum(Ezn.*uLong);
   Izq0=sum(Ezn.*uLong0);
   
   fLi=find(diff(uLong)>0);
   fLu=find(diff(uLong)<0)+1;
   
   uL0=zeros(size(uLong));
   for kz=1:length(fLi)
    fi1=fLi(kz)+1:fLu(kz)-1;
    uLi=uL0;
    uLi(fi1)=1;
    Izqi(kz)=sum(Ezn.*uLi)/Izt;
    Edu=Ezf.*uLi;
    IEmax(kz)=max(Edu);
   end
   
   NQW_ef=Izqtot/Izq0;
   Gaqw=Izq0/Izt;
uL(1)=Izq0/Izt;   
uL(2)=Izqtot/Izt;
uL(2+[1:length(fLi)])=Izqi;   


Z0=377;
Zref=Z0/rr;
Po=real(Ez.*conj(Hz))/Zref;
Pqw=max(Po)-min(Po);
Putile=abs(Po(1));
Pback=abs(Po(end));
etaP=Putile/Pqw;

% legame guadagno-potenza

% old
d_qw=La*1e-4;  %in cm
n_qw=real(na);
k_gain=d_qw*n_qw/Z0;

%new
Zref=Z0/rr;
Po=real(Ez.*conj(Hz))/Zref;
Pqw=max(Po)-min(Po);
Putile=abs(Po(1));
Pback=abs(Po(end));
etaP=Putile/Pqw;

Teff=rr/n_qw*Putile*Zref/IEmax(2);


Pgain=k_gain*g0*max(IEmax);

Pgain/Pqw;


%' ferma fiez NUOVO', keyboard
%' ferma fiez', keyboard

return

figure, 
plot(zet,Po,'.'),

figure, 
subplot(211)
plot(zet,Po), 
subplot(212)
semilogy(zet,abs(Po)), 
title('Pointing'), pausak


figure, plot(zet,Ez,zet,Hz), 
title('parti reali'), pausak
figure, plot(zet,imag(Ez),zet,imag(Hz)), 
title('parti immaginarie'), pausak

figure, plot(zet,abs(Ez),zet,abs(Hz)), 
title('Moduli'), pausak