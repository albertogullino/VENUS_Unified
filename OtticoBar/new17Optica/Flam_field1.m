function [fie,nz,zet,Gaqw,NQW_ef,az]=Flam_field1(L_in,n_i1,iat,fiQW,la0,rr,kt,ibast,par_grat,g0,Nx_in,a_i);



if isfield(par_grat,'itetm')==1
iem=par_grat.itetm;
else
iem=1;
end

L_i=L_in/1000;
n_i=n_i1(:,1);

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
Lsot=L_i(1);

iwi=1;

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
na=n_i(iat)+j*GA;

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
 if ibast<iat
  iBa=ibast-1;
 end
end

sNx=Nx_in/1000;
fie=[0; 1];
%zed=linspace(0,Lsop,Nx+1);
%' Lsop', keyboard
Nx=20;
zed=linspace(0,Lsop,Nx+1);
zed=zed(2:end);
fie=[exp(-j*n1*k0*zed)*fie(1); exp(j*n1*k0*zed)*fie(2)];
ze=ones(1,Nx)*diff(zed(1:2));
nz=ones(Nx,1)*n_i1(1,:);
az=ones(Nx,1)*a_i(1,:);

%'nz Flam_filed', keyboard

   uLong=[zeros(1,Nx)];

   uLongT=[zeros(1,Nx)];

%'nz Flam_filed', keyboard

%' sopra', keyboard
for in=1:length(nx)
 Li=Lx(in);
 if Li>0
 sNx=Nx_in/1000;
 zed=[0:sNx:Li];
 Nx=length(zed);
 %in
 %pausak

 if in==iBa 
  if length(par_grat.r1)>1 
   Li=sum(par_grat.th)/1000;
   zed=[0:sNx:Li];
   Nx=length(zed);
  end
  [Te,Tm,Nx,neq]=Tgrating(la,Li,rr,kk,par_grat,Nx_in);  
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
  ze=[ze 0];
  fie=[fie fie(:,end)];
  nz=[nz; nin];
  az=[az; ain];
%  'ain', keyboard
  [fie,ze,nz,az]=fieldz_siyi(dx,nin,Ti,fie,ze,nz,Nx,az,ain);  
  uLong=[uLong zeros(1,Nx+1)];
%   'qui campi', keyboard
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
  ze=[ze 0];
  fie=[fie fie(:,end)];
  nz=[nz; nin];
  az=[az; ain]; 
 [fie,ze,nz,az]=fieldz_siyi(dx,nin,Ti,fie,ze,nz,Nx,az,ain);
  if  fiqwi(in)==1
   uLong=[uLong ones(1,Nx+1)];
   %' uL =1 ', keyboard
  else
   uLong=[uLong zeros(1,Nx+1)];
  end
% 'ze ', pausak
 end
% 'Ti ', pausak
end %>0
end


%'att prima ', keyboard

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
 Nxq=10; 
 dx=Li/Nxq;
 Mv=dx*be*M;
 Ti=expm(Mv); 
 nin=nato;
 ain=anato;
  ze=[ze 0];
  fie=[fie fie(:,end)];
  nz=[nz; nin];
  az=[az; ain];  
 [fie,ze,nz,az]=fieldz_siyi(dx,nin,Ti,fie,ze,nz,Nxq,az,ain);  
%' att 1',keyboard
   uLong0=uLong*0;
   uLong0=[uLong0 ones(1,Nxq+1)];
   uLong=[uLong ones(1,Nxq+1)];
   
%' att ',keyboard
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
 
 if Li>0
 sNx=Nx_in/1000;
 zed=[0:sNx:Li];
 Nx=length(zed); 
% Li=Lx(in);
 if in==iBa 
  [Te,Tm]=Tgrating(la,Li,rr,kk,par_grat);  
  if iem==1
   Ti=Te;
  else
   Ti=Tm;
  end
  ze=[ze Li];
  nz=[nz; zeros(size(nz(end,:)))];
  az=[az; zeros(size(nz(end,:)))];
  
  fie=[fie Ti*fie(:,end)];
   uLong0=[uLong0 zeros(1,Nx+1)];
   uLong=[uLong zeros(1,Nx+1)];
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
%  Mt=Li*be*M;
%  Tt=expm(Mt);  
 nin=nxto(in,:);
 ain=axto(in,:);
  ze=[ze 0];
  fie=[fie fie(:,end)];
  nz=[nz; nin];
  az=[az; ain];
 [fie,ze,nz,az]=fieldz_siyi(dx,nin,Ti,fie,ze,nz,Nx,az,ain);
   uLong0=[uLong0 zeros(1,Nx+1)]; 
 if  fiqwi(in)==1
   uLong=[uLong ones(1,Nx+1)];
  else
   uLong=[uLong zeros(1,Nx+1)];
  end 
 end
 %'sotto'
 %in,
 %pausak
 %'Ti ', pausak
 end
end

Nx=20;
zep=linspace(0,Lsot,Nx+1);
zeu=[zep(1:end)];
zea=[0 ones(1,Nx)*diff(zeu(1:2))];
zev=ze;
ze=[ze zea];

nzu=ones(size(zea'))*n_i1(end,:);
azu=zeros(size(zea'))*a_i(end,:);
%nz=[nz ones(size(zeu))*n3];
nz=[nz; nzu];
az=[az; azu];


fied=[exp(-j*n3*k0*zeu)*sum(fie(:,end)); exp(j*n3*k0*zeu)*0];
fie=[fie fied];

   uLongT=[uLongT ones(1,length(uLong0)-length(uLongT))];

   uLong0=[uLong0 zeros(1,Nx+1)];
   uLong=[uLong zeros(1,Nx+1)];
   uLongT=[uLongT zeros(1,Nx+1) ];


%'qui campo u', keyboard

Ezn=abs(sum(fie).*nz(:,1).'.*uLongT).^2.*ze;
I=abs(sum(fie)).^2;


   
   Izt=sum(Ezn);
   Izqtot=sum(Ezn.*uLong);
   Izq0=sum(Ezn.*uLong0);
   
   NQW_ef=Izqtot/Izq0;
   Gaqw=Izq0/Izt;
   

zet=cumsum(ze);
%' ferma fiez', keyboard
%' ferma fiez', keyboard
%' ferma fiez', keyboard
