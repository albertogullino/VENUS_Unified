function [T22]=T_zero(var,L_in,n_i,iat,icav,rr,kt,ibast,par_grat,ifp);

la0=var(1);
g0=var(2);


if exist('g0')==1
  G0=g0;
else
 G0=0;
end

if ~exist('ifp')
  ifp=-4;
end

if isfield(par_grat,'itetm')==1
iem=par_grat.itetm;
else
iem=1;
end
if iem==1
 ite=1;
else
 ite=0;
end 

L_i=L_in/1000;

if exist('rr')
 n=rr;
else
 n=pi*50;
 n=2.5;
end

if exist('kt')
 kk=abs(kt);
else
 ite=0;
 kk=0;
end

la=la0;
k0=2*pi/la;
be=j*k0*n;
na=real(n_i(iat));


La=L_i(iat);


n1=n_i(1);
n3=n_i(end);

 ns=n_i(2:iat-1);
 nd=n_i(iat+1:end-1);
 
 Ls=L_i(2:iat-1);
 Ld=L_i(iat+1:end-1);
 
 
 %while n1==ns(1)
 % ns=ns(2:end);
 % Ls=Ls(2:end);
 %end
 
 %while n3==nd(end)
 % nd=nd(1:end-1);
 % Ld=Ld(1:end-1);
%end 

%'contr na prima', keyboard


if exist('balle')

%if ~exist('g0')

 Luv=flipud(Ls);
 nuv=flipud(ns);
%   [GGe2,GGm2,TTe2,TTm2]=gaperd(KKs,freq,lambda,Lvbr,nvbr,...
%                   Lbb,nbb,nstratid,rfd,rr,iLP,Luvb,nuvb,ring,fapes); 
 [Gs,Gm]=gaperdm(kk,0,la,[],[],[],[],0,n1,n,0,Luv,nuv,na);
 if ite==0
  Gs=Gm;
 end
 Luv=Ld;
 nuv=nd;
 [Gd,Gm]=gaperdm(kk,0,la,[],[],[],[],0,n3,n,0,Luv,nuv,na);
  if ite==0
   Gd=Gm;
 end
 kc=k0*ns(end);
 Lc=sum(L_i(icav(1):icav(end)));
 Iat=La/2+sin(kc*La)/(2*kc);
 Ito=Lc/2+sin(kc*Lc)/(2*kc);
 Gat=Iat/Ito;
 Gat=2*La/Lc;
 g0=-Gat*1e4/(2*La)*log(abs(Gs*Gd));
end

g0=0;
g_gamma=g0;


%'contr na', keyboard
GA=g0*1e-4/k0;
na=real(n_i(iat))+j*GA/2;


%na=n_i(iat);


%n2=ns(1);
%nu=nd(end);
%Ga1=-(n2-n)/(n2+n);
%t1=sqrt(1-Ga1^2);
%Van1=[1 Ga1 ; Ga1 1]/t1;
%Van1i=[1 -Ga1 ; -Ga1 1]/t1;
%
%Ga2=-(nu-n)/(nu+n);
%t2=sqrt(1-Ga2^2);
%Van2=[1 Ga2 ; Ga2 1]/t2;
%Van2i=[1 -Ga2 ; -Ga2 1]/t2;

    bb=conj(sqrt(1-kk.^2));
    ZEv=(1./bb);
     if ite==0
      ZEv=(bb);
     end
%    Ideltad=([ZEv; ZMv])/2;    
%    ZMv=real(mr.*bb);

T=eye(2);
nx=ns;
Lx=Ls;

 iBa=0;
if ibast>0
% 'ibast', keyboard
 if ibast<iat
  iBa=ibast-1;
 end
end

%' sopra', keyboard
for in=1:length(nx)
 Li=Lx(in);
 if in==iBa 
% 'iBa in gamLu Car', keyboard
  [Te,Tm]=TgratingCar(la,Li,rr,kk,par_grat);  
  if iem==1
   Ti=Te;
  else
   Ti=Tm;
  end
  %' Lu Grat', keyboard
 else
 ni=nx(in);
 del=(ni^2-n^2)/(2*n^2)*ZEv;
 M=[-(bb+del) -del; del (bb+del)];
  if ite==0  
   delz=(1-(n/ni)^2)/2*ZEv*(kt/bb)^2;
   Mz=[-delz delz; -delz delz];
   M=M+Mz;
  end
 Mv=Li*be*M;
 Ti=expm(Mv); 

 end

 Tiv(:,:,in)=Ti;
 T=Ti*T;
 if in>=1 & in<=3 & ibast==0
%  in
%  keyboard
 end
  if in==-12 | in>=3600
   in
   T
  'Ti ', pausak
  end
end

Tg=1;
if length(ibast)>0
if ibast==0
for kgg=1:3
 Tg=Tiv(:,:,kgg)*Tg;
end
else
 Tg=Tiv(:,:,iBa);
end
end

%' ferma', keyboard

Ts=T;

T=eye(2);
nx=nd;
Lx=Ld;

 iBa=0;
if ibast>0
% 'ibast', keyboard
 if ibast>iat
  iBa=ibast-length(Ls)-2;
 end
end
%iBa
%pausak
for in=1:length(nx)
 Li=Lx(in);
 if in==iBa 
  [Te,Tm]=TgratingCar(la,Li,rr,kk,par_grat);  
%  iem
%  pausak
  if iem==1
   Ti=Te;
  else
   Ti=Tm;
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
 Mv=Li*be*M;
 Ti=expm(Mv); 
 end


 T=Ti*T;
end

%' fine sotto vecchio', keyboard
Td=T;


 GA=0.5*G0*1e-4/k0;
 ni=na+j*GA;
 
 Li=La;
 del=(ni^2-n^2)/(n^2)*ZEv;
 M=[-(bb+del/2) -del/2; del/2 (bb+del/2)];
  if ite==0  
   delz=(1-(n/ni)^2)/2*ZEv*(kt/bb)^2;
   Mz=[-delz delz; -delz delz];
   M=M+Mz;
  end 
 Mv=Li*be*M;
 TanG=expm(Mv); 


Tn=Td*TanG*Ts;

if kk==0
Gas=-(n1-n)/(n1+n);
Tes=2*sqrt(n*n1)/(n1+n);
Gad=-(n3-n)/(n3+n);
Ted=2*sqrt(n*n3)/(n3+n);
else
[Gas,Gm,Tes,Tmi]=gaperdm(kk,0,la,[],[],[],[],0,n1,n,0,[],[],n);
 if ite==0
  Gas=Gm;
  Tes=Tmi;
 end
[Gad,Gm,Ted,Tmi]=gaperdm(kk,0,la,[],[],[],[],0,n3,n,0,[],[],n);
 if ite==0
  Gad=Gm;
  Ted=Tmi;
 end
end

if iem==2
Vi=1/Tes*[1 -Gas; -Gas 1];
else
Vi=1/Tes*[1 Gas; Gas 1];
end
Vu=1/Ted*[1 -Gad; -Gad 1];
Vi=1/Tes*[1 Gas; Gas 1];


Mt=Vu*Tn*Vi;

Mn=Gad*Tn(1,1)*Gas+Gad*Tn(1,2)-Tn(2,1)*Gas-Tn(2,2);
%Mi=Gad*Ti(1,1)*Gas+Gad*Ti(1,2)-Ti(2,1)*Gas-Ti(2,2);

T22=log10(abs(Mt(2,2)));
T22=log10(abs(Mn));  % OK per TE
%'T22'

%keyboard
if ifp==-100
figure(111),
hold on
plot3(la0,G0,T22,'o'), view(3)
end
%'in Flam_gamLu ', keyboard

return


[Gei,Gmi,Tei,Tmi]=gaperdm(KK,0,lambda,[],[],[],[],0,r_ini,rr,0,[],[],rr);
Geia=repmat(Gei./Tei,length(mbv),1);
Gmia=repmat(Gmi./Tmi,length(mbv),1);
Gti=[Geia; Gmia];
Teia=repmat(1./Tei,length(mbv),1);
Tmia=repmat(1./Tmi,length(mbv),1);
Tti=[Teia; Tmia];
tnD=-diag(Gti);
tD=diag(Tti);
Vi=[tD tnD; tnD tD];


[Gei,Gmi,Tei,Tmi]=gaperdm(KK,0,lambda,[],[],[],[],0,rr,rr,0,[],[],r_outi);
Geu=Gei;
Gmu=Gmi;
Teu=Tei;
Tmu=Tmi;
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


%Tg1=Vu*Ttotc*Vi;
%Tg2=Vu*Ttots*Vi;

% scambio per girare verso
Tg1=inv(Vi)*Ttotc*inv(Vu);





%%%%%%%%%%%%%%%%%%%

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
 Mv=Li*be*M;
 Tan0=expm(Mv); 


Mi=[-1 -1; 1 1]*ZEv/2;
  if ite==0  
   delz=(n/ni)^4*ZEv*(kt/bb)^2;
   Mz=[-delz delz; -delz delz]/2;
   Mi=Mi+Mz;
  end 

autov=-G0*1e-4*La*na/rr;
Tai=autov*Tan0*Mi;

Tatot=Tan0+Tai
TTot=Td*Tatot*Ts
Tver=Td*TanG*Ts;
MnV=Gad*TTot(1,1)*Gas+Gad*TTot(1,2)-TTot(2,1)*Gas-TTot(2,2)
