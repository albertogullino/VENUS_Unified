function [lam,g_gamma]=Flam_gamPack(L_in,n_i,iat,fst,la0,rr,kt,ibast,par_grat,g0);

if isfield(par_grat,'itetm')==1
iem=par_grat.itetm;
else
iem=1;
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
 if kt>0
  ite=1;
 else 
  ite=0;
 end 
else
 ite=0;
 kk=0;
end

la=la0;
k0=2*pi/la;
be=j*k0*n;
na=n_i(iat);


La=L_i(iat);


n1=n_i(1);
n3=n_i(end);

 ns=n_i(2:iat-1);
 nd=n_i(iat+1:end-1);

 Ls=L_i(2:iat-1);
 Ld=L_i(iat+1:end-1);

 fss=fst(2:iat-1,:);
 fsd=fst(iat+1:end-1,:);
 
 while n1==ns(1)
  ns=ns(2:end);
  Ls=Ls(2:end);
  fss=fss(2:end,:);
 end
 
 while n3==nd(end)
  nd=nd(1:end-1);
  Ld=Ld(1:end-1);
  fsd=fsd(1:end-1,:);
end 

%'contr na prima', keyboard



g0=0;
g_gamma=g0;


%'contr na', keyboard
GA=g0*1e-4/k0;
na=n_i(iat)+j*GA/2;


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
fx=fss;

 iBa=0;
if ibast>0
% 'ibast', keyboard
 if ibast<iat
  iBa=ibast-1;
 end
end

%' sopra', keyboard
in=1;
while in<=length(nx)
 irep=fx(in,1);
 npair=abs(fx(in,2));
 if irep==0
  irep=1;
 end
 Tloc=1;
 for inr=1:irep
 Li=Lx(in);
  if in==iBa 
 %'iBa in gamLu', keyboard
  [Te,Tm]=Tgrating(la,Li,rr,kk,par_grat);  
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
 Tloc=Ti*Tloc;
 in=in+1;
 %Tloc
 %pausak
end  % irep
%Ts=T;
 T=Tloc^npair*T;
% T
 if in>=1300 
 in
 'Ti ', pausak
 end
end

%' fine sopra', keyboard

Ts=T;

T=eye(2);
nx=nd;
Lx=Ld;
fx=fsd;

 iBa=0;
if ibast>0
% 'ibast', keyboard
 if ibast>iat
  iBa=ibast-length(Ls)-2;
 end
end
%iBa
%pausak

%' sopra', keyboard
in=1;
while in<=length(nx)
 irep=fx(in,1);
 npair=abs(fx(in,2));
 if irep==0
  irep=1;
 end
 Tloc=1;
 for inr=1:irep
 Li=Lx(in);
  if in==iBa 
 %'iBa in gamLu', keyboard
  [Te,Tm]=Tgrating(la,Li,rr,kk,par_grat);  
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
 Tloc=Ti*Tloc;
 in=in+1;
 %Tloc
 %pausak
end  % irep
%Ts=T;
 T=Tloc^npair*T;
% T
 if in>=1300 
 in
 'Ti ', pausak
 end
end

%'fine sotto', keyboard
Td=T;

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
 Tan=expm(Mv); 


Mi=[-1 -1; 1 1]*ZEv/2;
  if ite==0  
   delz=(n/ni)^4*ZEv*(kt/bb)^2;
   Mz=[-delz delz; -delz delz]/2;
   Mi=Mi+Mz;
  end 

Tai=Tan*Mi;

%Tn=Van2i*Td*Tan*Ts*Van1;
%Ti=Van2i*Td*Tai*Ts*Van1;

Tn=Td*Tan*Ts;
Ti=Td*Tai*Ts;

% per questi gamma, sempre valutarli da interno a esterno
% per questo Gai ha segno n_i scambiati

%Gas=-(n1-n2)/(n2+n1);
%Gad=-(n3-nu)/(n3+nu);

if kk==0
Gas=-(n1-n)/(n1+n);
Gad=-(n3-n)/(n3+n);
else
[Gas,Gm]=gaperdm(kk,0,la,[],[],[],[],0,n1,n,0,[],[],n);
 if ite==0
  Gas=Gm;
 end
[Gad,Gm]=gaperdm(kk,0,la,[],[],[],[],0,n3,n,0,[],[],n);
 if ite==0
  Gad=Gm;
 end
end

Mn=Gad*Tn(1,1)*Gas+Gad*Tn(1,2)-Tn(2,1)*Gas-Tn(2,2);
Mi=Gad*Ti(1,1)*Gas+Gad*Ti(1,2)-Ti(2,1)*Gas-Ti(2,2);
P=inv(Mi)*Mn;
lamf=P*n/(2*na*La*1e-4);
lam=2*lamf+g0;

%if kk~=0
%'in Flam_gamPack ', keyboard
%end

