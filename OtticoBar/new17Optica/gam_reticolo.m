function [Ge,Gm,Tre,Trm]=gam_ret(L_in,n_i,la0,n_in,n_out,rr,kt,ibast,par_grat);

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

n3=n_out;
 nd=n_i;
 Ld=L_i;
 

    bb=conj(sqrt(1-kk.^2));
    ZEv=(1./bb);
     if ite==0
      ZEv=(bb);
     end
%    Ideltad=([ZEv; ZMv])/2;    
%    ZMv=real(mr.*bb);


Te=eye(2);
Tm=eye(2);
nx=nd;
Lx=Ld;

  iBa=ibast;

for in=1:length(nx)
 Li=Lx(in);
 if in==iBa 
  [Tie,Tim]=TgratingGam(la,Li,rr,kk,par_grat);  
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
 Tie=expm(Mv); 
 Tim=Tie; 
 end

 Te=Tie*Te;
 Tm=Tim*Tm;
end


if kk==0
Gad=-(n3-n)/(n3+n);
else
[Gad,Gm]=gaperdm(kk,0,la,[],[],[],[],0,n3,n,0,[],[],n);
 if ite==0
  Gad=Gm;
 end
end




r_in=n_in;
r_out=n_out;
KK=kt;
lambda=la0;

%%%%%%%
[Gei,Gmi,Tei,Tmi]=gaperdm(KK,0,lambda,[],[],[],[],0,r_in,rr,0,[],[],rr);

[Geu,Gmu,Teu,Tmu]=gaperdm(KK,0,lambda,[],[],[],[],0,rr,rr,0,[],[],r_out);

Gti=-[Gei];
o=ones(size(Gti));
oM=diag(o);
Gtu=-[Geu];
gM=-diag(Gtu);
Vu=[oM gM; gM oM]/Teu;

gM=-diag(Gti);
Vi=[oM gM; gM oM]/Tei;

Tret=Vu*Te*Vi;


Gti=-[Gmi];
o=ones(size(Gti));
oM=diag(o);
Gtu=-[Gmu];
gM=-diag(Gtu);
Vu=[oM gM; gM oM]/Tmu;

gM=-diag(Gti);
Vi=[oM gM; gM oM]/Tmi;

Trmt=Vu*Tm*Vi;
%%%%%%%

Ge=-Tret(2,1)/Tret(2,2);
Gm=-Trmt(2,1)/Trmt(2,2);


Ge1=(Gad*Te(1,1)-Te(2,1))/(-Gad*Te(1,2)+Te(2,2));
Gm1=(Gad*Tm(1,1)-Tm(2,1))/(-Gad*Tm(1,2)+Tm(2,2));

Tre=Tret(1,:)*[1; Ge];
Trm=Trmt(1,:)*[1; Gm];
%'conty', keyboard