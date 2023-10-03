%load qua
%load ott
%load pro
close all

KK    = Ppol.KK    ;
KKt   = Ppol.KKt   ;
KKv   = Ppol.KKv   ;
numodi= Ppol.numodi;
pasnu = Ppol.pasnu ;
lbv   = Ppol.lbv   ;
kcav0 = Ppol.kcav0 ;
x     = Ppol.x     ;
fian  = Ppol.fian  ;
Nx    = Ppol.Nx    ;
mbv   = Ppol.mbv   ;


nube=Ppol.nube;
a0ref=Ppol.a0ref;
xroJ=Ppol.xroJ;
axtot=Ppol.aitot;
aytot=Ppol.aitot;
icampi=Ppol.icampi;


kcav0=Ppol.kcav0;
pola0=Ppol.pola0;
iLP=Ppol.iLP;
fimaxi=Ppol.fimaxi;
lfi_inp=Ppol.lfi_inp;
ldap=Ppol.ldap;
Pus=Ppol.Pus;
isi=Ppol.isi;
rr=Ppol.rr;
rfu=Ppol.rfu;
mbv1=Ppol.mbv1;
pimu=Ppol.pimu;
Fint=Ppol.Fint;
iFFsect=Ppol.iFFsect;
iFF=Ppol.iFF;
iFFte=Ppol.iFFte;
iant=Ppol.iant;
pes=Ppol.pes;
z0c=Ppol.z0c;
lambda=Ppol.lambda;
nk1max=Ppol.nk1max;
avero=max(max(axtot));
bvero=max(max(aytot));
alim=max(KK);
nK_dis=nk1max;
numodiacc=numodi-1;
Cug=Ppol.Cug;
aax=Ppol.aax;

istrumix=1;
ifp=-10;
idyn=0;
iraff=0;
isalva=0;
i2D=3;
if iLP==0
 iztm=1;
else
 iztm=0;
end

if length(pola0)==2
 Azvet1= Ppol.Az1;
 Azvet2= Ppol.Az2;
 Azvetf1= Ppol.Azf1;
 Azvetf2= Ppol.Azf2;
 Gvet1 = Ppol.G1 ;
 Gvet2 = Ppol.G2 ;
 alvet1 = Ppol.A1 ;
 alvet2 = Ppol.A2 ;
 ipolar=2;
elseif length(pola0)==1
 if pola0==-1
  ipolar=-1;
  Azvet1= Ppol.Az1;
  Azvetf1= Ppol.Azf1;
  Gvet1 = Ppol.G1 ;
  alvet1 = Ppol.A1 ;
 else
  ipolar=1;
  Azvet2= Ppol.Az2;
  Azvetf2= Ppol.Azf2;
  Gvet2 = Ppol.G2 ;
  alvet2 = Ppol.A2 ;
 end
end

imod=0;
imm=1;
%for pola=[-1 1]
for pola=pola0

  if pola==-1
   Azvet=Azvet1;
   Azvetf=Azvetf1;
   Gvet=Gvet1;
   alvet=alvet1;
  else
   Azvet=Azvet2;
   Azvetf=Azvetf2;
   Gvet=Gvet2;
   alvet=alvet2;
  end
  ifps=ifp;
  if ifp>=-2
   ifp=-1;
  end
  isav=0;
  camv_new
  if ifp>=-3 | ifp==-10 | ifp==-4
   fsaf=figure,
   vg=0.5;
   puf=1:length(fso);
   subplot(211), plot(fou(puf,:)',aou(puf,:)',fso,fso*0,'wo'), grid;
   dista=(max(gou(puf,:))-min(gou(puf,:)))/min(gou(puf,:));
   title(' polar=2 ')
   if dista>10
    subplot(212), semilogy(fou(puf,:)',gou(puf,:)'/vg,fso,gso/vg,'wo');
   else
    subplot(212), plot(fou(puf,:)',gou(puf,:)'/vg,fso,gso/vg,'wo');
   end
  end

end

