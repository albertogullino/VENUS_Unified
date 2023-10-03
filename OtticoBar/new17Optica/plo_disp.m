%load qua
%load ott
%load pro
close all

enlar=0;
nsub2=1;
enlar=300;
nsub2=2;
desi=1;

if iLP==0
 if desi==1
  sinc=+1;
  pograp=[ 40 100 280+enlar 850];
  pogram=[ 40 -50 280+enlar 850];
  else
  sinc=-1;
  pograp=[ 710 100 280+enlar 850];
  pogram=[ 710 -50 280+enlar 850];
 end
else
 if desi==1
  sinc=+1;
  pograp=[150 300 600 650];
  pogram=[150 50 600 650];
  else
  sinc=-1;
  pograp=[650 300 600 650];
  pogram=[650 50 600 650];
 end
end

Plot=PPlot{loop};

Ppol=PPol{loop};
KK    = Ppol.KK    ;
KKt   = Ppol.KKt   ;
numodi= Ppol.numodi;
pasnu = Ppol.pasnu ;
lbv   = Ppol.lbv   ;
kcav0 = Ppol.kcav0 ;
x     = Ppol.x     ;
fian  = Ppol.fian  ;
Nx    = Ppol.Nx    ;
mbv   = Ppol.mbv   ;
%ibm   = Ppol.ibm   ;
%ibp   = Ppol.ibp   ;
ibm=[1:2:19];
ibp=ibm+2;
%Azvet1= Ppol.Az1;
%Azvet2= Ppol.Az2;
%Gvet1 = Ppol.G1 ;
%Gvet2 = Ppol.G2 ;


sinc=0;
nmod=length(PPlot{1}.parmod);
nmasce=nmod/2;
imm=1;
ifp=-10;
iEz=0;
iFF=0;
for km=1:nmasce
 for imo=1:2
  imod=(imo-1)*nmasce+km;
  Ap=Plot.Ap{imod,imm};
  ApQ=Plot.ApQ{imod,imm};
  vpa=Plot.parmod{imod,imm};
  glosout=vpa(1);
  lares=vpa(2);
  ze=vpa(3);
  XP=Plot.XP{imod,imm};
  Yp=Plot.YP{imod,imm};
  E2xo=Plot.E2xo{imod,imm};
  E2xp=Plot.E2xp{imod,imm};
  if iLP==0
   E2yo=Plot.E2yo{imod,imm};
   E2yp=Plot.E2yp{imod,imm};
  end
  Ef=Plot.Ef{imod,imm};
  Cug=Plot.Cug{imod,imm};
  aax=.8*max(max(XP));

%  graph
%  pausak
  sub_cont

 end
  sinc=sinc+1;
end







