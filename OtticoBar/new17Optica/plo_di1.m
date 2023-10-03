%load qua
%load ott
%load pro
%close all

iaoff=1;
enlar=0;
nsub2=1;
enlar=300;
nsub2=2;
desi=1;

if loop==1
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

end
PPlot=Plom;
Plot=PPlot{loop};
Fs=1;



nmod=length(PPlot{1}.parmod);
nmasce=nmod/2;
imm=1;
ifp=-10;
iEz=0;
iFF=0;
%for km=1:nmasce
for km=1:1
 for imo=1:2
  imod=(imo-1)*nmasce+km;
  Ap=abs(Plot.Ap{imod,imm});
  Ap=Ap/max(Ap);
  ApQ=abs(Plot.ApQ{imod,imm});
  ApQ=ApQ/max(ApQ);
  vpa=Plot.parmod{imod,imm};
  glosout=vpa(1);
  lares=vpa(2);
  ze=vpa(3);
  XP=Plot.XP{imod,imm};
  YP=Plot.YP{imod,imm};
  E2xo=Plot.E2xo{imod,imm};
  E2xp=Plot.E2xp{imod,imm};
  if iLP==0
   E2yo=Plot.E2yo{imod,imm};
   E2yp=Plot.E2yp{imod,imm};
  end
  Ef=Plot.Ef{imod,imm};
  Cug=Plot.Cug{imod,imm};
  aax=.8*max(max(XP));

  graph
  pausak
%  sub_cont

 end
%  sinc=sinc+1;
end







