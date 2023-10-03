close all
 co1(1:2,1)='.r';
 co1(1:2,2)='sg';
 co1(1:2,3)='oc';
 co1(1:2,4)='db';
 co1(1:2,5)='pm';
 co1(1:2,6)='hw';
clear Gav Wav

wa=figure;
hold on

iinte=0;
%ipo=input(' pol x/y [0/1]  ');
%iar=input(' array 4 o 6 [4/6]  ');
ipri=1;
 enlar=300;
 if iLP==0
  plar=850;
 else
  plar=550;
 end
 pograp=[ 40 100 280+enlar plar];
 pogram=[ 40 -50 280+enlar plar];
iaoff=1;
for loop=loopind

for imodi=1:2

Plot=PPlot{loop,imodi};
%nmod=length(PPlot{2}.parmod);
%nmasce=nmod/2;
imm=1;
ifp=-10;
if iLP==1
 iEz=0;
else
 iEz=1;
end
iEz=0
iFF=1;
iFFte=1;
sinc=1;
nsub2=2;

%nmasce=10

 Cug=Plot.Cug;
 aax=Ppol.a0ref;
 XP=Plot.XP;
 YP=Plot.YP;
 X=Plot.X;
 Y=Plot.Y;
% cce=Ppol.cce;
 ibar=1;
 if iLP==0
  ipolari=input(' vuoi polar [0/1] ');
 else
%  clear ipolari
  ipolari=0;
 end
du=size(Plot.parmod);
clear gdu fdu
for km=1:du(1)
 du=Plot.parmod{km};
 gdu(km,1)=du(1);
 fdu(km,1)=du(2);
end
dus=ones(size(gdu));
[fis,ifis]=sort(gdu);
mod=find(fis>=3*fis(1));
if length(mod)>1
 mop=1:mod(1)-1;
 moplo=ifis(mop);
else
 moplo=ifis;
end

parlop=spacv(loop)
figure(wa),
plot(1e3*fdu(moplo),gdu(moplo),co1(:,loop))
Wav(1:length(moplo),loop)=1e3*fdu(moplo);
Gav(1:length(moplo),loop)=gdu(moplo);

imm=1;
pos=20;
for imod=moplo'
  Ap=abs(Plot.Ap{imod,imm});
  Ap=Ap/max(Ap);
  ApQ=abs(Plot.ApQ{imod,imm});
  ApQ=ApQ/max(ApQ);
  vpa=Plot.parmod{imod,imm};
  glosout=vpa(1);
  lares=vpa(2);
  ze=vpa(3);
  pola=vpa(4);
  X=Plot.X{imod,imm};
  Y=Plot.Y{imod,imm};
  XP=Plot.XP{imod,imm};
  YP=Plot.YP{imod,imm};
  E2xo=Plot.E2xo{imod,imm};
  E2xp=Plot.E2xp{imod,imm};
  if iLP==0
   E2zo=Plot.E2zo{imod,imm};
   E2yo=Plot.E2yo{imod,imm};
   E2yp=Plot.E2yp{imod,imm};
  end
  Ef=Plot.Ef{imod,imm};
  Cug=Plot.Cug{imod,imm};
  aax=.8*max(max(XP));
  axli=2*.8*max(max(X));

  graph

end
end % imodi
end % loop

fi=find(Gav==0);
Gav(fi)=NaN;
Wav(fi)=NaN;
