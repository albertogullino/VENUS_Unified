close all
iinte=0;
%ipo=input(' pol x/y [0/1]  ');
%iar=input(' array 4 o 6 [4/6]  ');
 enlar=300;
 pograp=[ 40 100 280+enlar 850];
 pogram=[ 40 -50 280+enlar 850];
iaoff=1;
iar=narv(loop);
alma=iar;
palma=fix(alma/4);

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
Azvet1= Ppol.Az1;
Azvet2= Ppol.Az2;
Gvet1 = Ppol.G1 ;
Gvet2 = Ppol.G2 ;
alvet1 = Ppol.A1 ;
alvet2 = Ppol.A2 ;

Plot=PPlot{loop};
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

 Cug=Ppol.Cug;
 aax=Ppol.aax;
 XP=Ppol.XP;
 YP=Ppol.YP;
 X=Plot.X;
 Y=Plot.Y;
 cce=Ppol.cce;
 ibar=1;
 if iLP==0
  ipolari=input(' vuoi polar [0/1] ');
 else
  ipolari=0;
 end
du=size(Plot.parmod);
for km=1:du(1)
 du=Plot.parmod{km};
 gdu(km,1)=du(1);
end
dus=ones(size(gdu));
[fis,ifis]=sort(gdu);
mod=find(fis>=3*fis(1));
mop=mod(1:mod(1)-1);
moplo=ifis(mop)

imm=1;
pos=20;
for pola=[-1 1]
%for pola=[1 -1]
 if pola==-1
  ipo=1;
 else
  ipo=2;
 end
 for imo=1:nmasce
% for imo=1:3
% for imo=[1 23]
  imod=(ipo-1)*nmasce+imo;
  Ap=abs(Plot.Ap{imod,imm});
  Ap=Ap/max(Ap);
  ApQ=abs(Plot.ApQ{imod,imm});
  ApQ=ApQ/max(ApQ);
  vpa=Plot.parmod{imod,imm};
  glosout=vpa(1);
  lares=vpa(2);
  ze=vpa(3);
  X=Plot.X{imod,imm};
  Y=Plot.Y{imod,imm};
%  XP=Plot.XP{imod,imm};
%  Yp=Plot.YP{imod,imm};
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


 if pola==-1
  Exdu=Ppol.Ex1;
  Eydu=Ppol.Ey1;
  poy=250;
 else
  Exdu=Ppol.Ex2;
  Eydu=Ppol.Ey2;
%  Exdu1=Ppol.Ex2;
%  Eydu1=Ppol.Ey2;
  poy=50;
 end

  Exdu=E2xo;
  if iLP==0
   Eydu=E2yo;
  end

% MY1=max(max(Eydu1));
% MX1=max(max(Exdu1));
% if abs(MX1)>abs(MY1)
%  N1=MX1;
%  cht=' Ex ';
%  Estro1=Exdu1/N;
%  Esmal1=Eydu1/N;
% else
%  N1=MY1;
%  cht=' Ey ';
%  Estro1=Eydu1/N1;
%  Esmal1=Exdu1/N1;
% end
% figure
% map_fnew(XP,YP,imag(Estro),aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
% figure
% map_fnew(XP,YP,imag(Estro1),aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
% figure
% map_fnew(XP,YP,real(Esmal+Esmal1),aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)


 if iLP==0
 MY=max(max(Eydu));
 MX=max(max(Exdu));
 if abs(MX)>abs(MY)
  N=MX;
  cht=' Ex ';
  Estro=Exdu/N;
  Esmal=Eydu/N;
 else
  N=MY;
  cht=' Ey ';
  Estro=Eydu/N;
  Esmal=Exdu/N;
 end

 Ptot1=abs(Estro).^2;
 Ptot2=abs(Esmal).^2;
 Ptot=Ptot1+Ptot2;
 if iLP==0
  Ptot=abs(E2xo.^2+E2yo.^2);
 else
  Ptot=abs(E2xo.^2);
 end
 Ptot=Ptot/max(max(Ptot));
 if iinte==1
 figure
 pograpd=[700   500   400   400];
 set(gcf,'Position',pograpd)
 titl=['Total power: dominant field ',cht];
 map_fnew(XP,YP,Ptot,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)


  pausak
 end
 end  %iLP
 if ipolari==1
  polari
 end %ipolari

pausak
end  %ipo
end  %ipo
return




 figure
 Ptot=abs(Estro).^2+abs(Estro).^2;alfav=linspace(0,2,11);
 titl='Total power';supo0=30;
 map_fnew(XP,YP,Ptot,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)

for ia=1:length(alfav)
 alfa=alfav(ia);
 supo=supo0*ia;
 salfa=sin(alfa*pi/180);
 calfa=cos(alfa*pi/180);
 Pup=abs(Esmal*salfa).^2+abs(calfa*Estro).^2+2*salfa*calfa*real(Esmal.*conj(Estro));

 alfa=-alfav(ia);
 salfa=sin(alfa*pi/180);
 calfa=cos(alfa*pi/180);
 Pum=abs(Esmal*salfa).^2+abs(calfa*Estro).^2+2*salfa*calfa*real(Esmal.*conj(Estro));

 figure
 pograpd=[0+supo   poy   550   700];
 set(gcf,'Position',pograpd)
 subplot(2,1,1)
 titlm=num2str(alfa);
 titlp=num2str(-alfa);
 map_fnew(XP,YP,Pum,aax,Cug.x,Cug.y,Cug.z,titlm,ibar,iaoff)

 subplot(2,1,2)
 map_fnew(XP,YP,Pup,aax,Cug.x,Cug.y,Cug.z,titlp,ibar,iaoff)
% pausak

end


% s1=abs(Ex1.^2)-abs(Ey2.^2);
% s2=-(conj(Ex1).*Ey2+conj(Ey2).*Ex1);
% s3=j*(conj(Ex1).*Ey2-conj(Ey2).*Ex1);


firat=Esmal.*conj(Estro)/abs(N^2);

Xi=imag(firat);
Phi=real(firat);
%fiI=find(isinf(abs(Xi)));
%Xi(fiI)=0;
%fiI=find(isinf(abs(Phi)));
%Phi(fiI)=0;
%
%fiI=find(abs(Xi)>.05);
%Xi(fiI)=0;
%fiI=find(abs(Phi)>.05);
%Phi(fiI)=0;

 subplot(2,3,3)
 titl=' Phi';
 map_fnew(XP,YP,Phi,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)

 titl='Xi';
 subplot(2,3,6)
 map_fnew(XP,YP,Xi,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)




 titl=' Alpha';
 figure
 map_fnew(XP,YP,Phi*180/pi,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)

 figure
 pograpd=[300   200   550   700];
 set(gcf,'Position',pograpd)
 subplot(2,1,1)
 titl=' + ';
 map_fnew(XP,YP,Xi.^2+(Phi-max(max(Phi))).^2,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)

 subplot(2,1,2)
 titl=' - ';
 map_fnew(XP,YP,Xi.^2+(Phi-min(min(Phi))).^2,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)

%        figure
%        map_fnew(XP,YP,E2xo,aax,Cug.x,Cug.y,Cug.z,'',0,1)
%        axis([-1 1 -1 1]*10)
%        colormap gray
%        figure
%        map_fnew(XP,YP,E2xp,aax,Cug.x,Cug.y,Cug.z,'',0,1)
%        axis([-1 1 -1 1]*10)
%        colormap gray

