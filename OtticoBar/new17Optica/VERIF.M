ipoi=0;
ieh=1;
iex=1;
isalf=0;    % 1 salta primi plots
ifita=1;
npadd=5;
iztm=1;
isce=0;


%FF=PPlot{1}.FF{1};
close all
iLP=   FF.iLP;
Mvefm0=FF.Mvefm0;
if iLP==0
 Mvefp0=FF.Mvefp0;
 Mvegm0=FF.Mvegm0;
 Mvegp0=FF.Mvegp0;
 Mvez0 =FF.Mvez0;
 Mvhz0 =FF.Mvhz0;
 besp=FF.besp;
 besz=FF.besz;
end
Pusas=FF.Pusas;
KKt=FF.KKt;
KKv=FF.KKv;
numodiacc=FF.numodiacc;
lambda=FF.lambda;
rr=FF.rr;
ze=FF.ze;
lambdam=lambda*(1+ze);
pimu=Ppol.pimu-1;
numo=length(pimu);
nudiv=numo-numodiacc;
dP=find(diff(KKt(Pusas))<0 | diff(Pusas')>1)+1;
dPai=[1; dP];
dPaf=[dP-1; length(Pusas)];
Psep=[dPai dPaf];

Nz=FF.Nz;
segem=FF.segem;
iredmat=FF.iredmat;
fmlstot=FF.fmlstot;
iauto=FF.iauto;
Litot=FF.Litot;
aitot=FF.aitot;
nitot=FF.nitot;
aiat=FF.aiat;
YP=FF.YP;
XP=FF.XP;
besm=FF.besm;
iFFsect=FF.iFFsect;
xvero=FF.xvero;
faX=1;
xvero=linspace(0,xvero(end)*faX,length(xvero));

Pusas=FF.Pusas;
nv=FF.nv;
Gas=FF.Gas;
Gad=FF.Gad;
Znor=FF.Znor;
Pus=   FF.Pus;

pes=FF.pes;
%pes=1;
Acoz=FF.Acoz;
zi=FF.zi;
if iredmat==1
 Pusff=Pusas;
else
 Pusff=Pus;
end

%keyboard

 puuu1=1:length(Gas);
 puuu2=length(Gas)+1:2*length(Gas);
 Liz0=1250/1000;
% Liz0=250/1000;

 Zsotto=5;
 Zma=80;
 dzc=5;
 offset=0;
 zcut=-([-(Zsotto+floor(zi(end)/dzc)*dzc):dzc:Zma]+offset);
 zcut=-([0:dzc:Zma]+offset);
% zcut=[0 3.5];

%' in Fz_tot.m ', keyboard

%  pvet=pola;
%  PUA=pua;



sM=size(Mvefm0);
sm2=sM(2);
if90=(sm2-1)/4+1;

 dz=Liz0;
 if iFFsect==1
  rff=nv(1,1);
  rffb=nv(end,1);
 else
  rffb=nv(1,1);
  rff=nv(end,1);
 end

 if iFFsect==1
  segdir=-1;
  sectb=[0:dz:Zsotto];
  sectb=sectb(2:end);
  sectu=[0:dz:Zma];
  sectu=-fliplr(sectu(1:end));
  ztot=[sectu zi sectb+zi(end)];
  Nztot=[ones(size(sectu))*rff Nz ones(size(sectb))*rffb];
 else
  segdir=1;
  sectb=[0:dz:Zsotto];
  sectb=-fliplr(sectb(2:end));
  sectu=[0:dz:Zma];
  sectu=sectu(2:end);
  ztot=[sectb zi zi(end)+sectu];
  Nztot=[ones(size(sectb))*rff Nz ones(size(sectu))*rffb];
 end

% zcut=sectu;
 fatEzv=(rr./Nz).^2;
 ztot=sectu;

 fatEz=fatEzv(1);

 kref=2*pi/lambdam*rr;
 dK=max(diff(KKv))/15;
 dK=max(diff(KKv))/1;



%sopra
 beta0=2*pi/lambdam*rff;
 sqr=sqrt(1-(KKt(Pusas)*rr/rff).^2);
 rrfa=rr/rff;
 lim=1-(rrfa*KKv).^2;
 kacce=find(lim<0);
 kaccea=[kacce; kacce+length(lim)];
 SQ=sqrt(lim);
 Znorn=[1./SQ; SQ];
 Ynorn=1./Znorn;
 Ynorn(kaccea)=0;
 Ynorn=Ynorn(Pusff);

 fkttn=KKt./[SQ; SQ];
 fkttn(kaccea)=0;
 fkttn=fkttn(Pusff);

 fatcut=ones(size(KKt));
% fatcut(kaccea)=0;
 fatcut=fatcut(Pusff);


 Are0=diag(pes)*(Acoz(puuu1,:)+Acoz(puuu2,:));
 Arh0=segdir*diag(pes)*(Acoz(puuu1,:)-Acoz(puuu2,:));

 Are=Are0;
 Arez=(-2*j)*(-segem)*fatEz*diag(fkttn)*Arh0;

 Arh=diag(Ynorn)*Arh0*diag(Nz);
 Arhz=(-2*j)*diag(Ynorn.*fkttn)*Are0*diag(Nz);

 Aue1=fatcut.*Are(:,1);
 Auez1=fatcut.*Arez(:,1);
 Auh1=fatcut.*Arh(:,1);
 Auhz1=fatcut.*Arhz(:,1);

% dK=max(diff(KKv));
 KKtdu=KKt(Pusas);

 Aue21=[];
 Aue22=[];
 Auez21=[];
 Auez22=[];
 Auh21=[];
 Auh22=[];
 Auhz21=[];
 Auhz22=[];

 Kfit=[];
 besp=[];
 besm=[];
 besz=[];
 Mvefm=[];
 Mvefp=[];
 Mvegm=[];
 Mvegp=[];
 Mvez =[];
 Mvhz =[];
 kin=0;
 kint=0;
 kref=2*pi/lambdam*rr;
 while kin<numo
  ialt=0;
  kin=kin+1;
  kint=kint+1;
  ibm=pimu(kin)-1;
  ibp=pimu(kin)+1;
  ibz=pimu(kin);
%  puk=Psep(kin,1):Psep(kin,2);
%  puk=Pusas(Psep(kin,:));
  puk=(Psep(kint,:));
  Ke=KKtdu(puk);
  if Ke(1)==min(KKv)
   Ke(1)=0;
  end
  Npn=ceil(diff(Ke)/dK)+1;
  Npo=diff(puk)+1;
%  keyboard
  if Npo<Npn
   Kad=linspace(Ke(1),Ke(2),Npn)';
   Np=Npn;
  else
   Kad=linspace(Ke(1),Ke(2),Npo)';
   Np=Npo;
  end

%  [Ke' Np], pausak
  Kfit=[Kfit; Kad];

   subfi

  if length(find(kin==nudiv))==1
   ialt=1;
   kint=kint+1;
%   puk=Pusas(Psep(kint,:));
%   Ke=KKv(puk);
   puk=(Psep(kint,:));
   Ke=KKtdu(puk);
%   [kint Ke' Np], pausak
   Np=ceil(diff(Ke)/dK)+1;
   Kad=linspace(Ke(1),Ke(2),Np)';
   Kfit=[Kfit; Kad];
   subfi
  end
 end


  besp=[besp; besp];
  besm=[besm; besm];
  besz=[besz; besz];

  Auh=[Auh21; Auh22];
  Aue=[Aue21; Aue22];
  Auhz=[Auhz21; Auhz22];
  Auez=[Auez21; Auez22];



 KKt=[Kfit;  Kfit];
 sqr=sqrt(1-(KKt*rr/rff).^2);

% keyboard




 clear Ae Aez  Ah Ahz
 for ks=1:length(sectu)
  Lu=sectu(ks);
  Adu=(Aue.*exp(-j*beta0*Lu*sqr));
  Ae(:,ks)=Adu;
  Adu=(Auh.*exp(-j*beta0*Lu*sqr));
  Ah(:,ks)=Adu;
  Adu=(Auez.*exp(-j*beta0*Lu*sqr));
  Aez(:,ks)=Adu;
  Adu=(Auhz.*exp(-j*beta0*Lu*sqr));
  Ahz(:,ks)=Adu;
 end
 Aeu=Ae;
 Ahu=Ah;
 Aezu=Aez;
 Ahzu=Ahz;


 if iFFsect==1
  Azue=[Aeu];
  Azuh=[Ahu ];
  Azuez=[Aezu ];
  Azuhz=[Ahzu ];
 else
  Azue=[  Aeu];
  Azuh=[  Ahu];
  Azuez=[ Aezu];
  Azuhz=[ Ahzu];
 end



  if iLP==1

   if exist('if90')
    bes_an=diag(Mvefm(:,if90))*besm;
   else
    bes_an=besm;
   end
   Esr=(bes_an'*Azu);
   Es=abs(Esr/max(max(Esr))).^2;

  else  %vettoriale

   Mfm=[ Mvefm;  segem*Mvefm];
   Mfp=[ Mvefp; -segem*Mvefp];

   if exist('if90')
    bes_amf=diag(Mfm(:,if90))*besm;
    bes_apf=diag(Mfp(:,if90))*besp;
   else
    if0=1;
    bes_amf=diag(Mfm(:,if0))*besm;
    bes_apf=diag(Mfp(:,if0))*besp;
   end


   Exz=Azue(:,end).'*(bes_apf+bes_amf);


   Aue=Azue(:,end);


   Mfme=diag(Aue)*[ Mvefm;  segem*Mvefm];
   Mfpe=diag(Aue)*[ Mvefp; -segem*Mvefp];


   Ext=besp'*Mfpe+besm'*Mfme;




  end

 figure, plot(xvero,abs(Exz),xvero,abs(Ext(:,if90)))
