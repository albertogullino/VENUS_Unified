 Zsotto=30;
 Zma=40;
dzc=10;
offset=2;
offset
pausak
assi=6;     % 0 se automatico
iFFexc=0;   % cambia lato far_field
fadK=1;
faX=1;
ipoi=0;
ieh=1;
iex=1;
isalf=0;    % 1 salta primi plots
ifita=1;
npadd=5;
iztm=1;
isce=0;

lP=input(' par = ');

for lp=lP
FF=PPlot{lp}.FF{1};
%FF.numodiacc=8;
%close all
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
%figure, plot(dPai,'g.'), hold on, plot(dPaf,'r.'), plot(Psep,'o'),
%pausak
if iLP==0
 puP=1:length(Pusas)/2;
 puS=1:length(Psep)/2;
else
 puP=1:length(Pusas);
 puS=1:length(Psep);
end
%figure, plot(KKv(Pusas(puP)),'.'), hold on,
%plot(Psep(puS,:),KKv(Pusas(Psep(puS,:))),'o')
%pausak

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
xvero=xvero*faX;

Pusas=FF.Pusas;
nv=FF.nv;
Gas=FF.Gas;
Gad=FF.Gad;
Znor=FF.Znor;
Pus=   FF.Pus;
if isfield(FF,'ikr')==0
 ikr=2;
else
 ikr=FF.ikr;
end

pes=FF.pes;
%pes=1;
Acoz=FF.Acoz;
zi=FF.zi;
fatEzv=(rr./Nz).^2;

if iredmat==1
 Pusff=Pusas;
else
 Pusff=Pus;
end

%keyboard

 puuu1=1:length(Gas);
 puuu2=length(Gas)+1:2*length(Gas);
 Liz0=1250/1000;
 Liz0=250/1000;
 Liz0=500/1000;

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
 if iFFexc==1
  if iFFsect==1
   iFFsect=0;
  else
   iFFsect=1;
  end
 end

 if iFFsect==1
  rff=nv(1,1);
  rffb=nv(end,1);
  isezE=1;
 else
  rffb=nv(1,1);
  rff=nv(end,1);
  zcut=-zcut;
  isezE=length(fatEzv);
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
 ztot=sectu;

 fatEz=fatEzv(isezE);

 kref=2*pi/lambdam*rr;
% dK=max(diff(KKv))/15;
 dK=max(diff(KKv))/fadK;



%sopra
 beta0=2*pi/lambdam*rff;
 sqr=sqrt(1-(KKt(Pusas)*rr/rff).^2);
 rrfa=rr/rff;
 if iLP==0
  lim=1-(rrfa*KKv).^2;
  SQ=sqrt(lim);
  Znorn=[1./SQ; SQ];
  Ynorn=1./Znorn;
  Ynorn=Ynorn(Pusff);

  fkttn=KKt./[SQ; SQ];
  fkttn=fkttn(Pusff);
 else
  Ynorn=1;
  fkttn=KKt(Pusff);
 end


 Are0=diag(pes)*(Acoz(puuu1,:)+Acoz(puuu2,:));
 Arh0=segdir*diag(pes)*(Acoz(puuu1,:)-Acoz(puuu2,:));

 Are=Are0;
 Arez=(-2*j)*(-segem)*fatEz*diag(fkttn)*Arh0;

 Arh=diag(Ynorn)*Arh0*diag(Nz);
 Arhz=(-2*j)*diag(Ynorn.*fkttn)*Are0*diag(Nz);
 fatcut=1;

 Aue1=fatcut.*Are(:,isezE);
 Auez1=fatcut.*Arez(:,isezE);
 Auh1=fatcut.*Arh(:,isezE);
 Auhz1=fatcut.*Arhz(:,isezE);

% keyboard

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
    if iLP==0
     ibm=pimu(kin)-1;
     ibp=pimu(kin)+1;
     ibz=pimu(kin);
    else
     ibm=pimu(kin);
    end
    puk=(Psep(kint,:));
    Ke=KKtdu(puk);
    if Ke(1)==min(KKv)
     Ke(1)=0;
    end
    Npn=ceil(diff(Ke)/dK)+1;
    Npo=diff(puk)+1;
    if Npo<Npn
     Kad=linspace(Ke(1),Ke(2),Npn)';
     Np=Npn;
    else
     Kad=linspace(Ke(1),Ke(2),Npo)';
     Np=Npo;
    end

    Kfit=[Kfit; Kad];

     subfi

    if length(find(kin==nudiv))==1 & numo>1 & ikr==2
     ialt=1;
     kint=kint+1;
     puk=(Psep(kint,:));
     Ke=KKtdu(puk);
     Np=ceil(diff(Ke)/dK)+1;
     Kad=linspace(Ke(1),Ke(2),Np)';
     Kfit=[Kfit; Kad];
     subfi
    end
   end

% figure, subplot(211), plot(Kfit,'g.'),
% title(' new K ')
% subplot(212), plot(KKt(Pusas(puP)),'r.'),
% title(' old K ')
% pausak

  if iLP==0
   besp=[besp; besp];
   besz=[besz; besz];
   besm=[besm; besm];

   Aue=[Aue21; Aue22];
   Auh=[Auh21; Auh22];
   Auhz=[Auhz21; Auhz22];
   Auez=[Auez21; Auez22];
   KKt=[Kfit;  Kfit];
  else
   Aue=Aue21;
   KKt=Kfit;
  end



 sqr=segdir*sqrt(1-(KKt*rr/rff).^2);

% keyboard




 clear Ae Aez  Ah Ahz
 for ks=1:length(sectu)
  Lu=sectu(ks);
  Adu=(Aue.*exp(-j*beta0*Lu*sqr));
  Ae(:,ks)=Adu;
  if iLP==0
   Adu=(Auh.*exp(-j*beta0*Lu*sqr));
   Ah(:,ks)=Adu;
   Adu=(Auez.*exp(-j*beta0*Lu*sqr));
   Aez(:,ks)=Adu;
   Adu=(Auhz.*exp(-j*beta0*Lu*sqr));
   Ahz(:,ks)=Adu;
  end
 end
 Aeu=Ae;
 if iLP==0
  Ahu=Ah;
  Aezu=Aez;
  Ahzu=Ahz;
 end


 if iFFsect==1
  if iLP==0
   Azue=[Aeu];
   Azuh=[Ahu ];
   Azuez=[Aezu ];
   Azuhz=[Ahzu ];
  else
   Azu=[Aeu];
  end
 else
  if iLP==0
   Azue=[  Aeu];
   Azuh=[  Ahu];
   Azuez=[ Aezu];
   Azuhz=[ Ahzu];
  else
   Azu=[  Aeu];
  end
 end



  if iLP==1

   if exist('if90')
    bes_an=diag(Mvefm(:,if90))*besm;
   else
    bes_an=besm;
   end
   Esr=(bes_an'*Azu);
   Edu=abs(Esr.'/max(max(Esr))).^2;
   tit=' LP ';

  else  %vettoriale

   Mfm=[ Mvefm;  segem*Mvefm];
   Mgm=[ Mvegm;  segem*Mvegm];
   Mfp=[ Mvefp; -segem*Mvefp];
   Mgp=[-Mvegp;  segem*Mvegp];
   Mez =[ 0*Mvez;         Mvez];
   Mhz =[   Mvhz;       0*Mvhz];

   if exist('if90')
    bes_amf=diag(Mfm(:,if90))*besm;
    bes_amg=diag(Mgm(:,if90))*besm;
    bes_apf=diag(Mfp(:,if90))*besp;
    bes_apg=diag(Mgp(:,if90))*besp;
    bes_aez=diag(Mez(:,if90))*besz;
    bes_ahz=diag(Mhz(:,if90))*besz;
   else
    if0=1;
    bes_amf=diag(Mfm(:,if0))*besm;
    bes_amg=diag(Mgm(:,if0))*besm;
    bes_apf=diag(Mfp(:,if0))*besp;
    bes_apg=diag(Mgp(:,if0))*besp;
    bes_aez=diag(Mez(:,if0))*besz;
    bes_ahz=diag(Mhz(:,if0))*besz;
   end


   Ex=Azue.'*(bes_apf+bes_amf);
   Ey=Azue.'*(bes_apg+bes_amg);
   Hy=Azuh.'*(bes_apf+bes_amf);
   Hx=-Azuh.'*(bes_apg+bes_amg);

   if iztm==1
    Ez=Azuez.'*bes_aez;
    Hz=Azuhz.'*bes_ahz;
   else
    Ez=0;
    Hz=0;
   end

   Pvectz=(Ex.*conj(Hy)-Ey.*conj(Hx));
   Pvectx=(Ey.*conj(Hz)-Ez.*conj(Hy));
   Pvecty=-((Ex.*conj(Hz)-Ez.*conj(Hx)));
   Pointingr=sqrt(real(Pvectx).^2+real(Pvecty).^2+real(Pvectz).^2);
   Pointingi=sqrt(imag(Pvectx).^2+imag(Pvecty).^2+imag(Pvectz).^2);

%   Es=Pointing';
%   Es=Es'/max(max(Es));
   if ipoi>=1
    if ipoi==2
     Edu=Pointingi;
    else
     Edu=Pointingr;
    end
%    Edu=Pvectz;
      tit= 'Pointing ';
   else
    if iex==1
     if ieh==1
      Edu=Ex;
      tit= 'Ex ';
     else
      Edu=Hx;
      tit= 'Hx ';
     end
    elseif iex==2
     if ieh==1
      Edu=Ey;
      tit= 'Ey ';
     else
      Edu=Hy;
      tit= 'Hy ';
     end
    elseif iex==3
     if ieh==1
      Edu=Ez;
      tit= 'Ez ';
     else
      Edu=Hz;
      tit= 'Hz ';
     end
    end

    Edu=abs(Edu/max(max(Edu))).^2;
   end
  end
  Es=Edu.';


 fiz=1:length(ztot);
 [du,pdu1]=min(abs(xvero-aiat/2));
 [du,pdu2]=min(abs(xvero-aiat));

% fiz0=find(ztot==0);
% fiz0=fiz0(1);
% figure
 if exist('soloplot')
  if soloplot==1
   return
  end
 end

if isalf==0
 isez=[1 pdu1 pdu2];
% isez=1;
%' wui', keyboard
 figure, plot(ztot(fiz),Es(isez,fiz)),
% figure, semilogy(ztot(fiz),Es(isez,fiz)),
 xlabel(' zeta toward device (micron)')
 ylabel(' Es ')
 title([tit,' in dis-fz: 3 radial coord: 0, Ox/2, Ox '])
 pausak

 if length(lP)==1

% map(Es(:,fiz)*diag(1./max(Es(:,fiz))),ztot(fiz),xvero)
%  map(log10(abs(Es(:,fiz))),ztot(fiz),xvero)
  map((abs(Es(:,fiz))).^2,ztot(fiz),xvero)

 axis([ztot((1)),ztot((end)),0,max(xvero)])


 hold on,
 if aitot(1,1)~=0 & nitot(1,2)~=0
  xi1=xvero;
  fif=find(xi1<aitot(1,1));
  xi1(fif)=NaN;
 else
  xi1=xvero;
 end
 xi0=xvero;
 zi1=ones(size(xvero))*0;
 zi0=ones(size(xvero));
 firef=find(iauto(:,1)==2)-1;
 firef=[firef find(iauto(:,1)==3)-2];
 fiox=find(nitot(:,2)<2 & imag(nitot(:,2))==0 & nitot(:,2)>1);
 firef=[firef fiox];
 for k=firef
  xdu=xi0;
  if aitot(k,1)~=0
   fif=find(xdu<aitot(k,1));
   xdu(fif)=NaN;
  end
  xi1=[xi1; xdu];
  fmul=fmlstot(1:k,2);
  fif=find(fmul==0);
  fmul(fif)=1;
  zi1=[zi1; zi0*sum(Litot(1:k).*fmul)];
 end
 hi=ones(size(xi1));
 hl=plot3(zi1',xi1',hi','w');
 set(hl,'linewidth',2)
 pausak
 end
end  % isalf

 if length(lP)==1
%  map(log10(abs(Es*diag(1./Es(1,:)))),ztot,xvero)
%  map(log10(abs(Es)),ztot,xvero)
%  title(' LOG scale ')
% axis([ztot((1)),ztot(end),0,max(xvero)])
% pausak

 fatEzv=(rr./Nztot).^2;


 h=figure;
 if assi==0
 assi=xvero(end);
 end
 for k=1:length(zcut)
  [du,ip]=min(abs(zcut(k)-ztot));



  if iLP==1

   Aue=Azu(:,ip);
   Mfme=diag(Aue)*Mvefm;
   Esr=besm'*Mfme;
   Edu=abs(Esr/max(max(Esr))).^2;

  else  %vettoriale

   Aue=Azue(:,ip);
   Auez=Azuez(:,ip);
   Auh=Azuh(:,ip);
   Auhz=Azuhz(:,ip);

   Mfme=diag(Aue)*[ Mvefm;  segem*Mvefm];
   Mgme=diag(Aue)*[ Mvegm;  segem*Mvegm];
   Mfpe=diag(Aue)*[ Mvefp; -segem*Mvefp];
   Mgpe=diag(Aue)*[-Mvegp;  segem*Mvegp];

   Mfmh=diag(Auh)*[ Mvefm;  segem*Mvefm];
   Mgmh=diag(Auh)*[ Mvegm;  segem*Mvegm];
   Mfph=diag(Auh)*[ Mvefp; -segem*Mvefp];
   Mgph=diag(Auh)*[-Mvegp;  segem*Mvegp];

   Mez =diag(Auez)*[ 0*Mvez;         Mvez];
   Mhz =diag(Auhz)*[  Mvhz;        0*Mvhz];

   Ex=besp'*Mfpe+besm'*Mfme;
   Ey=besp'*Mgpe+besm'*Mgme;

   Hy=besp'*Mfph+besm'*Mfmh;
   Hx=-(besp'*Mgph+besm'*Mgmh);


   if iztm==1
    Ez=besz'*Mez;
    Hz=besz'*Mhz;
   else
    Ez=0;
    Hz=0;
   end
   Pvectz=(Ex.*conj(Hy)-Ey.*conj(Hx));
   Pvectx=(Ey.*conj(Hz)-Ez.*conj(Hy));
   Pvecty=-((Ex.*conj(Hz)-Ez.*conj(Hx)));
   Pointingr=sqrt(real(Pvectx).^2+real(Pvecty).^2+real(Pvectz).^2);
   Pointingi=sqrt(imag(Pvectx).^2+imag(Pvecty).^2+imag(Pvectz).^2);

   if ipoi>=1
    if ipoi==2
     Edu=Pointingi;
    else
     Edu=Pointingr;
    end

   else
    if iex==1
     if ieh==1
      Edu=Ex;
     else
      Edu=Hx;
     end
    elseif iex==2
     if ieh==1
      Edu=Ey;
     else
      Edu=Hy;
     end
    elseif iex==3
     if ieh==1
      Edu=Ez;
     else
      Edu=Hz;
     end
    end

    Edu=(Edu/max(max(Edu))).^2;
   end
%
  end

  Esez{k}=Edu;
%  map(log10(abs(Edu)),YP*faX,XP*faX,h,-1)
  map((abs(Edu)),YP*faX,XP*faX,h,-1)

%  axis([-1 1 -1 1]*assi(k))
  title(' real part' )
  axis([-1 1 -1 1]*assi)
  axis equal
 if ipoi==0 & length(zcut)<4
  pausak
  map(imag(Edu),YP,XP)
  title(' imag part' )
%  axis([-1 1 -1 1]*assi(k))
  axis([-1 1 -1 1]*assi)
  axis equal
  pausak
 end
 colorbar off
% colormap gray
  if k==1
  title([' Section = ',num2str(-zcut(k)),' micron'])
  else
  title([num2str(-zcut(k))])
  end
  pausak
%  keyboard
 end
 end %lP==1
end
