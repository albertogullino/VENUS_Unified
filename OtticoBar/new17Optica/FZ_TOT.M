
%' in Fz_tot.m ', keyboard
%clear
close all
%load toto1
%load toto2
%load GRd45
nox=1.6;
imo=1;
FF=PPlot{1}.FF{imo};
%' FF ', keyboard
ipoi=0;
ieh=1;
iex=imo;
isalf=0;    % 1 salta primi plots
ifita=1;
npadd=5;
%load pro
iztm=1;

Nz=FF.Nz;
segem=FF.segem;
iredmat=FF.iredmat;
fmlstot=FF.fmlstot;
iauto=FF.iauto;
lambda=FF.lambda;
Litot=FF.Litot;
aitot=FF.aitot;
nitot=FF.nitot;
aiat=FF.aiat;
KKv=FF.KKv;
KKt=FF.KKt;
YP=FF.YP;
XP=FF.XP;
ze=FF.ze;
rr=FF.rr;
besm=FF.besm;
iFFsect=FF.iFFsect;
xvero=FF.xvero;
Pusas=FF.Pusas;
nv=FF.nv;
Gas=FF.Gas;
Gad=FF.Gad;
Znor=FF.Znor;
Pus=   FF.Pus;
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
 Liz0=250/1000;

 Zsotto=1;
 Zma=40;
 dzc=5;
 offset=0;
 zcut=-([-(floor(zi(end)/dzc)*dzc):dzc:Zma]+offset);
% zcut=[7:.05:10];


%  pvet=pola;
%  PUA=pua;



sM=size(Mvefm0);
sm2=sM(2);
if90=(sm2-1)/4+1;

 dz=Liz0;
 lambdam=lambda*(1+ze);
 if iFFsect==1
  rff=nv(1,1);
  rffb=real(nv(end,1));
 else
  rffb=nv(1,1);
  rff=nv(end,1);
 end

 if iFFsect==1
  segdir=-1;
  sectb=[0:dz:Zsotto];
  sectb=sectb(2:end);
  sectu=[0:dz:Zma];
  sectu=-fliplr(sectu(2:end));
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
 fatEzv=(rr./Nz).^2;

 fatEz=fatEzv(1);


%sopra
 beta0=2*pi/lambdam*rff;
 sqr=sqrt(1-(KKt(Pusas)*rr/rff).^2);
 rrfa=rr/rff;
 lim=1-(rrfa*KKv).^2;
 SQ=sqrt(lim);

 if iLP==0
  Znorn=[1./SQ; SQ];
  Ynorn=1./Znorn;
  Ynorn=Ynorn(Pusff);

  fkttn=KKt./[SQ; SQ];
  fkttn=fkttn(Pusff);
 else
  fkttn=KKt(Pusff);
  Ynorn=1;
 end

% fatcut=ones(size(KKt));
% fatcut=fatcut(Pusff);
fatcut=1;


 Are0=diag(pes)*(Acoz(puuu1,:)+Acoz(puuu2,:));
 Arh0=segdir*diag(pes)*(Acoz(puuu1,:)-Acoz(puuu2,:));

 Are=Are0;
 Arez=(-2*j)*(-segem)*fatEz*diag(fkttn)*Arh0;

 Arh=diag(Ynorn)*Arh0*diag(Nz);
 Arhz=(-2*j)*diag(Ynorn.*fkttn)*Are0*diag(Nz);

 Aue=fatcut.*Are(:,1);
 Auez=fatcut.*Arez(:,1);
 Auh=fatcut.*Arh(:,1);
 Auhz=fatcut.*Arhz(:,1);


 clear Ae Aez  Ah Ahz
 for ks=1:length(sectu)
  Lu=sectu(ks);
  fae=j*beta0*Lu*sqr;
  fifa=find(real(fae)>0);
  fae(fifa)=-fae(fifa);
  Axp=exp(fae);

  Adu=(Aue.*Axp);
  Ae(:,ks)=Adu;
  Adu=(Auh.*Axp);
  Ah(:,ks)=Adu;
  Adu=(Auez.*Axp);
  Aez(:,ks)=Adu;
  Adu=(Auhz.*Axp);
  Ahz(:,ks)=Adu;
 end
 Aeu=Ae;
 Ahu=Ah;
 Aezu=Aez;
 Ahzu=Ahz;

% buffer
 fatEz=fatEzv(end);
 beta0=2*pi/lambdam*rffb;
 sqr=sqrt(1-(KKt(Pusas)*rr/rffb).^2);
 lim=1-(KKv*rr/rffb).^2;

 if iLP==0
  SQ=sqrt(lim);
  Znorn=[1./SQ; SQ];
  Ynorn=1./Znorn;
  Ynorn=Ynorn(Pusff);

  fkttn=KKt./[SQ; SQ];
  fkttn=fkttn(Pusff);
 else
  fkttn=KKt(Pusff);
  Ynorn=1;
 end

% fatcut=ones(size(KKt));
% fatcut=fatcut(Pusff);


 Are=diag(fatcut)*Are0;
 Arez=-segem*2*(-j)*fatEz*diag(fatcut.*fkttn)*Arh0;

 Arh=diag(fatcut.*Ynorn)*Arh0;
 Arhz=-2*j*diag(Ynorn.*fkttn)*Are0;

 Aue=Are(:,end);
 Auez=Arez(:,end);
 Auh=Arh(:,end);
 Auhz=Arhz(:,end);


 clear Ae Aez  Ah Ahz
 for ks=1:length(sectb)
  Lu=sectb(ks);
  fae=-j*beta0*Lu*sqr;
  fifa=find(real(fae)>0);
  fae(fifa)=-fae(fifa);
  Axp=exp(fae);

  Adu=(Aue.*Axp);
  Ae(:,ks)=Adu;
  Adu=(Auh.*Axp);
  Ah(:,ks)=Adu;
  Adu=(Auez.*Axp);
  Aez(:,ks)=Adu;
  Adu=(Auhz.*Axp);
  Ahz(:,ks)=Adu;
 end
 Aeb=Ae;
 Ahb=Ah;
 Aezb=Aez;
 Ahzb=Ahz;

if iLP==0
 if iFFsect==1
  Azue=[Aeu Are Aeb];
  Azuh=[Ahu Arh Ahb];
  Azuez=[Aezu Are Aezb];
  Azuhz=[Ahzu Arh Ahzb];
 else
  Azue=[Aeb Are Aeu];
  Azuh=[Ahb Arh Ahu];
  Azuez=[Aezb Are Aezu];
  Azuhz=[Ahzb Arh Ahzu];
 end
else
 if iFFsect==1
  Azu=[Aeu Are Aeb];
 else
  Azu=[Aeb Are Aeu];
 end
 Azue=Azu;
end



  if iLP==1

   if exist('if90')
    bes_an=diag(Mvefm0(:,if90))*besm;
   else
    bes_an=besm;
   end
   Esr=(bes_an'*Azu);
   Es=abs(Esr/max(max(Esr))).^2;
   tit=' LP ';

  else  %vettoriale

   Mfm=[ Mvefm0;  segem*Mvefm0];
   Mgm=[ Mvegm0;  segem*Mvegm0];
   Mfp=[ Mvefp0; -segem*Mvefp0];
   Mgp=[-Mvegp0;  segem*Mvegp0];
   Mez =[ 0*Mvez0;         Mvez0];
   Mhz =[  Mvhz0;        0*Mvhz0];

   if exist('if90')
    bes_amf=diag(Mfm(:,if90))*besm;
    bes_amg=diag(Mgm(:,if90))*besm;
    bes_apf=diag(Mfp(:,if90))*besp;
    bes_apg=diag(Mfp(:,if90))*besp;
    bes_aez=diag(Mez(:,if90))*besz;
    bes_ahz=diag(Mhz(:,if90))*besz;
   else
    if0=1;
    bes_amf=diag(Mfm(:,if0))*besm;
    bes_amg=diag(Mgm(:,if0))*besm;
    bes_apf=diag(Mfp(:,if0))*besp;
    bes_apg=diag(Mfp(:,if0))*besp;
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

    Edu=(Edu/max(max(Edu)));
   end
   Es=abs(Edu.');
  end


 fiz=1:length(ztot);
 [du,pdu1]=min(abs(xvero-aiat/2));
 [du,pdu2]=min(abs(xvero-aiat));

 fiz0=find(ztot==0);
 fiz0=fiz0(1);
% figure
 if exist('soloplot')
  if soloplot==1
   return
  end
 end

if isalf==0
 isez=[1 pdu1 pdu2];
% isez=1;
 figure, semilogy(ztot(fiz),Es(isez,fiz)),
 pausak
 figure, plot(ztot(fiz),Es(isez,fiz)),
 xlabel(' zeta toward device (micron)')
 ylabel(' Es ')
 title([tit,' in dis-fz: 3 radial coord: 0, Ox/2, Ox '])
 pausak
 Nor=max(abs(Es(1,:)));

isa=1;
if isa==0
 Ex=Ex.';
 Ey=Ey.';
 Ez=Ez.';
 Hx=Hx.';
 Hy=Hy.';
 Hz=Hz.';
 fiz=[187 188];
 fiz=[267-40:268+40];
 figure, plot(ztot(fiz),real(Ex(isez,fiz)),ztot(fiz),imag(Ex(isez,fiz))),
 title(' Real, imag Ex')
 pausak
 figure, plot(ztot(fiz),real(Hy(isez,fiz)),ztot(fiz),imag(Hy(isez,fiz))),
 title(' Real, imag Hy')
 pausak
 figure, plot(ztot(fiz),real(Ex(isez,fiz)),ztot(fiz),real(Hy(isez,fiz))),
 title(' Real ')
 pausak
 figure, plot(ztot(fiz),imag(Ex(isez,fiz)),ztot(fiz),imag(Hy(isez,fiz))),
 title(' Imag ')
 pausak
 figure, plot(ztot(fiz),real(Ex(isez,fiz).*conj(Hy(isez,fiz)))),
% figure, plot(ztot(fiz),imag(Ex(isez,fiz).*conj(Hy(isez,fiz)))),
 title(' Ex*Hy')
 pausak
end

% map(Es(:,fiz)*diag(1./max(Es(:,fiz))),ztot(fiz),xvero)
  map(log10(abs(Es(:,fiz))),ztot(fiz),xvero)

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
 fiox=find(nitot(:,2)==nox & imag(nitot(:,2))==0 & nitot(:,2)>1);
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
end  % isalf

% fiz0=find(ztot==0);
% figure
% plot(xvero,Es(:,fiz0)*diag(1./max(Es(:,fiz0))))
% pausak
%' dis_fzt', keyboard

 fatEzv=(rr./Nztot).^2;

 h=figure;
% assi=linspace(5,27,9)
 assi=15;
 for k=1:length(zcut)
%  subplot(3,3,k)
  [du,ip]=min(abs(zcut(k)-ztot));

 fatEz=fatEzv(ip);
 rfip=Nztot(ip);
 beta0=2*pi/lambdam*rfip;
 sqr=sqrt(1-(KKt(Pusas)*rr/rfip).^2);
 lim=1-(KKv*rr/rfip).^2;

 if iLP==0
  SQ=sqrt(lim);
  Znorn=[1./SQ; SQ]./Nztot(ip);
  Ynorn=1./Znorn;
  Ynorn=Ynorn(Pusff);

  fkttn=KKt./[SQ; SQ];
  fkttn=fkttn(Pusff);
 else
  fkttn=KKt(Pusff);
  Ynorn=Nztot(ip);
 end



   Aue=Azue(:,ip);

  if iLP==1

   Mfme=diag(Aue)*Mvefm0;
   Esr=besm'*Mfme;
%   Edu=abs(Esr/max(max(Esr))).^2;
   Edu=abs(Esr).^2/Nor;

  else  %vettoriale
   Auez=Azuez(:,ip);
   Auh=Azuh(:,ip);
   Auhz=Azuhz(:,ip);

   Mfme=diag(Aue)*[ Mvefm0;  segem*Mvefm0];
   Mgme=diag(Aue)*[ Mvegm0;  segem*Mvegm0];
   Mfpe=diag(Aue)*[ Mvefp0; -segem*Mvefp0];
   Mgpe=diag(Aue)*[-Mvegp0;  segem*Mvegp0];

   Mfmh=diag(Auh)*[ Mvefm0;  segem*Mvefm0];
   Mgmh=diag(Auh)*[ Mvegm0;  segem*Mvegm0];
   Mfph=diag(Auh)*[ Mvefp0; -segem*Mvefp0];
   Mgph=diag(Auh)*[-Mvegp0;  segem*Mvegp0];

   Mez =diag(Auez)*[ 0*Mvez0;         Mvez0];
   Mhz =diag(Auhz)*[  Mvhz0;        0*Mvhz0];

%   Mvegp=diag(Az)*Mgp;
%   Ex=besp'*Mvefp+besm'*Mvefm;
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

%    Edu=(Edu/max(max(Edu)));
    Edu=(Edu/sqrt(Nor));
   end
%
  end

  Esez{k}=Edu;
%  map(log10(abs(Edu)),YP,XP,h,-1)
  map((abs(Edu)),YP,XP,h,-1)
%  axis([-1 1 -1 1]*assi(k))
  title(' real part' )
  axis([-1 1 -1 1]*assi)
  axis equal
% if ipoi==0 & length(zcut)<4
%  pausak
%  map(imag(Edu),YP,XP)
%  title(' imag part' )
%  axis([-1 1 -1 1]*assi)
%  axis equal
%  pausak
% end
  if k==1
  title([' Section = ',num2str(-zcut(k)),' micron'])
  else
  title([num2str(-zcut(k))])
  end
  pausak
 end

