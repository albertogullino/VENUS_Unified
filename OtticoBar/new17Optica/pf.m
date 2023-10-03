

if idyn==1
   sE=size(Eqw);
   if i2D==0
    modcv=[1:sE(2)];
   else
    if length(sE)==3
     modcv=[1:sE(3)];
    else
     modcv=1;
    end
   end
%   modc0=[1:sE(end)];
    modc=modcv;
   if modc0~=0
    if length(modcv)>length(modc0)
     modc=modc0;
    end
   end
%    modcu=modcv;
    modcu=modc;

%' modc ', keyboard



%   Gam_u =Pa.Gam_u;
%   Wvmo =Pa.W;
%   fa4mo =Pa.famod;
%   czmo =Pa.cz;
%   fqwmo=Pa.fqw;
%   Lfmo =Pa.Lf;
%   Tfmo =Pa.Tf;
%   Tfimo=Pa.Tfinf;
%   gtmo =Pa.gt;
%   fatqw=Pa.fatqw;
%   d_at =Pa.dat;
%   NQW  =Pa.NQW;
%   rr   =Pa.rr;
%   avero=Pa.avero;

   Tefv=Pa.trasm;


   rr   =Pa.rr;
   rvg   =Pa.rmed;
%   fqwmo=Pa.fqw;
   d_at =Pa.dat*1e-6;    % in metri
   NQW=Pa.NQW;
   dat_tot=d_at*NQW;
   tyPmodu=Pa.type;
%   'calopt rr', keyboard

   tyPmod=tyPmodu(modc);

   taut=Pa.taut;
   tauu=Pa.tauu;
   taub=Pa.taub;
%   alca=Pa.alca;
%   pvol=Pa.pvol;

   confztot=Pa.confztot(2);
   confz=Pa.confztot(1);

%   losm=Pa.losm/confz*confztot/NQW;
   losm=Pa.losm;


%   fatqw=mean(fqwmo);
%   losm=;
%   Lefv=Lfmo;
%   Tefv=gtmo;
%   Tefiv=Tfimo;

%   confz=czmo;
%   confztot=czmo.*NQW.*fqwmo;
   frism=delf(modc);
%global ijos tauN0 hbar c htom
global ijos  hbar c_light htom tauN0


   om0=2*pi*1e6*c_light/lambda;
   htom=hbar*om0;
   lamod=lambda./(1+frism);

%' calopt dope ', keyboard

% set discretizzazione spaziale

  differen

% end set discretizzazione spaziale
%' % end set discretizzazione spaziale '

'cam2 ', keyboard




   if i2Ddyn==1

    for k=modc
     F=abs(reshape(Eqw(ifidif,ifisho,k),length(ifidif),lfs)).^2;
     F1=abs(reshape(Eqw(:,ifisho,k),length(xro),lfs)).^2;
     Sint=spacintv(F1,xdx',lfp,pesf,Nsett);
     cam2(k,:)=reshMv(F)'/Sint;
     Fx=abs(reshape(Eout.x(ifidif,ifisho,k),length(ifidif),lfs)).^2;
     Fy=abs(reshape(Eout.y(ifidif,ifisho,k),length(ifidif),lfs)).^2;
     Sint=spacintv(Fx+Fy,xdx1',lfp,pesf,Nsett);
     Fx=reshape(Eout.x(ifidif,ifisho,k),length(ifidif),lfs);
     Fy=reshape(Eout.y(ifidif,ifisho,k),length(ifidif),lfs);
     Eo_x(:,:,k)=Fx/sqrt(Sint);
     Eo_y(:,:,k)=Fy/sqrt(Sint);
    end

   else
     Eo_x=0;
     Eo_y=0;
     if i2D==3
'Stoet', keyboard
      for k=modc
       F=abs(reshape(Eqw(ifidif,ifisho,k),length(ifidif),lfs)).^2;;
       F1=abs(reshape(Eqw(:,ifisho,k),length(xro),lfs)).^2;;
       Sint1=spacintv(F1,xdx',lfp,pesf,Nsett);
       Sint=spacintv(F,xdx1',lfp,pesf,Nsett);
       cdp=pesf*F';
       cdp1=pesf*F1';
       cam2(k,:)=cdp/Sint1*Nsett/(2*pi);
       cam2v(k,:)=cdp1/Sint*Nsett/(2*pi);
      end

     else

       F=abs(Eqw(ifidif,:)).^2;;
       F1=abs(Eqw).^2;;
       Sint=diag(1./(2*pi*xdx*F));
       Sint1=diag(1./(2*pi*xdx1*F1));
       cam2=Sint1*F';
       cam2v=Sint*F1';

     end
   end

      ifi=0;
      omP0=losm(modc)';
      fPGA=confztot;
%      Ttot=Tefv(modc)';  % per uscita solo dall'alto
%      Nu_tot=fPGA.*tauu(modc)'/dat_tot;  % per uscita solo dall'alto
%      fPdif_old=tauN0*4*1e9./(c/rr*hbar*om0.*(1+frism).*Ttot)*1e-24;
%      fPdif=tauN0*1e9./(hbar*om0.*(1+frism))*1e-24.*Nu_tot;
      Lef=Pa.Lf;
      CP=2*Lef/(c_light/rvg);
      Ttc=CP./tauu;
%      Pp1=(tauu./Lef)';
%      fPdif=tauN0*2*Pp1(modc)*1e9./(hbar*om0.*(1+frism))*1e-24;
%      fPdif=tauN0*Pp1(modc)*1e9./(hbar*om0.*(1+frism))*1e-24;
      Teff=Pa.Teff;
      fatqw=Pa.fatqw;
      vr=c_light/rvg*1e2;
      Ppn=1./(Teff'*vr)*fatqw;
      
      vr1=c_light/rr*1e2;
      Ppn=1./(Teff'*vr1);
      Ppn=1./(Teff'*vr1)*fatqw*NQW;
      Ppn=1./(Teff'*vr1)*fatqw;

      
      
      Ppnli=1./(Teff)';
%       fPdif=tauN0*Ppn(modc)*1e9./(hbar*om0.*(1+frism))*1e-24/1000;
      fPdif=Ppn(modc)./(hbar*om0.*(1+frism)); %-- tibaldi
      global fPnli
%      'fPdif', keyboard

%      Teq=1./(1./Tefv+1./Tefiv);
%      taP=2*Lefv/(c/rr.*Teq).*fatqw;
%      taP=2*Lefv/(c/rr.*Teq).*fatqw;

      taP=tauu(modc)';
      fPnli=1e-3*Ppnli(modc)./(hbar*om0.*(1+frism)).*taP;
%            'fPdif', keyboard
      fPES=hbar*om0*(1+frism)./taP.*confztot;
%      fPES=tauN0*hbar*om0*(1+frism)./taP*1e3;
%      fPES=tauN0*hbar*om0*(1+frism)./taP*1e3.*fatqw;
      eta_eff=(taut./tauu)';
      Tper{1}=taut';
      Tper{2}=tauu';
      Tper{3}=taub';
      Tper{4}=pvol';
      Tper{5}=alca';
%' qui pvol calotp taut ', keyboard
%' qui pvol calotp taut ', keyboard
%' qui pvol calotp taut ', keyboard

if ifp>=0
' %% fine set costanti',
keyboard
end

% function [cam2,gain,lamod,modcu,xro,fian,fPdif,fPES,fPGA,tyPmod,omP0,...
% eta_eff,Tper,Ppol,Plot,Pa,Eo_x,Eo_y]=calopt

else  % idyn

 [cam2,xro,gain,lamod,modc,fPdif,fPES,fPGA,omP0]=deal(0);

end

%'fine calopt', keyboard
