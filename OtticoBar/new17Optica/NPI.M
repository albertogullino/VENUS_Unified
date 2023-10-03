   if idyn>=1
    kmo=2;
    Andu=Afz(:,kmo)+Afzf(:,kmo);
    cam_val;
    Eanuo=Etot;
   end
  Int_at=xdx*Eanuo*defi;
  Int_pow=xdx*pow*defi;

  Teff=(Int_b-Int_f)/Int_at

  keyboard
 mi=pi*4e-7;
 eps0=8.8541e-12;
 c_light=1/sqrt(mi*eps0);
 c=1/sqrt(mi*eps0);
 Z0=120*pi;
 q=1.602e-19;
 vr=c_light/rr;

 CP=2*Lef/(vr);
 Ttc=CP./tauu;
 Pp1=(tauu./Lef)'

  Ppn=1/(Teff*vr)*fatqw/NQW

% fPdif=tauN0*2*Pp1(modc)*1e9./(hbar*om0.*(1+frism))*1e-24;
 fPdif=tauN0*Pp1(modc)*1e9./(hbar*om0.*(1+frism))*1e-24;





  if N==0
   Ec_at(fiGa)=0;
  else
   Nga=[N; zeros(length(x)-length(ro_in),1)];
   Nmat=repmat(Nga,1,length(fian));
%   NNor=xdx*Nga*2*pi;
%   Are=pi*avero^2;
   NNor=1;
   Are=1;
   Ec_at=Ec_at.*Nmat/NNor*Are;
  end

  Int_at=xdx*Ec_at*defi;
  Gam_at=Int_at/Int_cb;



  cmi=c*1e6;
  omega=cmi*k0*(1+fris);
  Perdvol=-omega*sum(Perd)/confz;
  if iplan==0
   Glostotal=2*gg0*Gam_at+Perdvol*(1-iPERD_vol);
  else
   Glostotal=2*gg0+Perdvol*(1-iPERD_vol);
  end
  Lsov(imod)=Glostotal*confz;

  Ec_su=Etcsum;
  Int_su=xdx*Ec_su*defi;

  fa4=Int_su/Int_cb*fatqw;
  Gam_v(imod)=Gam_at;

  if i2D==3
   Ec_at=Etf(:,:,2);
  else
   Ec_at=Etf(:,2);
  end
  Int_cf=xdx*Ec_at*defi;

%   disp(' Completamente planare '), keyboard
  if ifp==-10
%   keyboard
   pausak
  end

  if ifp>-4 | ifp==-10
   disp('  ')
   disp('  ')
   disp(' Completamente planare ')
   disp('  ')
  end
  rg=rr;
  rg=nmean;
  CP=c/rg/(2*Lef*confz);
  CP0=c/rg/(2*Lef);
  Gsup=fTras*CP
  Ginf=fTrasinf*CP
  Gper=(fPinf+fPsup)*CP/2
  Gplanar=Gsup+Ginf+Gper

  if ifp>-4 | ifp==-10
   disp('  ')
   disp('  ')
   disp('  ')
   disp(' Planare corretto  da integrali')
   disp('  ')
  end

  if i2D==3
   Ef=Etf(:,:,1);
   Eb=Etb(:,:,1);
  else
   Ef=Etf(:,1);
   Eb=Etb(:,1);
  end
  Efu=Ef;
  Ebu=Eb;

  Int_f=xdx*Ef*defi;
  Int_b=xdx*Eb*defi;

  Goutsup=(Int_b-Int_f)*CP/Int_cb
%  Goutsup=(Int_b+Int_f)*CP/Int_cb




  if i2D==3
   Ef=Etf(:,:,3);
   Eb=Etb(:,:,3);
  else
   Ef=Etf(:,3);
   Eb=Etb(:,3);
  end
  Int_f=xdx*Ef*defi;
  Int_b=xdx*Eb*defi;
  Goutinf=(Int_f-Int_b)*CP/Int_cf
%  Goutinf=(Int_f+Int_b)*CP/Int_cf

  if istmet==0
   istmet1=[];
  else
   istmet1=istmet;
  end
  if length(istmet1)==1
   if i2D==3
    Ef=Etf(:,:,istmet);
    Eb=Etb(:,:,istmet);
   else
    Ef=Etf(:,istmet);
    Eb=Etb(:,istmet);
   end
   Int_f=xdx*Ef*defi;
   Int_b=xdx*Eb*defi;
   Gtransmet=-(Int_f-Int_b)*CP/Int_cb;
   Perdring=Gtransmet-Goutsup
  else
   Perdring=0;
  end

  Perdrad=Goutsup+Goutinf;

  Perdvol_u=-omega*Perd(1)/confz;
  Perdvol_d=-omega*Perd(2)/confz;
  Gapprox=(Perdvol*iPERD_vol+Perdrad+Perdring);
   disp(' confronto approx  numerico ')
   [Gapprox glosout*Gam_at]

  if ifp==-10
%   keyboard
   pausak
  end

  if ifp>-4 | ifp==-10
   disp('  ')
   disp('  ')
   disp('  ')
   disp('[total= pred_vol + out_sup + out_inf]');
   disp('  ')
   disp(' Planar ');
   [Gplanar Gper Gsup Ginf ]

   disp('  ')
   disp(' integral method ');
   [Gapprox  Perdvol Goutsup Goutinf]

   disp(' Confronto completamente planare: verifica '),
  end


