    TL_ef
%  CP=c/rr/(2*Lef*confz);
  CP=c/nmean/(2*Lef*confz);

  if i2D==3
   Ec_at=Etb(:,:,2);
  else
   Ec_at=Etb(:,2);
  end
  Ec_atb=Ec_at;

  Int_c=xdx*Ec_at*defi;
  Ec_at(fiGa)=0;
  Ec_aa=Ec_at;
  Int_at=xdx*Ec_at*defi;
  Gam_at=Int_at/Int_c;

  if i2D==3
   Ec_at=Etf(:,:,2);
  else
   Ec_at=Etf(:,2);
  end
  Ec_atf=Ec_at;

  Int_cf=xdx*Ec_at*defi;
  G_sup=(1-Int_cf/Int_c)*CP

  if i2D==3
   Ec_at=Etb(:,:,3);
  else
   Ec_at=Etb(:,3);
  end

  Int_i=xdx*Ec_at*defi;

  if i2D==3
   Ec_at=Etf(:,:,3);
  else
   Ec_at=Etf(:,3);
  end

  Int_if=xdx*Ec_at*defi;
  G_inf=(Int_if/Int_i-1)*CP



  Ec_su=Etcsum;
  Int_su=xdx*Ec_su*defi;

  fa4=Int_su/Int_cb*fatqw;
  Gam_v(imod)=Gam_at;


  if ifp>-4 | ifp==-10
   disp('  ')
   disp('  ')
   disp(' Completamente planare ')
   disp('  ')
  end

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

  Int_f=xdx*Ef*defi;
  Int_b=xdx*Eb*defi;

%  Goutsup=-(Int_f-Int_b)*CP/Int_cf
  Goutsup=-(Int_f-Int_b)*CP/Int_c


  if i2D==3
   Ef=Etf(:,:,3);
   Eb=Etb(:,:,3);
  else
   Ef=Etf(:,3);
   Eb=Etb(:,3);
  end
  Int_f=xdx*Ef*defi;
  Int_b=xdx*Eb*defi;
%  Goutinf=(Int_f-Int_b)*CP/Int_c
  Goutinf=(Int_f-Int_b)*CP/Int_cf

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
   Gtransmet=-(Int_f-Int_b)*CP/Int_c;
   Perdring=Gtransmet-Goutsup
  else
   Perdring=0;
  end

  Perdrad=Goutsup+Goutinf;
  cmi=c*1e6;
  omega=cmi*k0*(1+fris);
  Perdvol=-omega*sum(Perd)/confz
  Gapprox=(Perdvol+Perdrad+Perdring);
   disp(' confronto approx  numerico ')
   [Gapprox glosout*Gam_at]

  if ifp==-10
   keyboard
  end

  if ifp>-4 | ifp==-10
   disp('  ')
   disp('  ')
   disp('  ')
   disp('[total pred_vol out_sup  out_inf]');
   disp('  ')
   disp(' Planar ');
   [Gplanar Gper Gsup Ginf ]

   disp('  ')
   disp(' integral method ');
   [Gapprox  Perdvol Goutsup Goutinf]

   disp(' Confronto numerico, planare, integrale '),
   [glosout*Gam_at Gplanar Gapprox]
  end

