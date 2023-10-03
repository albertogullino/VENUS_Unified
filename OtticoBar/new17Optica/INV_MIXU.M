%'invmixu', keyboard
if iplan==0

  Kosi=KA;
  if iztm==1
   Koszi=KAz;
   KOt=Kosi*diag(pes);
   KOz=Koszi*diag(pes);
%   Ksu=KOt+KOz;
%   Kdi=KOt-KOz;
%   P=[Ksu Kdi; -Kdi -Ksu];
   P=[KOt+KOz  KOt-KOz; -(KOt+KOz)  -(KOt-KOz)];
  else
   KOt=Kosi*diag(pes);
   P=[KOt KOt; -KOt -KOt];
  end
  Moi=P(Pust,Pust);

 DGp=-g0*1e9;  %in ampiezza
 if ipilat==0
  dea=(real(ra)^2-rr^2)/rr^2+2j*imag(ra)/rr;
  ck=-j*kcav*d*1e6;
  pI=ck*(dea*Idelta+deltany*Ideltaa/4);  % fattore 4 per la normalizzazione

    global del_n_ag ianti_gui

    if ~exist('KA_ag')
     KA_ag=KA;
    end
    iplat=-10;
     if ianti_gui==1
      iplat=0;
      global del_n_ag
      del_eps_ag=(2*ra*del_n_ag+del_n_ag.^2)/rr^2;
      del_eps_ag=2*del_n_ag/rr;
      pIn=ck*del_eps_ag*KA_ag*diag(pes);
%      ' del_eps ', keyboard
      if iLP==0
       pInz=ck*del_eps_ag*KAz_ag*diag(pes);
      else
       pInz=0;
      end
      pI=pI+pIn;
     else
      pInz=0;
     end

  pIz=ck*dea*Ideltaz+pInz;
  Belo=-j*be*kcav*d*1e6*fatd;


%  Ksu=pI+pIz*flp;
%  Kdi=pI-pIz*flp;
%  P=[Ksu Kdi; -Kdi -Ksu];
  P=[pI+pIz*flp  pI-pIz*flp; -(pI+pIz*flp)  -(pI-pIz*flp)];
%  Pe=P+diag(([Belo; -Belo]));
%  iama=[];
%  Mon=expm_mio(Pe,iplat,iama,Pust);
  P=P+diag(([Belo; -Belo]));
  iama=[];
  Mon=expm_mio(P,iplat,iama,Pust);

 else

  fiat0=find(aitot~=0);
  aitd0=aitot(fiat0);
  fiat=find(aitd0==apilat);
  istatt=0;
  ipri=fiat(1);
  istr=1;
  Li=d*1e6;
  ni=niat;
  ai=aiat;
  bi=aral.y.a';
  pai=aral.p.a';
  shai=shav.a';
  tyari=aral.ar.a';

 if iany==0
  Kan=0;
  ani=zeros(size(Li));
 elseif iany==1
  Kan=Kosan;
  ani=aniat*fattany;
 elseif iany==2
  Kan=Kosan*diag(pes);
  ani=aniat*fattany;
 end


 if ianys~=0
  if ianys-iany==0
   ani=ani+aniats;
   ani1=0;
   Kan1=0;
  else
   if ianys==1
    Kan1=Kosan1;
    ani1=aniats;
   elseif ianys==2
    Kan1=Kosan1*diag(pes);
    ani1=aniats;
   end
  end
 else
  Kan1=0;
  ani1=0;
 end

  istatt

if inuo==1
 fisem=find(ai-bi~=0);
 asem=zeros(size(ai));
 asem(fisem)=ai(fisem);
 ilaymem=find(diff([asem 0])<0 & shai>1);
 if pol==pvet(1)
  layskip=zeros(size(ai));
 else
  filay=find(asem~=0 & shai>1);
  layskip=ones(size(ai));
  layskip(filay)=0;
 end

else

 layskip=zeros(size(ai));
end

disp(' attivo '), pausak
  if layskip==0
   eval(emme)

   istatt=0;


   Mon=Oo;
  end
disp(' dopo attivo '), pausak

%  eval(emme)
%  disp(' verifica att. '), keyboard
 end

else  %iplan =1
%  disp(' passo 3')
% DGp=0;
% ra0=real(ra);
 DGp=-g0*1e9;  %in ampiezza
 ra0=(ra);
% dea=(ra0^2-rr^2)/rr^2;
 dea=2j*imag(ra0)/rr;
 ck=-j*kcav*d*1e6;
 pI=ck*dea*Idelta;
 pIz=ck*dea*Ideltaz;
 Belo=-j*be*kcav*d*1e6;
% Ksu=pI+pIz*flp;
% Kdi=pI-pIz*flp;
% P=[Ksu Kdi; -Kdi -Ksu];

  P=[pI+pIz*flp  pI-pIz*flp; -(pI+pIz*flp)  -(pI-pIz*flp)];

% Pe=[Ksu Kdi; -Kdi -Ksu];
%
% ded=(real(ra0)^2-rr^2)/rr^2;
% pI=ck*ded*Idelta;
% pIz=ck*ded*Ideltaz;
% Ksu=pI+pIz*flp;
% Kdi=pI-pIz*flp;
% Pd=[Ksu Kdi; -Kdi -Ksu];
%
% Mon=diag(exp([Belo; -Belo]))+Pe+Pd;
%
% Mon=Mon(Pust,Pust);
%

 ded=(real(ra0)^2-rr^2)/rr^2;
 pI=ck*ded*Idelta;
 pIz=ck*ded*Ideltaz;
% Ksu=pI+pIz*flp;
% Kdi=pI-pIz*flp;
% P=P+[Ksu Kdi; -Kdi -Ksu];
 P=P+[pI+pIz*flp  pI-pIz*flp; -(pI+pIz*flp)  -(pI-pIz*flp)];

 P=diag(exp([Belo; -Belo]))+P;

 Mon=P(Pust,Pust);

 if iztm==1
  KOt=Idelta;
  KOz=Ideltaz;
%  Ksu=KOt+KOz;
%  Kdi=KOt-KOz;
%  P=[Ksu Kdi; -Kdi -Ksu];
  P=[KOt+KOz  KOt-KOz; -(KOt+KOz)  -(KOt-KOz)];
 else
  KOt=Idelta;
  P=[KOt KOt; -KOt -KOt];
 end

 Moi=P(Pust,Pust);
end

'MON', keyboard
