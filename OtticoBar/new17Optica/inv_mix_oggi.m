%' prima di inv_mix ' ,keyboard
%if exist('gpla')
if ~exist('gpla')
 gpla=gth1d*NQW_ef;
end
 %if gpla<2000
 gplam=gpla;
 %else
 %gplam=0;
 %end
 g0=gplam/2;          % include anche perdite QW
% g0=0
if ifr==1
ra=niat(1);
ra_centrale=ra;
global del_n_ag ianti_gui
 if ianti_gui==1
  ra_centrale=niat(1)+del_n_ag*NQW_ef;
 end 
% DGp=-g0*1e9;  %in ampiezza
 pqw=-imag(niat(1))*k0*1e4; 
 
 DGp=g0;  %in ampiezza
 DGp_step=g0;  %in ampiezza
 r_attivo=niat(1); 
% ra_centrale=real(ra_centrale)+j*g0/k0*1e-4; 
 ra_centrale=real(ra_centrale)+j*(DGp_step)/k0*1e-4; 
 if iplan==1
  ra=ra_centrale;
 end
 niat(1)=ra;
end
%g0=-gpla/2;

%' contre g0', keyboard


if iplan~=1
%
%' prima di inv_mix ' ,keyboard
  Kosi=KA;
  if iztm==1
   Koszi=KAz;
   KOt=Kosi*diag(pes);
   KOz=(rr/ra_centrale)^4*Koszi*diag(pes);
%           KOt(1:20,:)=0; 
%           KOt(:,1:20)=0; 
%           KOz(1:20,:)=0; 
%           KOz(:,1:20)=0; 
%'passatp KOt', keyboard
   Ksu=KOt+KOz;
   Kdi=KOt-KOz;
   P=[Ksu Kdi; -Kdi -Ksu];
  else
   KOt=Kosi*diag(pes);
   P=[KOt KOt; -KOt -KOt];
  end
  Moi=P(Pust,Pust);
%  'Moi', keyboard
if length(find(niat>0))>1
% ipilat=1;
 apilat=aiat;
 ' da sistemare il case attivo con altro .... per ora ae posto infinito '
%  'Moi'
%keyboard 
end 

 if ipilat==0 

  dea1=(ra_centrale^2-ra^2)/rr^2;
  deaz=1-(rr/ra)^2; 
%  deaz1=2*ra_centrale/ra*(rr/ra)^2; 
  deaz1=0;


  dea=(ra^2-rr^2)/rr^2;
  drc=ra_centrale-ra;
  dea_c=(2*ra+drc)*drc/rr^2;
  deaz_c=(rr/ra_centrale)^4*dea_c;    
  ck=-j*kcav*d*1e6;
  pI=ck*(dea*Idelta+deltany*Ideltaa/4);  % fattore 4 per la normalizzazione
  
  pIn=ck*dea_c* KA*diag(pes);

  pIn1=pI+ck*dea1* KA*diag(pes);

  
    if ~exist('KA_ag')
     KA_ag=KA;
     KAz_ag=KAz;
    end
    iplat=-10;
     if ianti_gui==1
      iplat=0;
     end 

%      ' del_eps ', keyboard
%      keyboard
      if iLP==0
       pInz=ck*deaz_c*KAz*diag(pes);
       pInz1=ck*deaz1*KAz*diag(pes);
      else
       pInz=0;
       pInz1=0;
      end
      
      
      pI=pI+pIn;
      pIz=ck*deaz*Ideltaz+pInz;
      
      pInz1=ck*deaz*Ideltaz+pInz1;


  Belo=-j*be*kcav*d*1e6*fatd;


  Ksu=pI+pIz*flp;
  Kdi=pI-pIz*flp;
  
  Ksu1=pIn1+pInz1*flp;
  Kdi1=pIn1-pInz1*flp;
    
  
  Pk=[Ksu Kdi; -Kdi -Ksu];

  Pk1=[Ksu1 Kdi1; -Kdi1 -Ksu1];
%  Pe=P+diag(([Belo; -Belo]));
%  iama=[];
%  Mon=expm_mio(Pe,iplat,iama,Pust);
  P=Pk+diag(([Belo; -Belo]));
  P1=Pk1+diag(([Belo; -Belo]));
%  'ferma Pk', keyboard
%   [du,D]=eig(P);
% Ood = du * diag(exp(diag(D))) / du;
%  [edu,eiga]=eig(P);
%  evet=diag(eiga);
%  ei=imag(evet);
%  fi=find(ei<0);
%  eim=ei(fi);
%  figure, plot(ei(fi))
%  fi=find(ei>0);
%  eip=ei(fi);
%  figure, plot(ei(fi))

  iama=[];
  Monv=expm_mio(P,iplat,iama,Pust);
  Mon1=expm_mio(P1,iplat,iama,Pust);
  Mon=Monv;
  Mon=Mon1;
  %'Mon prima', keyboard  
  Moi0=Moi;
  Moi=Mon*Moi0;       % esce dal calcolo esatto della linearizzazione
%  Moi=Moi0*Mon;       % esce dal calcolo esatto della linearizzazione
%  'Moi, Mon', keyboard  
 else

'questo ramo in inv_mix', keyboard
  fiat0=find(aitot~=0);
  aitd0=aitot(fiat0);
  fiat=find(aitd0==apilat(1));
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
 ilaymem=find(diff([asem; 0])<0 & shai>1);
 if pol==pvet(1)
  layskip=zeros(size(ai));
 else
  filay=find(asem~=0 & shai>1);
  layskip=ones(size(ai));
  layskip(filay)=0;
 end

else

 layskip=0;
 ilayfast=0;
end

%disp(' attivo '), pausak
'attivo', layskip
%disp(' attivo '), keyboard
%iatttivo=1;
if ifr==3
%  'passo emme', keyboard
end
  if layskip==0
%  'passo emme', keyboard
   eval(emme)

   istatt=0;


   Mon=Oo;
  end
%disp(' dopo attivo '), pausak

%  eval(emme)
   if ifp==-11
    disp(' verifica att. '), keyboard
   end 
 end
 if length(niat)>1
  ipilat=0;
  apilat=aiat;
 %  'Moi'
 %keyboard 
  end

else  %iplan =1

  if idyn==1
    dea=(real(ra)^2-rr^2)/rr^2+2j*imag(ra)/rr;
    ck=-j*kcav*d*1e6;
    pI=ck*(dea*Idelta+deltany*Ideltaa/4);  % fattore 4 per la normalizzazione
    global del_n_ag ianti_gui

    if ~exist('KA_ag')
     KA_ag=KA;
     KAz_ag=KAz;
    end
    iplat=-10;
     if ianti_gui==1
%      iplat=0;
%      global del_n_ag
%      del_eps_ag=(2*ra*del_n_ag+del_n_ag.^2)/rr^2;
      del_eps_ag=-2*del_n_ag/rr;
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


    Ksu=pI+pIz*flp;
    Kdi=pI-pIz*flp;
    P=[Ksu Kdi; -Kdi -Ksu];
    P=P+diag(([Belo; -Belo]));
    iama=[];
    Mon=expm_mio(P,iplat,iama,Pust);

  else %idyn

   DGp=-g0*1e9;  %in ampiezza
   DGp=-g0;  %in ampiezza
   
   ra0=(ra);
   ck=-j*kcav*d*1e6;
   Belo=-j*be*kcav*d*1e6;   
   ded=(ra0^2-rr^2)/rr^2;
   dedz=1-(ra/rr)^2;
   pI=ck*ded*Idelta;
   pIz=ck*dedz*Ideltaz;
   Ksu=pI+pIz*flp;
   Kdi=pI-pIz*flp;
   Psal=[Ksu Kdi; -Kdi -Ksu];

   P=diag([Belo; -Belo])+Psal;
   
   Mon=expm(P(Pust,Pust));
%' qui mon', keyboard
  end %idyn

 if iztm==1
  KOt=Idelta;
  KOz=(rr/ra).^4*Ideltaz;
  Ksu=KOt+KOz;
  Kdi=KOt-KOz;
  P=[Ksu Kdi; -Kdi -Ksu];
 else
  KOt=Idelta;
  P=[KOt KOt; -KOt -KOt];
 end

 Moi=Mon*P(Pust,Pust);
 %'Altro Moi', keyboard
end


if ifiez==1
%' prima di invpro ' ,keyboard
 Mois=Moi;
 Mons=Mon;
 return
end
if ifp==-10
 % disp(' prima di invpro '), keyboard
end
  global ired_ret
%  'cont Pured', keyboard
 if length(ired_ret)==1
  if ired_ret==1
  
   Tts=Tt;
   Tbs=Tb;
   Ga1s=Ga1;
   Ga2s=Ga2;
   Gass=Gas;
   Gads=Gad;
   Mons=Mon;
   Mois=Moi;
   Moi=Moi(Pured,Pured);
   Mon=Mon(Pured,Pured);
   Tt=Tt(Pured,Pured);
   Tb=Tb(Pured,Pured);
   if pol==pvet(1)
    Ga1=Ga1(Pured1,Pured1);
    Ga2=Ga2(Pured1,Pured1);
    Gas=Gas(Pured1);
    Gad=Gad(Pured1);
    pu1=1:length(Gas);
    pu2=length(Gas)+1:2*length(Gas);
   end 
  end
 end
  Mn=Tb*Mon*Tt;
  Mi=Tb*Moi*Tt;
  
%Tn=Td*Tan*Ts;
%Ti=Td*Tai*Ts;  

%  verifica conversione interfacce
%   I=eye(size(Ga1));
%   V1=[I Ga1; Ga1 I]; 
%   V2=[I -Ga2; -Ga2 I]; 
%   MnR=V2*Mn*V1;

%  'Mi, Mn', keyboard
 if isfield(Ps,'istop_fie')==1
  if Ps.istop_fie==1
%   ' Mn Mi ', keyboard
  end
 end 
  if ifr>=2
%   ' Mn Mi ', keyboard
  end

if ifp==-11
  disp(' prima di invpro '), keyboard
end
%  clear P Pe

%' prima invpro', keyboard

  invpro

if ifp>=-1,  pausak, end
