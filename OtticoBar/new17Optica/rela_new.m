'entro realanew', keyboard
imem=1;   %=1 velocizza le matrici diagonali
%imem=0;   %=1 velocizza le matrici diagonali
icmem=0;

% set global IdeOo pMu0u pMu0u1 lKA nk1max nures

if iLP==0
 lKA=nk1max*numodi*2;
else
 lKA=nk1max*numodi;
end

IdeOo=zeros(2*lKA,2*lKA);

%  ld=[0 1 2 numodi+[0 1 2] ]*nk1max;
  nures=numodi*4;
  numoa=0:numodi-1;
  ld=[numoa numodi+numoa]*nk1max;
  ldu=[ld ld+lKA];
  pMu0u=[];
  for kdu=1:nures
   pMu0u=[pMu0u ldu+(kdu-1)*2*lKA*nk1max];
  end

  pMu0u1=[0 2*lKA^2; lKA lKA+2*lKA^2];


istatt=0;
flp=1-iLP1;
mr=1;
Ide=diag(ones(1,2*length(be)));


% devo calcolare le matrici di trasmissione sopra e sotto
%
% sopra
ictr=0;

 Li=Litn;
 ni=nitn;
 ai=aitn;
 bi=aral.y.t';
 pai=aral.p.t';
 shai=shav.t';
 icrit=zeros(size(Li));
 ficrip=find(iauto(:,1)==-10)-1;
 icrit(1:ficrip)=1;
 
'sopra icrit ', keyboard
%'sopra ', keyboard
%'sopra ', keyboard
 
 istfie=istfi.t;
 tyari=aral.ar.t';

 fmlsi=fmlst;
 T=Ide;
 Tsave=Ide;
 Tmir=Ide;
 Tss=Ide;
 ilaymem=find(diff(ai)~=0);


 if iany==0
  Kan=0;
  ani=zeros(size(Li));
 elseif iany==1
  Kan=Kosan;
  ani=deltanyn*fattany;
 elseif iany==2
  Kan=Kosan*diag(pes);
  ani=deltanyn*fattany;
 end

 if ianys~=0
  if ianys-iany==0
   ani=ani+deltanyns;
   ani1=zeros(size(ani));
   Kan1=0;
  else
   if ianys==1
    Kan1=Kosan1;
    ani1=deltanyns;
   elseif ianys==2
    Kan1=Kosan1*diag(pes);
    ani1=deltanyns;
   end
  end
 else
  Kan1=0;
  ani1=zeros(size(ani));
 end

disp(' rala sopra'), keyboard

if igau==4
% disp('rel')
% keyboard
 fmol=fmlsi(:,2);
 fi=find(fmol==0 & real(ni(1,:))'~=1);
 fmol(fi)=ones(size(fi));
 zrif=Li.*fmol;
end


 istr=1;
 isomv=0;

 itutmir0=1;
 if igau==4
  itutmir=itutmir0;
 else
  itutmir=0;
 end
 icaco=1;
 ifplatot=1;
% disp(' prima del while '), keyboard


 while istr<=length(Li)
   if ifp>=1
    disp(' while sopra: istr'), [istr]
   end
  icop=fmlsi(istr,1);
  ncop=fmlsi(istr,2);

  istr0=istr;
  if itutmir==1
   if ncop~=0
    Ncop=ncop;
   else
    Ncop=1;
   end
  else
   Ncop=1;
  end
  iziv=1;
  ncopF=1;
  istrctv=0;
  ncopco=0;

  ifplaco=1;

  for istrct=1:Ncop

   istr=istr0;

   if ncopco==0
    Tc=Ide;
   end
   ncopco=ncopco+1;
   for istrc=1:icop
    if ifp>=1
     disp(' for: istr, istrc'), [istr istrc], pausak
    end
    if igau==4
     isomv=isomv+Li(istr);
     [du,izi]=min(abs(isomv-zedis(2:length(zedis))));
      KTe=reshape(KTemp(:,:,izi),si);
      KTez=reshape(KTempz(:,:,izi),si);

    end
%   disp(' [istrct istrc izi icaco istr]'),
%   [istrct istrc izi icaco istr]
    if icaco==1

     eval(emme)

     if istrc==1
      Tcm=Tc;
     end
    end
    istr=istr+1;
   end  %istrc
   if igau==4 & itutmir==1
    if Ncop>1
     if izi~=iziv | istrct==Ncop
      icaco=1;
     else
      icaco=0;
     end
     if icaco==1
      ncopF=ncopco;
      ncopco=0;
     end
    else
     icaco=1;
     ncopco=0;
    end
    iziv=izi;
    istrctv=istrct;
   else
    icaco=1;
    ncopF=ncop;
   end

   if icaco==1

    if ncop<=1
     if imem==0
      T=Tc*T;
     else
      T=prodmat(Tc,T,ifplatot);
     end
    else

     if ncopF-fix(ncopF)==0
      if imem==0
       Tco=Tc^ncopF;
      else
       Tco=powermat(Tc,ncopF,ifplaco);
      end
%       ncopF
     else
      if imem==0
       Tco=Tcm*Tc^fix(ncopF);
      else
       Pow=powermat(Tc,fix(ncopF),ifplaco);
       Tco=prodmat(Tcm,Pow,ifplaco);
      end
     end

     if istr==pucavi & pucavi>1
      Tss=Tco;
     end

     if ncop>1
      if imem==0
       Tmir=Tco*Tmir;
      else
       Tmir=prodmat(Tco,Tmir,ifplatot);
      end
     end
     if imem==0
      T=Tco*T;
     else
      T=prodmat(Tco,T,ifplatot);
     end
    end
   end  %icaco
%   disp(' istrct '), istrct
%   pausak
  end  %istrct
  if length(find(istfie==istr))==1
   icmem=icmem+1;
   Tmef(:,:,icmem)=T;
   disp(' memorizzo  Tmef')
   keyboard
  end
 end   %istr

 Tmirup=Tmir;
 if length(ilaymem)==0
  Tt=T;
 else
  Tt=T*Tsave;
 end

%%%%   sotto

 Li=Lib;
 ni=nib;
 ai=aib;
 bi=aral.y.b';
 pai=aral.p.b';
 shai=shav.b';
 tyari=aral.ar.b';
 istfie=istfi.b;
 icmem=0;
 icrit=zeros(size(Li));
 fmlsi=fmlsb;
 T=Ide;
 ilaymem=find(diff(ai)~=0);


 if iany==0
  Kan=0;
  ani=zeros(size(Li));
 elseif iany==1
  Kan=Kosan;
  ani=deltanyb*fattany;
 elseif iany==2
  Kan=Kosan*diag(pes);
  ani=deltanyb*fattany;
 end

 if ianys~=0
  if ianys-iany==0
   ani=ani+deltanybs;
   ani1=zeros(size(ani));
   Kan1=0;
  else
   if ianys==1
    Kan1=Kosan1;
    ani1=deltanybs;
   elseif ianys==2
    Kan1=Kosan1*diag(pes);
    ani1=deltanbns;
   end
  end
 else
  Kan1=0;
  ani1=zeros(size(ani));
 end


if igau==4
% disp('rel')
% keyboard
 fmol=fmlsi(:,2);
 fi=find(fmol==0 & real(ni(1,:))'~=1);
 fmol(fi)=ones(size(fi));
 zrif=Li.*fmol;
 isomv=isomv+d*1e6;
end

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%disp(' relautan')
%keyboard

 Tmir=Ide;
 T=Ide;
 Tsave=Ide;
 istr=1;
 ifplatot=1;


 while istr<=length(Li)

   if ifp>=1
    disp(' while sopra: istr'), [istr]
   end
  icop=fmlsi(istr,1);
  ncop=fmlsi(istr,2);
  istr0=istr;
  if itutmir==1
   if ncop~=0
    Ncop=ncop;
   else
    Ncop=1;
   end
  else
   Ncop=1;
  end
  iziv=1;
  ncopF=1;
  istrctv=0;
  ncopco=0;

  ifplaco=1;

  for istrct=1:Ncop

   istr=istr0;
   if ncopco==0
    Tc=Ide;
   end
   ncopco=ncopco+1;
   for istrc=1:icop
    if ifp>=1
     disp(' for: istr, istrc'), [istr istrc], pausak
    end
    if igau==4
     isomv=isomv+Li(istr);
     [du,izi]=min(abs(isomv-zedis(2:length(zedis))));
      KTe=reshape(KTemp(:,:,izi),si);
      KTez=reshape(KTempz(:,:,izi),si);
    end
    if icaco==1

     eval(emme)

     if istrc==1
      Tcm=Tc;
     end
    end
%   disp(' [istrct istrc izi icaco istr]'),
%   [istrct istrc  icaco istr]
%   pausak
    istr=istr+1;
   end  %istrc
   if igau==4 & itutmir==1
    if Ncop>1
     if izi~=iziv | istrct==Ncop
      icaco=1;
     else
      icaco=0;
     end
     if icaco==1
      ncopF=ncopco;
      ncopco=0;
     end
    else
     icaco=1;
     ncopco=0;
    end
    iziv=izi;
    istrctv=istrct;
   else
    icaco=1;
    ncopF=ncop;
   end

   if icaco==1
%    if ncop>1
%     disp(' icaco ')
%     keyboard
%    end

    if ncop<=1
     if imem==0
      T=Tc*T;
     else
      T=prodmat(Tc,T,ifplatot);
     end
    else

     if ncopF-fix(ncopF)==0
      if imem==0
       Tco=Tc^ncopF;
      else
       Tco=powermat(Tc,ncopF,ifplaco);
      end
%       ncopF
     else
      if imem==0
       Tco=Tcm*Tc^fix(ncopF);
      else
       Pow=powermat(Tc,fix(ncopF),ifplaco);
       Tco=prodmat(Tcm,Pow,ifplaco);
      end
     end
     if ncop>1

      if imem==0
       Tmir=Tco*Tmir;
      else
       Tmir=prodmat(Tco,Tmir,ifplatot);
      end
     end
     if imem==0
      T=Tco*T;
     else
      T=prodmat(Tco,T,ifplatot);
     end
    end
   end  %icaco
%   disp(' istrct '), istrct
%   pausak
  end  %istrct


  if length(find(istfie==istr))==1
   icmem=icmem+1;
   Tmefb(:,:,icmem)=T;
%   disp(' memorizzo  Tmefb')
%   keyboard
  end

 end   %istr

 if length(ilaymem)==0
  Tb=T;
 else
  Tb=T*Tsave;
 end

 Tmirdw=Tmir;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Ga1=diag(Gas);
  Ga2=diag(Gad);

%% calcolate le matrici di trasmissione sopra e sotto, calcolo l'eq. agli autov.

 clear alphavv Gvav Anz Anzf CCd

   pu1=1:length(Gas);
   pu2=length(Gas)+1:2*length(Gas);

  inv_mix

% disp(' dopo inv_mix ')
% keyboard

   if isopro==1
    return
   end

     [alphavv,iav]=sort(al(igv));
     nmodv=[1:length(iav)];
     Nd1vu=igv(iav);
     Gvav=Gv(Nd1vu);

     if ifp>=-1
     disp(' fr , Gvav, alphav')
     [ Gvav alphavv]
     end

   if isopro>1
    return
   end

     iNc=0;

    for iN=Nd1vu'
      iNc=iNc+1;
      Ar=V(:,iN);
      Lam=Lamv(iN);
      Apr=[Gas.*Ar; Ar];
      CCd=Apr;
      ics=1;
      Ab=CCd(pu2);
      Af=CCd(pu1);
      Anz(:,iNc,ics)=Ab.*Trs;
      Anzf(:,iNc,ics)=Af.*Trs;

      CCd=Tt*Apr;
      ics=ics+1;
      Ab=CCd(pu2);
      Af=CCd(pu1);
      Anz(:,iNc,ics)=Ab;
      Anzf(:,iNc,ics)=Af;

      Ttotal=(Mn+Lam*Mi);
      CCd=Ttotal*Apr;
      ics=ics+1;
      Ab=CCd(pu2);
      Af=CCd(pu1);
      Anz(:,iNc,ics)=Ab.*Trd;
      Anzf(:,iNc,ics)=Af.*Trd;


      if exist('Tmef')
       sdu=size(Tmef);
       sTm=sdu(1:2);
      end


%  sopra
      if icsfit>0
       for kf=1:icsfit
        ics=ics+1;
        CCd=reshape(Tmef(:,:,kf),sTm)*Apr;
        Ab=CCd(pu2);
        Af=CCd(pu1);
        Anz(:,iNc,ics)=Ab;
        Anzf(:,iNc,ics)=Af;
       end
      end

      if exist('Tmefb')
       sdu=size(Tmefb);
       sTm=sdu(1:2);
      end

%  sotto
      if icsfib>0
       TQW=(Mon+Lam*Moi);
       AprQW=TQW*Tt*Apr;
       for kf=1:icsfib
        ics=ics+1;
        CCd=reshape(Tmefb(:,:,kf),sTm)*AprQW;
        Ab=CCd(pu2);
        Af=CCd(pu1);
        Anz(:,iNc,ics)=Ab;
        Anzf(:,iNc,ics)=Af;
       end
      end

    end   %iNc
    disp(' rela_new ')

%    keyboard
     if ifp>=1, pausak, end

if ifp>=0, pausak, end
