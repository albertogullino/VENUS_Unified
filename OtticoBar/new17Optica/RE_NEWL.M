' *****   re_newLLLLLL '
dpes=diag(pes);
iant=0;
if iant==0
 dpes1=1;
 dpes2=dpes;
else
 dpes2=1;
 dpes1=dpes;
end

iTsav=1;
%iTsav=0;
if length(pvet)==2
 pvetch=pvet;
else
 pvetch=[pvet -100];
end

ifast=1;
%ifast=0;
isem=1;   %=1 metodo veloce per matrici

if ispeed==0
 ifast=0;
 isem=0;   %=1 metodo veloce per matrici
 iTsav=0;   %=1 metodo veloce per matrici
end


lfishav=length(find(shavet(:,1)~=0));
sMAT=lfishav*length(Pust)^2*16;
Mb=2^20;
%Mb=2^10;
if length('if_only_out')>0
 Mb=0;
% iTsav=1;
end

if ifr==1 & pol==pvet(1)
 if ifp==-12
  ifp=-10;
 end
end

%' salva controllo ', keyboard
siMAX=20;
%siMAX=0.1;

if sMAT<=Mb*siMAX
 iTsav=0;
end
%'%%%%%%%%%% setto per prova iTsav '


if iTsav==1 & ifr==1 & pol==pvet(1) & sMAT>Mb*siMAX
  ck=clock;
  Dire=num2str(floor(prod(ck( find(ck~=0) ))*1e-4));
  Droo=cd;
  eval(['!md ',Droo,'\',Dire]);
  Dsav=[Droo,'\',Dire];
%  disp(' chain: layskip '), keyboard
%else
% iTsav=0;
end




%if iTsav==1 & ifr==1 & pol==pvet(1)
%  ck=clock;
%  Dire=num2str(floor(prod(ck( find(ck~=0) ))*1e-4));
%  Droo=cd;
%  eval(['!md ',Droo,'\',Dire]);
%  Dsav=[Droo,'\',Dire];
%end

ifastsa=-1;
if ifast==1
 if ifr==1
  ifastsa=1;
 else
  ifastsa=0;
 end
end

istatt=0;
flp=1-iLP1;
mr=1;
Ide=diag(ones(1,2*length(be)));

inuo=0;
lfishav=length(find(shavet(:,1)~=0));
if iany==0 & ianys==0
 lfishav=length(find(shavet(:,1)~=0));
% if lfishav<10
  inuo=1;
% end
end
inuo=1;

if ifast==1
 inuo=2;
end

% devo calcolare le matrici di trasmissione sopra e sotto
%
% sopra

[Tfas,Ffas]=eltime(Tfas,Ffas,-5);

 Li=Litn;
 ni=nitn;
 ai=aitn;
 bi=aral.y.t';
 pai=aral.p.t';
 shai=shav.t';
 istfie=istfi.t;
 tyari=aral.ar.t';
 fmlsi=fmlst;


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

%disp(' re_new')
%keyboard

if igau==4
% disp('rel')
% keyboard
 fmol=fmlsi(:,2);
 fi=find(fmol==0 & real(ni(1,:))'~=1);
 fmol(fi)=ones(size(fi));
 zrif=Li.*fmol;
end



 itutmir0=1;
 if igau==4
  itutmir=itutmir0;
 else
  itutmir=0;
 end
% disp(' sopra '), keyboard

% catena

  if iTsav==0
   if length(pvet)>1
    if pol==pvet(2)
%      ' qui problema ', keyboard
      Tstor=Tstort;
    end
   end
  end



%'prima di Tt', keyboard
if ifp~=-4
disp(' ++++ Top ')
end
 itop=1;
 if iTsav==0
  if ifastsa==0
   if pol==-1
    Tstof=Tstortm;
   else
    Tstof=Tstortp;
   end
  end
 else
   if pol==-1
    nTstof='Tstortm';
   else
    nTstof='Tstortp';
   end
 end
  if ifr>1 | pol==pvetch(2)
   icousav=icousavto;
  end

% pack

 isomv=0;

 chain_i

  if ifr==1 & pol==pvet(1)
   icousavto=icousav;
   clear icousav
  end


 if length(if_only_out)==0
  Tt=Tdu;
 end

 if ifp==-11 | ifp==0
%  'Tt', keyboard
   ' dopo Tt', keyboard
 end

% if pol==pvet(1)
  if iTsav==0
   Tstort=Tstor;
  end
% end

% clear Tdu Tstor
clear Tstor

 if exist('Tmeduf')
  Tmef=Tmeduf;
  clear Tmeduf
 end

% Tmirup=Tmir;

[Tfas,Ffas]=eltime(Tfas,Ffas,-5);

 inv_mixu

 if length(if_only_out)~=0
  Mn=Mon*Tdu;
  Mi=Moi*Tdu;
  if ifp==0
  ' Moi, Mon ', keyboard
  end
 end



%disp(' re_new per velocita')
%keyboard

%%%%   sotto
[Tfas,Ffas]=eltime(Tfas,Ffas,0);

 Li=Lib;
 ni=nib;
 ai=aib;
 bi=aral.y.b';
 pai=aral.p.b';
 shai=shav.b';
 tyari=aral.ar.b';
 istfie=istfi.b;

 fmlsi=fmlsb;
 T=IdeOon;
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

 if iTsav==0
  if length(pvet)>1
   if pol==pvet(2)
    Tstor=Tstorb;
   end
  end
 end

%'prima di Tb', keyboard
if ifp~=-4
disp(' ++++ Bottom')
disp('   ')
end
 itop=0;

  if iTsav==0
   if ifastsa==0
    if pol==-1
     Tstof=Tstorbm;
    else
     Tstof=Tstorbp;
    end
   end
  else
   if pol==-1
    nTstof='Tstorbm';
   else
    nTstof='Tstorbp';
   end
  end


  if ifr>1 | pol==pvetch(2)
   icousav=icousavbo;
  end


% pack
 chain_i


  if ifr==1 & pol==pvet(1)
   icousavbo=icousav;
   clear icousav
  end


% if pol==pvet(1)
  if iTsav==0
   Tstorb=Tstor;
  end
% end

 if length(if_only_out)~=0
  Mn=Tdu*Mn;
  Mi=Tdu*Mi;
 end



% Tb=Tdu;
% 'Tb', keyboard
% clear Tdu Tstor

if ifp==-11 | ifp==0
 'dopo Tb', keyboard
end
 if exist('Tmeduf')
  Tmefb=Tmeduf;
  clear Tmeduf
 end
% catenas

%disp(' re_new')
%keyboard
% Tmirdw=Tmir;

[Tfas,Ffas]=eltime(Tfas,Ffas,0);
%disp('re_new'), keyboard

%   ilayfastto=ilayfastd;

if ifast==1 & iTsav==0
 if ifastsa==1
  fi=find(ilayfastto==1);
  ipsat=icousavto(fi);
  fi=find(ilayfastbo==1);
  ipsab=icousavbo(fi);

    if pol==-1
%    ' salva pol=-1 ', keyboard
     if length(ipsat)>0
      for ip=ipsat
       Tstortm{ip}=Tstort(:,:,ip);
      end
     else
       Tstortm=0;
     end
     if length(ipsab)>0
      for ip=ipsab
       Tstorbm{ip}=Tstorb(:,:,ip);
      end
     else
       Tstorbm=0;
     end
%     save polm Tstortm Tstorbm
%     clear Tstortm Tstorbm
    else
%    ' salva pol=1 ', keyboard
     if length(ipsat)>0
      for ip=ipsat
       Tstortp{ip}=Tstort(:,:,ip);
      end
     else
       Tstortp=0;
     end
     if length(ipsab)>0
      for ip=ipsab
       Tstorbp{ip}=Tstorb(:,:,ip);
      end
     else
       Tstorbp=0;
     end
%   save polp Tstortp Tstorbp
%   clear Tstortp Tstorbp
  end
 end
end
if ifp==-11
 ' qui per verifica ', keyboard
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Ga1=diag(Gas);
  Ga2=diag(Gad);

%% calcolate le matrici di trasmissione sopra e sotto, calcolo l'eq. agli autov.

 clear alphavv Gvav Anz Anzf CCd

   pu1=1:length(Gas);
   pu2=length(Gas)+1:2*length(Gas);

[Tfas,Ffas]=eltime(Tfas,Ffas,-20);

'prima di invpro'
mem_used
%keyboard

%  invpro
%  inv_mix
inv_lll


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

if ifp==-11 | ifp==0
 disp(' dopo inv_mix ')
 keyboard
end

   if isopro>1
    return
   end

  if iraffi==0

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
%      Anz(:,iNc,ics)=Ab.*Trs;
%      Anzf(:,iNc,ics)=Af.*Trs;
      Anz(:,iNc,ics)=Ab;
      Anzf(:,iNc,ics)=Af;

     if length(if_only_out)==0

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
%      Anz(:,iNc,ics)=Ab.*Trd;
%      Anzf(:,iNc,ics)=Af.*Trd;
      Anz(:,iNc,ics)=Ab;
      Anzf(:,iNc,ics)=Af;


%' BAK ', keyboard, keyboard

      if exist('Tmef')
       sdu=size(Tmef);
       sTm=sdu(1:2);
      end


%  sopra
%disp('icsfit'), keyboard

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

     end  %if_only_out

    end   %iNc

%    clear Tb Tt Moi Mon Tmefb Tmef Mn Mi

   end % iraffi

%' BAK ', keyboard, keyboard

clear Ttotal

% raf=nomeFs(1:end-4);
% tinum=now-731239; tim=round(tinum*100000);
% tich=num2str(tim);
% nfsave=[raf tich];
% if exist('Tmef') & exist('Tmefb')
%  eval(['save ',nfsave,' Tb Tt Moi Mon Mi Mn Tmefb Tmef'])
% elseif exist('Tmef')
%  eval(['save ',nfsave,' Tb Tt Moi Mon Mi Mn Tmef'])
% elseif exist('Tmefb')
%  eval(['save ',nfsave,' Tb Tt Moi Mon Mi Mn Tmefb '])
% else
%  eval(['save ',nfsave,' Tb Tt Moi Mon Mi Mn '])
% end

%disp(' save in re_new ')
%keyboard

[Tfas,Ffas]=eltime(Tfas,Ffas,-20);

%nomo=['stime',Ev_or_Od];
nomo=['stime',num2str(mm)];
eval([' save ',nomo,' Tfas Ffas '])
% disp('re_newLLLL plotim'), keyboard
plotim
if ifp==-11
 disp('re_new'), keyboard
end

%keyboard
if ifp>=0, pausak, end

'dopo di invpro'
mem_used
