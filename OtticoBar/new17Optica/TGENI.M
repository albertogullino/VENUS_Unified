'*****   tgeni '
 if ~exist('ifiez')
  ifiez=0;
 end
icfiez=1;
Amod=[Gas.*Ar; Ar];
Acoz(:,icfiez)=Amod;

ifast=1;
isem=1;   %=1 metodo veloce per matrici

if ispeed==0
 ifast=0;
 isem=0;   %=1 metodo veloce per matrici
 iTsav=0;   %=1 metodo veloce per matrici
end


lfishav=length(find(shavet(:,1)~=0));
sMAT=lfishav*length(Pust)^2*16;
Mb=2^20;

if ifr==1 & pol==pvet(1)
 if ifp==-12
  ifp=-10;
 end
end

%' salva controllo ', keyboard
siMAX=50;

if sMAT<=Mb*siMAX
 iTsav=0;
end

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
 anir=anyr.t;


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

% chain_fi
 cha_fi1


' dopo Top ', keyboard

 clear Tdu Tstor


 if exist('Tmeduf')
  Tmef=Tmeduf;
  clear Tmeduf
 end

% Tmirup=Tmir;

[Tfas,Ffas]=eltime(Tfas,Ffas,-5);

%disp(' re_new per velocita')
%keyboard

%  inv_mix
  icfiez=icfiez+1;
  Acoz(:,icfiez)=(Mons+Mois*Lami)*Acoz(:,icfiez-1);
  Aqw=Acoz(:,icfiez);
  Amod=Acoz(:,icfiez);
  zi(icfiez)=zi(icfiez-1)+d*1e6;

%  ' ifiez dopo invmix', pausak


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
 anir=anyr.b;

 fmlsi=fmlsb;
 T=1;
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
% chain_fi
 cha_fi1

iold=0;
if iold==1
 lsum=0;
 Lth=L_i;
 Lth(1)=0;
 sA=size(Acoz);
 for ks=1:sA(2)
  lsum=lsum+Lth(ks);
  Lu(ks)=lsum;
  Exdu=besm'*(Acoz(pu1,ks)+Acoz(pu2,ks));
  Esdu=abs(Exdu/max(Exdu)).^2;
  Es(:,ks)=abs(Exdu).^2;
%  Es(:,ks)=Esdu;
 end
end
  dis_fz



' dopo bottom ', keyboard
  return



  if ifr==1 & pol==pvet(1)
   icousavbo=icousav;
   clear icousav
  end


% if pol==pvet(1)
  if iTsav==0
   Tstorb=Tstor;
  end
% end

 Tb=Tdu;
% 'Tb', keyboard
 clear Tdu Tstor
 if ifiez==1
  Aend=Tb*Aqw;
 end

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
  inv_mix


clear Ttotal T

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
% disp('re_new plotim'), keyboard
plotim
if ifp==-11
 disp('re_new'), keyboard
end

%keyboard
if ifp>=0, pausak, end
