% 'entro re-new', keyboard
isemp_plan=0;
ideb=1;
if ifr>1
 ideb=0;
end 
% ideb=0;

if isfield(Pf,'isemp_plan')==1
 isemp_plan=Pf.isemp_plan;
end
imatTsto=0
%isemp_plan=0;

ide=1;
ide=0;
isto=0;

decision= ifr>0 & ide==1


if ifr==1 & pol==pvet(1) & ifp==-10
% 'entro re-new', keyboard
end
iatttivo=0;
%save renew
%'*****   emme ',  keyboard
igam=0;  %=1: gamma equivalente dal qw: =2 matrice da invertire
%clear icrit
 if ~exist('ifiez')
  ifiez=0;
 end
 if ifiez==1
  clear Acoz Nz
  icfiez=1;
  ncompar=-1;
  incdiv=1;
  clear zi
  zi(1)=0;
  Nz(1)=nv(1);
if ifp==-10
    ' passo fiez re_new critico', keyboard
end  
%  Amod=[Gas.*Ar; Ar];
%' keyboard crit  re_new', keyboard
if exist('Trcritu')
 if iscell(Trcritu)

  Are=Trcritu{icpo}*Ar;
 else
  Are=Trcritu*Ar;
 end
else
 Are=Ar;
end 
  Amod=[Gas.*Are; Are];  
  Amodi=[Ga1*Ar; Ar];  
  Amod=[Gas.*Are; Are];  
  Acoz(:,icfiez)=Amod;

if ~exist('ficrit')
 ficrit=[];
end

%' keyboard crit  re_new', keyboard 
  if length(ficrit)>0
%   ficri=ficrit;
   ficri=pucrit;    % correzione per tener conto dello specchio con molteplicit�
   Are=Ar;
   lfic=length(ficri);
   if lfic>1
   	if length(ficriall)==1
	    for kc=lfic:-1:2
	     Are=Trc(:,:,kc)*Are;
	     Apr=Grc(:,:,kc)*Are;
	     Amod=[Apr; Are];
	     Acozd(:,kc)=Amod;
	    end 
	    Acozd(:,end+1)=Amodi;
	    for kc=(1:lfic)
	     icfiez=icfiez+1;
	     Acoz(:,icfiez)=Acozd(:,kc+1);
	     Nz(icfiez)=nitn(1,ficri(kc));
	     zi(icfiez)=zi(icfiez-1)+Litn(ficri(kc));
	    end 

    
%'ICRIT per plot fine', keyboard 
    	else
    aggiunta_campi_Scatt   
    
%'ICRIT per plot fine', keyboard  
  	
    	[Gacritf,Trcritf,Trcf,Grcf]=Gam_critScatt_Field(Ga_inf,Mcrit,Layz,LzRatio,dovesonoC);
% 'fine critico VEDOOO' , keyboard          	
            lfic=length(Lz);
	    for kc=lfic:-1:2
	     Are=Trcf(:,:,kc)*Are;
	     Apr=Grcf(:,:,kc)*Are;
	     Amod=[Apr; Are];
	     Acozd(:,kc)=Amod;
	    end 
	    Acozd(:,end+1)=Amodi;
	    for kc=(1:lfic)
	     icfiez=icfiez+1;
	     Acoz(:,icfiez)=Acozd(:,kc+1);
%	     Nz(icfiez)=nitn(1,ficri(kc));
%	     zi(icfiez)=zi(icfiez-1)+Litn(ficri(kc));
	    end 
            Nz=[Nz Refz.'];	
            zi=[zi cumsum(Lz')];
    	end
    

    
   else
     icfiez=icfiez+1;
     kc=1;
     Acoz(:,icfiez)=Amodi;
     Nz(icfiez)=nitn(1,ficri(kc));
     zi(icfiez)=zi(icfiez-1)+Litn(ficri(kc));
   end
   %'quo', keyboard
   

   
   Amod_i=Acoz(:,end);
   %Amod_i=Acoz(:,length(ficri)+1);
   %Amod_i=[Gas.*(Ar); Ar];
   %icfiez=icfiez+1;
   % zi(icfiez)=zi(icfiez-1)+Litn(ficri(kc));   
   %Acoz(:,icfiez)=Amod_i;   
   Amod_ref=[Gacrit*Ar_sav; Ar_sav];
%   Amod=Amod_ref;
   Amod=Amod_i;
  %' Prima Litn fiez re_new critico', keyboard
  
  
   Litn_sav=Litn;
   Litn(ficri)=0;
  end
if ifp==-10
  ' dopo fiez re_new critico', keyboard
end  
%  ' dopo fiez re_new', keyboard

 else
  ioxi=0;
 end
dpes=diag(pes);
iant=0;
%iant=1;
%'passo iant', keyboard
% 
if iant==0
 dpes1=1;
 dpes2=dpes;
else
 dpes2=1;
 dpes1=dpes;
end
% dpes2=diag(sqrt(pes));
% dpes1=diag(sqrt(pes));
% 'diag', keyboard
iTsav=1;
iTsav=0
%keyboard

if length(pvet)==2
 pvetch=pvet;
else
 pvetch=[pvet -100];
end

%ifast=1;
%ifast=0;
global ifast
%'qui ifast', keyboard
isem=1;   %=1 metodo veloce per matrici
%ifast=1;
if ispeed==0
 ifast=0;
 isem=0;   %=1 metodo veloce per matrici
 iTsav=0;   %=1 metodo veloce per matrici
end
%if ispeed==-1
% ' moltoa maostnl attenzione: campbiato ispped '
% ifast=0;
% isem=1;   %=1 metodo veloce per matrici
% iTsav=1;   %=1 metodo veloce per matrici
%end

Mb=2^20;
lfishav=length(find(shavet(:,1)~=0));
sMAT=lfishav*length(Pust)^2*16/Mb;

%Mb=2^10;

if ifr==1 & pol==pvet(1)
 if ifp==-12
%  ifp=-10;
 end
end

%' salva controllo ', keyboard
siMAX=30;
siMAX=5000;  %max MBite
%siMAX=1;  %max MBite
' %siMAX=100;  %max MBite  in re_new '
%siMAX=500;  %max MBite
if isfield(Pf,'isemp_plan')==1
 if Pf.isemp_plan==1
  siMAX=1e6;
 end
end
%siMAX=.1;
global iNOsav
if length(iNOsav)>0
 if iNOsav==1
  siMAX=1e6;
 end
 if iNOsav==-1
   siMAX=.01;
 end
end

global Dsav

iTsav=1;
if sMAT<=siMAX
 iTsav=0;
end
if ifr==1 & pol==pvet(1)
%' controllo iTsav', keyboard
end
%' controllo iTsav', keyboard
if iTsav==1 & ifr==1 & pol==pvet(1)
%' controllo iTsav', keyboard
  ck=clock;
  Dire=num2str(floor(prod(ck( find(ck~=0) ))*1e-4));
  Droo=cd;
  eval(['!md ',Droo,'\',Dire]);
  Dsav=[Droo,'\',Dire];
  iDsav=1;
%  disp(' chain: layskip '), keyboard
%else
% iTsav=0;


end

global iDsav
   if ifr==1 
    %'re_new', keyboard
   end 
  if isfield(Ps,'igraef_new')==1 & length(iDsav)==0
   if Ps.igraef_new==4
    if ifr==1 & pol==pvet(1)
     iDsav=1;
     ck=clock;
     Dire=num2str(floor(prod(ck( find(ck~=0) ))*1e-4));
     Droo=cd;
     eval(['!md ',Droo,'\',Dire]);
     Dsav=[Droo,'\',Dire];
    end
   end
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
%ifastsa=1;

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

inuo=2;
% devo calcolare le matrici di trasmissione sopra e sotto
%
% sopra

[Tfas,Ffas]=eltime(Tfas,Ffas,-5);

 dovesono=1;  %sopra
 pumir=putm;
 Li=Litn;
 ni=nitn;

  global set_perd0
  if length(set_perd0)==0
   set_perd0=0;
  end
    'iper', keyboard
  if set_perd0==1
   nip=ni;
   fin=find(abs(imag(nip))<.1);
%   nip(fin)=nip(fin)-i*imag(nip(fin));
   nip(fin)=real(nip(fin));
   ni=nip;
  end
 ai=aitn;


 icrit=zeros(size(Li));
 if size(iauto,2)==3
  ficriall=find(iauto(:,3)==-10);
  if ficriall==length(iauto)
   ICRI=0;   % scatt totale
  end
   iatti=find(iauto(:,1)==2);
  if ficriall<iatti
     ICRI=1;   % scatt sopra
  end
 else 
  ficriall=[];
 end

 icriso=0;
%  'ICRI 0', keyboard
if isfield(Pf,'ICRIT')==1
 ICRIT=Pf.ICRIT;
else
 ICRIT=0;
end

if ICRIT==0
 ficriall=[];
end
 
 if length(ficriall)>0
  if dovesono==2
%   'CRIT non ancora implementato sotto', 
  ficrip=find(iauto(pub,3)==-10);
  if length(ficrip)==1
   icrit(ficrip:end)=1;
%   'ICRIT', keyboard   
  elseif length(ficrip)==2
   icrit(ficrip(1):end)=1;
  elseif length(ficrip)>2
   'Errore in CRIT: troppi CRIT in file .str', keyboard
  end
   icriso=1;
   'icrit sotto ',  %keyboard  

  else
  ficrip=find(iauto(put,3)==-10);
  if length(ficrip)==1
   icrit(1:ficrip)=1;
%'ICRIT', keyboard   
  elseif length(ficrip)==2
   pucri=1:ficrip(1);
   icrit(pucri)=1;
   fmlsdu=fmlst;
% calcolo tutti gli z da fare con Gam_crit

    
   
  elseif length(ficrip)>2
   'Errore in CRIT: troppi CRIT in file .str', keyboard
  end
   icriso=1;
   'icrit sopra ',  %keyboard  
  end


 end

 %' dopo icri', keyboard
 
 anir=anyr.t;
 bi=aral.y.t.';
 pai=aral.p.t.';
 shai=shav.t.';
 istfie=istfi.t;
 prosi=prosn;

 tyari=aral.ar.t.';
 fmlsi=fmlst;



 vfia=(shai(1,:)'.*ai(1,:)').*Li.*(1-anir)*isemp_plan;
 
shaiP=shai;
LiP=Li;
aiP=ai;
anirP=anir;
 
 fia=find(vfia~=0)';
 
if length(fia)>0
  iloo=fia;
  if fia(1)~=1 & Li(1)>0
   ilins=[0 fia];
  else
   ilins=[fia];
 end
 
% 'qui', keyboard
% 'qui', keyboard
else
 ilins=0;

end  %length(iloadd)

'dopo ifast', keyboard

if ifr==1
 if ideb==1
  ' ai sopra ', keyboard
 end 
end

 if iany==0
  Kan=0;
  ani=zeros(size(Li));
 elseif iany==1 | iany==3
  Kan=Kosan;
  ani=deltanyn*fattany;
 elseif iany==2
  Kan=Kosan*diag(pes);
  ani=deltanyn*fattany;
 end


if ~exist('Kan_gr')
 Kan_gr=0;
end 
 if iLP==0
 if exist('i_grap')
  if i_grap==1
   Kan_gr=Kanr;
  elseif i_grap==2
   Kan_gr=Kan_lim*diag(pes);
  end
 end
 end

%' Gattany', keyboard

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
 %fi=find(fmol==0 & real(ni(1,:))'~=1);
 fi=find(fmol==0);
 fmol(fi)=ones(size(fi));
 zrif=Li.*fmol;
end



 itutmir0=0;
 if igau==4
  itutmir=itutmir0;
 else
  itutmir=0;
 end
% se itutmir=1 analizza la Temp di tutto lo specchio, altrimenti no


 %disp(' sopra '), keyboard

% catena

  if iTsav==0
   if length(pvet)>1
    if pol==pvet(2)
%      ' qui problema ', keyboard
      Tstor=Tstort;
    end
   end
  end
if ifastsa==0
%      Tstor=Tstort;
end

%'prima di Tt', keyboard
if ifp~=-4
disp(' ++++ Top ')
end
IFPS=ifp;
ifp=-4;
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
   if isemp_plan==1
    istrv=istrvto;
   end 
  end

% pack

 isomv=0;
  Ga1=diag(Gas);
  Ga2=diag(Gad);
% 'Sopra ai ', keyboard  
if ifp==-10
% 'Sopra ai ', keyboard
end
if decision
 if ispeed==0
  save zero_1
 else
  save uno_1
 end
 ifr
 'salvato'
end


zsopra=0;

%'% parte per fare veloce gli strati piani, Sett 2015', keyboard

if ~exist('iloo')
 chain_i
else
 %'inizio parte semplificata fastu', keyboard
 chain_fastuu
 
 %'inizio parte semplificata', keyboard
 if ifr==1 & pol==pvet(1)
%  if ilins(1)~=1
%   ilins=[1 ilins];
%  end
% if ifr==1
  fipla=ilins(find(diff([ilins 0])~=1));
  ilinsi=ilins(find(diff([ilins 0])~=1))+1;
  ilinsu=[ilins(find(diff([ilins])~=1)+1)-1 length(Li)];
%  ilinsi=fipla+1;
%  ilinsu=[fipla(2:end)-1 length(Li)];
%  ilinsi=ilins+1;
%  ilinsu=[ilins(2:end)-1 length(Li)];
%   'inizio parte piana', keyboard
  Pstack.in=ilinsi; 
  Pstack.fin=ilinsu;
  Pstack.n=ni(1,1:end);
  Pstack.Li=Li;
  Pstack.rep=fmlsi(:,2);
  Pstack.nmo=numodiacc+1;
  Pstack.KK=KK;
  Pstack.rr=rr;
  Pstack.pol=pol;
  Pstack.Pus0=Pus0;
 end 
 
 if pol==pvet(1)
  if imatTsto==1
   Tstack=plan_stack0(kcav,Pstack); 
   ' verifica Matrici', keyboard

  else
%  'plan', keyboard
   Tstack=plan_stack(kcav,Pstack); 
  end 
 end
 
%   ' verifica Matrici', keyboard
 Tdu=1;
 ilopl=0;
%  if ilins(1)==0
  if fipla(1)==0
   ilopl=ilopl+1;
   Tdu=Tstack(:,:,ilopl);
  end
  %keyboard
  
 if numodiacc<3 
% ' dopo Sopra Tdu ', keyboard
 end
 
%  for klo=1:length(iloo)
%%  klo, pausak
%   Tdu=Tstor(:,:,klo)*Tdu;
%   if length(find(iloo(klo)==fipla))==1
%    ilopl=ilopl+1;
%    Tdu=Tstack(:,:,ilopl)*Tdu;
%   end
%  end  
    ilopl=0;
    Tdu=1;
 if iloo(1)>1
    filo=1:iloo(1)-1;
    if sum(Li(filo))>0
     ilopl=1;
     Tdu=Tstack(:,:,ilopl);
    end 
 end   
 klo=0; 
 rep=fmlsi(iloo,2);
 nstmir=fmlsi(iloo,1);
 while klo<length(iloo)
  klo=klo+1;
%  ind=iloo(klo);
  nsti=nstmir(klo);
%  pausak
   if nsti==1
%  klo, pausak
    Tdu=Tstor(:,:,klo)*Tdu;
   else
%     klo, pausak
    Tmir=1;
    NPair=rep(klo);
    pumir=[klo:klo+nsti-1];
    klos=klo;
    for pum=pumir
%     ind=iloo(pum);
     ind=pum;
     Tmir=Tstor(:,:,ind)*Tmir;
    end
    Tmir=Tmir^NPair;
    klo=klo+nsti-1;
    Tdu=Tmir*Tdu;
%    'mir', keyboard
   end
   if length(find(iloo(klo)==fipla))==1
    ilopl=ilopl+1;
    Tdu=Tstack(:,:,ilopl)*Tdu;
   end
  end   
  
 if ifr==1 & pol==pvet(1)
%  ' dentro semplificato ', keyboard
 end 
 
%   ' dentro semplificato fine', keyboard
end

if ifr==1 & pol==pvet(1)
 if ideb==1
  if numodiacc<3 
%  ' dopo Sopra chain ', keyboard
  end
 end
end
 
if igau==4
 zsopra=sum(zrif);
end


  if ifr==1 & pol==pvet(1)
   icousavto=icousav;
   clear icousav
    if isemp_plan==1 & exist('istrv')
       istrvto=istrv;
       clear istrv
    end
  end
 Tt=Tdu;

%   ' dopo Tt', keyboard
   if ideb==1
      ' dopo Tt', keyboard
    pausak
   end
   
   %'Dopo Tt', keyboard
   if isto==1, keyboard, end
 if ifp==-11 | ifp==0 | ifp==-12
%  'Tt', keyboard
   ' dopo Tt', keyboard
 end
%if ifr==1 & pol==-1
%if  pol==-1
%'ifr=  ', ifr
%   ' dopo Tt', keyboard
%end   

%  if (igam==1 | ifr==1) & iLP==0
  if (igam==1 | ifr==1) & pol==pvet(1) & ifp==-100
%   ' dopo Tt', keyboard
   iproga=1;
  if exist('ldapuu')
   King=KKt(Pusas(1:ldapuu(2)));
  else
   King=KKt(Pusas(1:ldapu(2)));
  end
  if length(Pusas)>130
   istop=0;
  else
   istop=1;
  end
  if ifp>=-10
   istop=0;
  end
  prif(Tt,iauto,dv,nv,fmlstot,rr,lambda,freq,King,Gas,iLP,istop,iem_like,fapeu)
%   if igam==1
%    keyboard
%   end
  end
%   keyboard

% if pol==pvet(1)
  if iTsav==0
   Tstort=Tstor;
  end
% end
if decision
' qui Tstort', keyboard
end

 clear Tdu Tstor

 if exist('Tmeduf')
  Tmef=Tmeduf;
  clear Tmeduf
 end

% Tmirup=Tmir;

[Tfas,Ffas]=eltime(Tfas,Ffas,-5);

%disp(' re_new per velocita')
%keyboard

if ifiez==1
%    ' qui zona attiva PRIMA ', keyboard
  if isav_Az==1
   inv_mix
  end
%    ' qui zona attiva ', keyboard
  icfiez=icfiez+1;
  Acoz(:,icfiez)=(Mons+Mois*Lami)*Acoz(:,icfiez-1);
  Aqw=Acoz(:,icfiez);
  Amod=Acoz(:,icfiez);
  Nz(icfiez)=ra;
  zi(icfiez)=zi(icfiez-1)+d*1e6;
 if ifp==-10
  ' qui per ultimo ', keyboard
 end 
end


%%%%   sotto
[Tfas,Ffas]=eltime(Tfas,Ffas,0);

 pumir=pubm;
 dovesono=2;  %sop1, sotto=2
 Li=Lib;
 ni=nib;
 ai=aib;
  if set_perd0==1
   nip=ni;
   fin=find(abs(imag(nip))<.1);
%   nip(fin)=nip(fin)-i*imag(nip(fin));
   nip(fin)=real(nip(fin));
   ni=nip;
  end 
 
 icrit=zeros(size(Li)); 
 
% 'ICRI 1', keyboard
 
 
  icriso=0;
  if length(ficriall)>0
   if dovesono==2
%    'CRIT non ancora implementato sotto', 
   ficrip=find(iauto(pub,3)==-10);
   if length(ficrip)==1
    Ldu=Li;
    Ldu(1:ficrip-1)=0;
    pucri=find(Ldu>0);
%    icrit(ficrip:end)=1;
    icrit(pucri)=1;
 %'ICRIT', keyboard   
   elseif length(ficrip)==2
    icrit(ficrip(1):end)=1;
   elseif length(ficrip)>2
    'Errore in CRIT: troppi CRIT in file .str', keyboard
   end
    icriso=1;
    'icrit sotto ',  %keyboard  
 
   else
   ficrip=find(iauto(put,3)==-10);
   if length(ficrip)==1
    icrit(1:ficrip)=1;
 %'ICRIT', keyboard   
   elseif length(ficrip)==2
    icrit(1:ficrip(2))=1;
   elseif length(ficrip)>2
    'Errore in CRIT: troppi CRIT in file .str', keyboard
   end
    icriso=1;
    'icrit sopra ',  %keyboard  
   end
 
 
  end
 
 %' dopo icri', keyboard
 
 prosi=prosb;
 anir=anyr.b;
 bi=aral.y.b.';
 pai=aral.p.b.';
 shai=shav.b.';
 fiai=find(ai(1,:)==0);
 shai(1,fiai)=0;
% 'qui shai', keyboard
 
 if isfield(aral.ar,'b')==1
  tyari=aral.ar.b.';
 else
  tyari=0;
 end

 istfie=istfi.b;

 fmlsi=fmlsb;
 T=1;
 ilaymem=find(diff(ai)~=0);
%' ilyamem', keyboard

 if iany==0
  Kan=0;
  ani=zeros(size(Li));
 elseif iany==1 | iany==3 
  Kan=Kosan;
  ani=deltanyb*fattany;
 elseif iany==2
  Kan=Kosan*diag(pes);
  ani=deltanyb*fattany;
 end

%  ' secondo Kan ', keyboard
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
iback_plan=1;
if length(Li)>0
iback_plan=0;
end

if iback_plan==0
 if iTsav==0
  if length(pvet)>1
   if pol==pvet(2)
    Tstor=Tstorb;
   end
  end
 end
if ifastsa==0
%      Tstor=Tstorb;
end
% ' qui Tstortb' , keyboard
%clear icrit

%'prima di Tb', keyboard
if ifp~=-4
disp(' ++++ Bottom')
disp('   ')
end
 itop=0;
ifp=IFPS;
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
   if isemp_plan==1
    istrv=istrvbo;
   end 
  end


if ifiez==1
  if length(ficrib)>0
   ficri=ficrib;
%'qui BOT Li', keyboard
   Li(ficri)=0;
  end
end

% pack
 chain_i

  if ifr==1 & pol==pvet(1)
   icousavbo=icousav;
   clear icousav
   if isemp_plan==1
    istrvbo=istrv;   
    clear istrv
   end 
  end




% if pol==pvet(1)
  if iTsav==0
   Tstorb=Tstor;
  end
% end

 Tb=Tdu;
 Tbfield=Tb;
else
 Tb=1;
 Gad_refresh
 Dt=diag(Trd);
 Ze=zeros(size(Dt));
 Tbfield=[Dt Ze; Ze Dt];
end

   ' dopo Tb'
   if ideb==1
    keyboard
   end 
    
   
  if Ps.ifpstop==1

 'dopo Tb', 
 figure, plot(zi,abs(Acoz))

 keyboard
 end

 isalta=1;
  if (igam==1 | ifr==1) & pol==pvet(1) & ifp==-10 & isalta==0
   ' dopo Tb', keyboard
  King=KKt(Pusas(1:ldapuu(2)));
  if length(Pusas)>130
   istop=0;
  else
   istop=1;
  end
  if ifp>=-10
   istop=0;
  end
  prifb(Tb,iauto,dv,nv,fmlstot,rr,lambda,freq,King,Gad,iLP,istop,iem_like,fapes,nstratid)
%   if igam==1
%    keyboard
%   end
  end


if ifiez==1
  if length(ficrib)>0

   Are=Acoz(end/2+1:end,end);
   Are0=Are;
%   'qui BOT campo', keyboard
   lfic=length(ficri);
   if lfic>1
    ckc=0;
    clear Acozd
    for kc=ficri'
     Are=Trc(:,:,kc)*Are;
     Apr=Grc(:,:,kc)*Are;
     Amod=[Apr; Are];
     ckc=ckc+1;
     Acozd(:,ckc)=Amod;
    end 
    for kc=(1:lfic)
     icfiez=icfiez+1;
     Acoz(:,icfiez)=Acozd(:,kc);
     Nz(icfiez)=ni(1,ficri(kc));
     zi(icfiez)=zi(icfiez-1)+Lib(ficri(kc));
    end 
   end
%   'quo', keyboard
   
  end
end

 %'Tb', keyboard
 clear Tdu Tstor
 if ifiez==1
  Aend=Tbfield*Aqw;
%  ' chiamo dis_fz in re_new ', keyboard
%  dis_fz
  FF.Acoz=Acoz;
  FF.zi=zi;
  FF.Nz=Nz;
  ' fine re_new per calcolo E(z) '
%  keyboard
  return
 end

% 'dopo Tb', keyboard
if ifp==-11 | ifp==0
 'dopo Tb', keyboard
end
 if exist('Tmeduf')
  Tmefb=Tmeduf;
  clear Tmeduf
 end
% catenas

% disp(' re_new'), keyboard
% Tmirdw=Tmir;

[Tfas,Ffas]=eltime(Tfas,Ffas,0);
%disp('re_new'), keyboard

%   ilayfastto=ilayfastd;

if ifast==1 & iTsav==0
% if ifastsa==1
  if isemp_plan==0
   fi=find(ilayfastto==1);
  else
   fi=1:length(icousavto);
  end
  ipsat=icousavto(fi);
  if exist('iloo')
   ipsat=1:length(iloo);
  end
  if iback_plan==0
   fi=find(ilayfastbo==1);
   ipsab=icousavbo(fi);
  end 

    if pol==-1
%    ' salva pol=-1 ', keyboard
     if length(ipsat)>0
      for ip=ipsat
       Tstortm{ip}=Tstort(:,:,ip);
      end
     else
       Tstortm=0;
     end

     if iback_plan==0
     if length(ipsab)>0
      for ip=ipsab
       Tstorbm{ip}=Tstorb(:,:,ip);
      end
     else
            Tstorbm=0;
     end
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




     if iback_plan==0
     if length(ipsab)>0
      for ip=ipsab
       Tstorbp{ip}=Tstorb(:,:,ip);
      end
     else
       Tstorbp=0;
     end
     end
%   save polp Tstortp Tstorbp
%   clear Tstortp Tstorbp
  end
% end
end
if ifp==-11
 ' qui per verifica ', keyboard
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% calcolate le matrici di trasmissione sopra e sotto, calcolo l'eq. agli autov.

 clear alphavv Gvav Anz Anzf CCd

   pu1=1:length(Gas);
   pu2=length(Gas)+1:2*length(Gas);

[Tfas,Ffas]=eltime(Tfas,Ffas,-20);
  
  
  inv_mix
  %inv_mix_oggi
  


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

%if ifiez==0 & ired_ret==0
if ifiez==1000 & ired_ret==0
%nomo=['stime',Ev_or_Od];
nomo=['stime',num2str(mm)];
eval([' save ',nomo,' Tfas Ffas '])
% disp('re_new plotim'), keyboard
plotim
if ifp==-11
 disp('re_new'), keyboard
end
end



%keyboard
if ifp>=0, pausak, end
