global Ps
if isfield(Ps,'igraef_new')==1
 igraef_new=Ps.igraef_new;
else
 igraef_new=-1;
end
%' entro reass /qui ', keyboard
% ishamem=find(ifield<0);
% sad=[1 0 -1];
% for k=1:3
%  fiir(k)=find(iauto(:,1)==k)+sad(k);
% end
%
% cfi=0;
% if length(ishamem)>0
%  for k=1:length(ishamem)
%   fi=find(fiir==ishamem(k));
%   if length(fi)==0
%    cfi=cfi+1;
%    fime(cfi)=ishamem(k);
%   end
%  end
%  if cfi>=1
%   istfiv=fime-1;
%  end
% end

nv=nv0;


  global set_perd0
  if length(set_perd0)==0
   set_perd0=0;
  end
  
  %'iper', keyboard
  if set_perd0==1
   nip=nv;
   fin=find(abs(imag(nip))<.1);
%   nip(fin)=nip(fin)-i*imag(nip(fin));
   nip(fin)=real(nip(fin));
   nv=nip;
  end 
   
%' nv', keyboard

rad=radii.a;

dvmi=dv/1000;


%'QUI reassign', keyboard
nmir.a=-sum(fst((find(fst(:,2)<0)),2));
iput=find(iauto(:,1)==1);
ipua=find(iauto(:,1)==2);
ipub=find(iauto(:,1)==3);

ipudd=find(iauto(:,2)==-4 | iauto(:,2)>0);
%fiC=[ipudd(1):ipudd(length(ipudd))]-1;
%fiQ=find(fst(:,2)==-1)-1;
fi1=find(iauto(:,1)==1)-1;
fiCdu=find(iauto(:,2)==-4)-fi1;
fiC=fiCdu(1):fiCdu(end);
%fiQ=find(fst(:,2)==-1);
%fim4=find(iauto(:,2)==-4)-fi1;
%ficavit=fim4-fi1;
%pcavit=ficavit(1):ficavit(2);
%'fiC',keyboard

put=iput+1:ipua-1;
putm=iput:fiC(1)-1;

pub=ipua+1:ipub-1;
pubm=fiC(end)+1:ipub-1;

pua=ipua;

%' qui', keyboard

% ifimem=find(ifield~=0);
 ifimem=find(ifield==-2);
 istfiv=zeros(size(ifield));
% istfiv(ifimem)=ifimem-1-iput;
% ifip=find(ifield==-2)-1-iput;
 istfiv(ifimem)=ifimem-iput;
 ifip=find(ifield==-2)-iput;
 if length(ifip)>0
  if ifip>0
   istfiv(ifip)=-istfiv(ifip);
  end
 end
 
if Ps.ifpstop==1    
% disp(' verifica'), keyboard
% disp(' verifica'), keyboard
end 

put_pl=2:ipua-1;
pub_pl=ipua+1:length(iauto)-1;

%'% setto variabili per struttura completamente planare'
%keyboard
%
% uscita
dt=dv(put_pl,1)/1000;
nt=nv(put_pl,1);
ft=fst(put_pl,:);
if exist('any')==1
 if iany==0 | iany==3
  Po=anyf(put_pl);
 else
  Po=any.t;
 end
 if length(Po)<length(nt)
  Po=[Po; zeros(length(nt)-length(Po),1)];
 end
end

fdu=find(diff(ft(:,1))~=0);
%'fdu', keyboard

if length(fdu)>=2
 ftdu=ft(:,2);
 nsmax=max(ftdu);
 fimax=find(ftdu==nsmax);
 ft1=ft(:,1);
 ft2=abs(ft(:,2));
 ft1(fimax)=1;
 ptn=expand(ft1,ft2,nsmax);
 
% 'expand', keyboard
 ftn=ftdu(ptn);
 fispe=find(ftn==nsmax);

 fisped=ptn';
 pu1d=1:fispe(1)-1;
 pure1=fisped(pu1d);
 pure1=fliplr(pure1);

 Lpla.t.o=dt(pure1);
 npla.t.o=nt(pure1);
 Nspla.t.o=0;
 if exist('Po')==1
  Poc.t.o=Po(pure1);
 end

 pu1d=fispe;
 pure1=fisped(pu1d);
 pure1=fliplr(pure1);
 Lpla.t.m=dt(pure1);
 npla.t.m=nt(pure1);
 Nspla.t.m=ft(pure1(1),2);
 if exist('Po')==1
  Poc.t.m=Po(pure1);
 end

 pu1d=fispe(end)+1:length(fisped);
 pure1=fisped(pu1d);
 pure1=fliplr(pure1);
 Lpla.t.i=dt(pure1);
 npla.t.i=nt(pure1);
 Nspla.t.i=0;
 if exist('Po')==1
  Poc.t.i=Po(pure1);
 end

%' Poc mod'
%keyboard

elseif length(fdu)==1
 pure1=1:fdu(1);
 pure1=fliplr(pure1);
 Lpla.t.m=dt(pure1);
 npla.t.m=nt(pure1);
 Nspla.t.m=ft(pure1(1),2);
 if exist('Po')==1
  Poc.t.m=Po(pure1);
 end

 pure1=fdu(1)+1:length(nt);
 pure1=fliplr(pure1);
 Lpla.t.i=dt(pure1);
 npla.t.i=nt(pure1);
 Nspla.t.i=0;
 if exist('Po')==1
  Poc.t.i=Po(pure1);
 end
 Lpla.t.o=[];
 npla.t.o=[];
 Nspla.t.o=[];
 Poc.t.o=[];
else
 Lpla.t.m=[];
 npla.t.m=[];
 Nspla.t.m=0;
 if exist('Po')==1
  Poc.t.m=[];
 end

 fii=find(nt==1);
 if length(fii)>0
  fiil=fii(end);
  pure1=fiil+1:length(nt);
  pure1=fliplr(pure1);
  Lpla.t.i=dt(pure1);
  npla.t.i=nt(pure1);
  if exist('Po')==1
   Poc.t.i=Po(pure1);
  end
 else
  Lpla.t.i=[];
  npla.t.i=[];
 end
  Nspla.t.i=0;

 Lpla.t.o=[];
 npla.t.o=[];
 Nspla.t.o=[];
 Poc.t.o=[];
end

% sotto
dt=dv(pub_pl,1)/1000;
nt=nv(pub_pl,1);
ft=fst(pub_pl,:);
%' qui '
%keyboard
fdu=find(diff(ft(:,1))~=0);

%'fdu back', keyboard

if exist('any')==1
 if iany==0 | iany==3
  Po=anyf(pub_pl);
 else
  Po=any.b;
 end
 if length(Po)<length(nt)
  Po=[Po; zeros(length(nt)-length(Po),1)];
 end
end
%' reassing ', keyboard

if length(fdu)>=2
 inu=1;  %nuovo metodo per trattare specchi con diversi tipi di paia.
 if inu==1
%%%%%%%%%

 ftdu=ft(:,2);
 nsmax=max(ftdu);
 fimax=find(ftdu==nsmax);
 ft1=ft(:,1);
 ft2=abs(ft(:,2));
 ft1(fimax)=1;
 ptn=expand(ft1,ft2,nsmax);
 ftn=ftdu(ptn);
 fispe=find(ftn==nsmax);

 fisped=ptn';
 pu1d=1:fispe(1)-1;
 pure1=fisped(pu1d);
% pu1=fliplr(pu1);

%%%%%%%%%%
 else



 nsmax=max(ft(:,2));
 fispe1=find(ft(:,2)<nsmax & ft(:,2)>1);
 ftdu=ft(:,2);
 if length(fispe1)>0
  fimul=[];
  for ik=1:ft(fispe1(1),2)
   fimul=[fimul; fispe1];
  end
  fisped=[[1:fispe1(1)-1] fimul' [fispe1(end)+1:length(ft)]];
  ftdu=ft(fisped,2);
  fispe=find(ftdu==nsmax);
 else
  fispe=find(ft(:,2)==nsmax);
  fisped=1:length(ft);
 end
 pure1=1:fispe(1)-1;
 end  %inu

 Lpla.b.i=dt(pure1);
 npla.b.i=nt(pure1);
 Nspla.b.i=0;
 if exist('Po')==1
  Poc.b.i=Po(pure1);
 end

 if inu==0
  pure1=fispe;
 else
  pu1d=fispe;
  pure1=fisped(pu1d);
 end

 Lpla.b.m=dt(pure1);
 npla.b.m=nt(pure1);
 Nspla.b.m=ft(pure1(1),2);
 if exist('Po')==1
  Poc.b.m=Po(pure1);
 end
 if inu==1
  pu1d=fispe(end)+1:length(fisped);
  pure1=fisped(pu1d);
 else
  pure1=fispe(end)+1:length(nt);
 end
 Lpla.b.o=dt(pure1);
 npla.b.o=nt(pure1);
 Nspla.b.o=0;
 if exist('Po')==1
  Poc.b.o=Po(pure1);
 end

elseif length(fdu)==1

 pure1=fdu(1)+1:length(nt);
 Lpla.b.m=dt(pure1);
 npla.b.m=nt(pure1);
 Nspla.b.m=ft(pure1(1),2);
 if exist('Po')==1
  Poc.b.m=Po(pure1);
 end
 pure1=1:fdu(1);
 Lpla.b.i=dt(pure1);
 npla.b.i=nt(pure1);
 Nspla.b.i=0;
 if exist('Po')==1
  Poc.b.i=Po(pure1);
 end

%'asso', keyboard
 Lpla.b.o=[];
 npla.b.o=[];
 Nspla.b.o=[];
 Poc.b.o=[];

else
 if length(fdu)~=0
  pure1=fdu(1)+1:length(nt);
  Lpla.b.m=dt(pure1);
  npla.b.m=nt(pure1);
  Nspla.b.m=ft(pure1(1),2);
  if exist('Po')==1
   Poc.b.m=Po(pure1);
  end

  pure1=1:fdu(1);
  Lpla.b.i=dt(pure1);
  npla.b.i=nt(pure1);
  Nspla.b.i=0;
  if exist('Po')==1
   Poc.b.i=Po(pure1);
  end
 else
  Lpla.b.i=[];
  npla.b.i=[];
  Nspla.b.i=[];
  Poc.b.i=[]; 
    Lpla.b.m=[];
    npla.b.m=[];
    Nspla.b.m=[];
  Poc.b.m=[]; 
 end
 Lpla.b.o=[];
 npla.b.o=[];
 Nspla.b.o=[];
 Poc.b.o=[];
end





pudu=1:iput;
if length(pudu)>1
 pupt=find(abs(iauto(pudu,2))==3);
 pup.t.last=pudu(pupt(1));
 pup.t.o=fliplr(pudu(pupt(2:length(pupt))));
 pupt=find(abs(iauto(pudu,2))==2);
 pup.t.m=fliplr(pudu(pupt));
 if length(pupt)>0
  imir=find(abs(iauto(pudu,2))==2);
  nmir.t=fst(pudu(imir(1)),2);
 else
  nmir.t=0;
 end
 pupt=find(abs(iauto(pudu,2))==1);
 pup.t.i=fliplr(pudu(pupt));
else
clear pup
 pup.t.last=pudu(1);
 pup.t.o=[];
 pup.t.m=[];
 nmir.t=0;
 pup.t.i=[];
end

if length(dv)~=ipub
 pudu=ipub:length(dv);
 pupt=find(abs(iauto(pudu,2))==3);
else
 pudu=ipub;
 pupt=1;
end

if length(pupt)>0
 pup.b.o=pudu(pupt(1:length(pupt)-1));
 pup.b.last=pudu(pupt(length(pupt)));
else
 pup.b.o=[];
 pup.b.last=[];
end

pupt=find(abs(iauto(pudu,2))==2);
if length(pupt)>0
 pup.b.m=pudu(pupt);
 imir=find(abs(iauto(pudu,2))==2);
 nmir.b=fst(pudu(imir(1)),2);
else
 pup.b.m=[];
 nmir.b=0;
end
pupt=find(abs(iauto(pudu,2))==1);
if length(pupt)>0
 pup.b.i=pudu(pupt);
else
 pup.b.i=[];
end


Lr.t=dvmi(put);
istfi.t=istfiv(put);


nr.t=nv(put,:);
xr.t=xm(put,:);
ar.t=rad(put,:);
frp.t=fst(put,:);
shav.t=shavet(put,:);
aral.y.t=radii.b(put,:);
aral.p.t=radii.c(put,:);
fisha=find(shav.t==6);
if length(fisha)>=0
 for kf=fisha'
  if igraef_new~=2
   aral.p.t(kf)=radii.array{put(kf)}{1};
  end 
 end
end

icop=0;

sar=size(radii.array);
 if length(sar)==2
  ltc=sar(2);
 else
  ltc=1;
 end

for duput=put
 icop=icop+1;
 for ktc=1:ltc
  du=radii.array{duput,ktc};
  if length(du)==0
   du=0;
  else
   du=radii.array{duput,ktc}{1};
   if length(du)==0
    du=0;
   end
  end
  aral.ar.t(icop,ktc)=du;
 end
end
%'reass', keyboard

Lr.b=dvmi(pub);
istfidu=istfiv(pub);
fidu=find(istfidu~=0);
istfi.b=zeros(length(pub),1);
istfi.b(fidu)=abs(istfidu(fidu))-length(put)-1;
zfidu=find(istfidu<0);
if length(zfidu)>0
 istfi.b(zfidu)=-istfi.b(zfidu);
end
nr.b=nv(pub,:);
xr.b=xm(pub,:);
ar.b=rad(pub,:);

%'reas'
%keyboard
frp.b=fst(pub,:);
shav.b=shavet(pub,:);
aral.y.b=radii.b(pub,:);
aral.p.b=radii.c(pub,:);
fisha=find(shav.b==6);
if length(fisha)>=0
 for kf=fisha'
  if igraef_new~=2
   aral.p.b(kf)=radii.array{pub(kf)}{1};
  end 
 end
end

icop=0;
%for duput=pub
% icop=icop+1;
% aral.ar.b(icop,1)=radii.array{duput,1};
%end
for duput=pub
 icop=icop+1;
 for ktc=1:ltc
  du=radii.array{duput,ktc};
  if length(du)==0
   du=0;
  else
   du=radii.array{duput,ktc}{1};
   if length(du)==0
    du=0;
   end
  end
  aral.ar.b(icop,ktc)=du;
 end
end

Lr.a=dvmi(pua);
nr.a=nv(pua,:);
xr.a=xm(pua,:);
ar.a=rad(pua,:);
%'cor', keyboard
if radii.b(pua,:)>0
 ar.cor=mean([rad(pua,:) radii.b(pua,:)]);
else
 ar.cor=rad(pua,:);
end
frp.a=fst(pua,:);
shav.a=shavet(pua,:);
aral.y.a=radii.b(pua,:);
aral.p.a=radii.c(pua,:);
%aral.ar.a=radii.array(pua,:);
%aral.ar.a=radii.array{pua,1};
for duput=pua
 for ktc=1:ltc
  du=radii.array{duput,ktc};
  if length(du)==0
   du=0;
  else
   du=radii.array{duput,ktc}{1};
   if length(du)==0
    du=0;
   end
  end

  aral.ar.a(duput,ktc)=du;
 end
end

%' qui plan ', keyboard


Lpl.t.o=dvmi(pup.t.o);
Lpl.t.m=dvmi(pup.t.m);
Lpl.t.i=dvmi(pup.t.i);

Lpl.b.o=dvmi(pup.b.o);
Lpl.b.m=dvmi(pup.b.m);
Lpl.b.i=dvmi(pup.b.i);

npl.t.last=nv(pup.t.last,1);
npl.t.o=nv(pup.t.o,1);
npl.t.m=nv(pup.t.m,1);
npl.t.i=nv(pup.t.i,1);

npl.b.last=nv(pup.b.last,1);
npl.b.o=nv(pup.b.o,1);
npl.b.m=nv(pup.b.m,1);
npl.b.i=nv(pup.b.i,1);




 if imet==1        %imet settato in set_stru
  sip=size(nr.t);
  righe=ones(sip);
%  fimet=find(imag(nr.t)<-.1 & righe<10);
%fimet=[];
 if sip(2)>1
%  fimet=find(imag(nr.t(:,2))<-.1);
  fimet=find(imag(nr.t(:,2))<-100);
 else
  fimet=[];
 end

%  if length(fimet)>20
  if length(fimet)>0
%   istoxd=righe(fimet);
%   istoxdu=max(istoxd);
%   istmet0=1;
%   fif=find(istfi.t==istoxdu);
%   if length(fif)==0
%    istfi.t(istoxdu)=istoxdu;
%   end
%   fif=find(istfi.t~=0);
%   istfi.t=istfi.t(fif);
   istfi.t=zeros(sip(1),1);
   istfi.t(fimet)=fimet;
  end

  sip=size(nr.b);
  riga=[1:sip(1)]';
  righe=riga*ones(1,sip(2));
  fimet=find(imag(nr.b)<-.1 & righe>length(nr.b)-5);
  if length(fimet)>20
   istoxd=righe(fimet);
   istoxdu=max(istoxd);
%   istmet0=istoxdu;
   fif=find(istfi.b==istoxdu);
   if length(fif)==0
    istfi.b(istoxdu)=istoxdu;
   end
   fif=find(istfi.b~=0);
   istfi.b=istfi.b(fif);
  end
 end
istmet=0;
% disp('reass'), keyboard
 icsfib=0;
 icsfit=0;
 ifiplo=0;
 fi=find(istfi.t~=0);
 if length(fi)>0
  icsfit=length(fi);
  ifip=find(istfi.t<0);
  if length(ifip)>0
   ifiplo=-istfi.t(ifip);
  end
  istfi.t=abs(istfi.t(fi));
 else
  istfi.t=0;
 end

%'icsfit', keyboard

 fi=find(istfi.b~=0);
 if length(fi)>0
  icsfib=length(fi);
  ifip=find(istfi.b<0);
  if length(ifip)>0
   ifiplo=-istfi.b(ifip);
  end
  istfi.b=abs(istfi.b(fi));
 else
  istfi.b=0;
 end
% 'REASSIGN qui'
% keyboard

 if ifiplo==0
  ifiplo=find(ifield==-2);
  if length(ifiplo)>0
   icaplo=4;
  else
   icaplo=2;
  end
 end

 icsfi=icsfit+icsfib;
%'end reassign ', keyboard