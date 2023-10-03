clear all
close all
%load any
load abs
'qui Losses_zeta entro', keyboard

%fiQ=find(fst(:,2)==-1)-1;

%'fiQ ', keyboard


xpl.b.last=xm(pup.b.last,1);
xpl.b.o=xm(pup.b.o,1);
xpl.b.m=xm(pup.b.m,1);
xpl.b.i=xm(pup.b.i,1);

xpl.t.o=xm(pup.t.o,1);
xpl.t.m=xm(pup.t.m,1);
xpl.t.i=xm(pup.t.i,1);


for kstra=1:length(nv0)
% idu=find(xm(kstra,:)~=0);
 idu=1;
 [du,ip]=max(xm(kstra,idu));
 n1(kstra,1)=nv0(kstra,idu(ip));
 c1(kstra,1)=xm(kstra,idu(ip));
end
%'dopo'
%keyboard

%n1=(nv0(:,1));
%c1=xm(:,1);

fima=find(dv>1000);
if length(fima)>0
% 'attenzione riduco dv per campo numerico' , keyboard
% dv(fima)=1000;
end

xd=[zeros(size(dv)) dv];
fl=abs(fst(:,2));
fim=find(fl>1);
%fim=find(fl>=1);
if length(fim)>0
 if length(find(diff(fim)>1))==0 & sum(diff(fim))>length(fim)-1
  fim=[fim; length(dv)-1];
 end
else
 fim=length(dv);
end
fw=length(fim);
%fw=length(fl)-1;

%'fim', keyboard

miup=fst(fim(1),2);
midw=fst(fim(end),2);

% paio usato per EO
poup=miup-1;
podw=2;

% paio usato per EO
%poup=miup;
%podw=1;

dto=[];
nto=[];
ato=[];
cto=[];
puo=[];

%fista=find(abs(Dov)>0);
%piv=fista(1);
piv=1;
ico=1;
%'fim ', keyboard
fsto=[];
while fw>0
 pd=fim(ico);
 nl=fst(pd,1);
 nlay=abs(fst(pd,2));
 %[pd piv nl nlay], pausak
 if piv-pd~=0
  pup0=piv:pd-1;
  ddu=dv(pup0);
  imz=find(ddu>0);
  pup0=pup0(imz);
  if length(find(pup0==pua))==1
   fi0=find(pup0==pua);
   leqw=length(fsto)+fi0;
  end
  ay=anyf(pup0);
  nto=[nto; n1(pup0)];
  cto=[cto; c1(pup0)];
  fsto=[fsto; fst(pup0,2)];
  dto=[dto; dv(pup0)];
  ato=[ato; dv(pup0).*ay];
  puo=[puo; pup0'.*ay];
 end
%'ICI 0 bir', keyboard

 pupe=pd:pd+nl-1;
%'ICI 0 bir', keyboard
%if nlay>1
% [pd piv nl nlay fw], pausak
%end 

 ay=anyf(pupe);
 if length(pupe)==0
  nlay=0;
 end
 for inp=1:nlay
  nto=[nto; n1(pupe)];
  cto=[cto; c1(pupe)];
  fsto=[fsto; fst(pupe,2)];
  dto=[dto; dv(pupe)];
  %'if nlay', keyboard
  if nlay==miup | nlay==midw
   if (inp==poup & nlay==miup) | (inp==podw & nlay==midw)
    ato=[ato; dv(pupe).*ay];
    puo=[puo; pupe'.*ay];
   else
    ato=[ato; -dv(pupe).*ay];
    puo=[puo; pupe'*0];
   end
  else
   if inp==fix(nlay/2)
    ato=[ato; dv(pupe).*ay];
    puo=[puo; pupe'.*ay];
   else
    ato=[ato; -dv(pupe).*ay];
    puo=[puo; pupe'*0];
   end
  end
 end
 
% 'pupe', pausak
 if length(pupe)>0
 piv=pupe(length(pupe))+1;
 fw=fw-nl;
 ico=ico+nl;
 else
  fw=0;
  piv=length(dv);
 end
%  'pupe dopo', pausak
end

%'ICI 1 bir', keyboard

if piv<length(dv)+1
 pup0=piv:length(dv);
 nto=[nto; n1(pup0)];
 cto=[cto; c1(pup0)];
 fsto=[fsto; fst(pup0,2)];
 dto=[dto; dv(pup0)];
 ay=anyf(pup0);
 ato=[ato; dv(pup0).*ay];
 puo=[puo; pup0'.*ay];
end

puac=puo;
puae=find(ato>0);
%puae=find(ato>=0);
%' puae ', keyboard


nd=[nto nto];
ccd=[cto cto];
%pudto=
xd=[zeros(size(dto)) dto];
xt=cumsum(reshape(xd',2*length(dto),1));
nt=real(reshape(nd',2*length(dto),1));
ct=(reshape(ccd',2*length(dto),1));
L_i=dto;
n_i=nd(:,1);
nto(leqw)=real(nto(leqw));

lmaxim=1e6;
%lmaxim=1e3;
fitro=find(dto>lmaxim);
if length(fitro)>0
 dto(fitro)=lmaxim;
end
dJ=cumsum([0; dto]);
fJ=([0; fsto]);
nJ=([nto(1); nto]);
xJ=([cto(1); cto]);
lJ=length(dJ);

%' hz set', keyboard
sthz=1;
if length(find(dto>2e4)>0)
sthz=3;
end


%' primapo fipe', keyboard
%' primapo fipe', keyboard

fipe=find(real(perm)<0);
perm(fipe)=-perm(fipe);

relPerm=conj(perm.^2);
%relPerm=real(relPerm);

xdu=[dJ(1:end-1) dJ(2:end)];
ydu=[nJ(2:end) nJ(2:end)];

x=reshape(xdu',length(xdu)*2,1);
y=reshape(ydu',length(xdu)*2,1);

fipe=find(real(y)<0);
y(fipe)=-y(fipe);

%'pastnew ICI 3 bir', keyboard
%lambda0=974.5/1e9

%'qui EO', keyboard
 x=ztot;
 Fi=3*abs(real(Ez/max(Ez)).').^2;

 z=ABS.zAbs;
 %z=max(zv)-fliplr(zv);
 Dov=dglo.Dop;
 
 hz=x*1000;
 uFunc=hz*0;
 uF0=hz*0;
 perm=[];
 xfit=[];
 for il=1:lJ-1
  fiz=find(hz>=dJ(il) & hz<dJ(il+1));
  if length(fiz)>0
   perm=[perm nJ(il+1)*ones(size(fiz))];
   xfit=[xfit xJ(il+1)*ones(size(fiz))];
   if fJ(il+1)==-1
    uFunc([fiz(1)-1 fiz])=1;
   end
   if fJ(il+1)==-1 & il==leqw
    uF0([fiz(1)-1 fiz])=1;
   end
  end  %if 
 end
% ' primapo fipe', keyboard
 
 xfit(length(hz))=xfit(length(hz)-1);
 perm(length(hz))=perm(length(hz)-1);
 uFunc(length(hz))=0;
 uF0(length(hz))=0;
 
 
 %' primapo fipe', keyboard
 
 fipe=find(real(perm)<0);
 perm(fipe)=-perm(fipe);
 
relPerm=conj(perm.^2).';
 
 %fmu(:,1)=mean(EO.Es,2);


 fista=find(abs(Dov)>0);
 finiz=1:fista(1)-1;


 dorigin=sum(dv(finiz));
 dorigin1=sum(dv(finiz(2:end)));
 x=(hz'-dorigin)*1e-3;
 xv=x;

% 'dorig', keyboard
  zo=z-z(2);
  fip=find(nz>1);
  xp=x(fip);
  Fip=Fi(fip);
 
  fmu_0=flipud(ABS.eleccentro);
  fmu=interp1(zo,fmu_0,xp);
  fiN=find(isnan(fmu)==1);
  fmu(fiN)=0;
  fmuE=fmu*ABS.e;
  
  fmu_0=flipud(ABS.holecentro);
  fmu=interp1(zo,fmu_0,xp);
  fiN=find(isnan(fmu)==1);
  fmu(fiN)=0;
  fmuH=fmu*ABS.h;
  
  Los=fmuE+fmuH;
  
  Dsu=flipud(ABS.Dope+ABS.Doph);
  Do=interp1(zo,Dsu,xp);

  
  Xsu=flipud(ABS.xmol);
  xmol=interp1(zo,Xsu,xp);

  nref_mol=real(nAlGaAs(lambda*1e-6,xmol));
  %z=x;
%  ' primapo fipe', keyboard

%z=x';

% pewa=spline(x,imag(relPerm),z);
%'qui', keyboard
  figure, semilogy(xp,Los,'.',xp,Fip/max(Fip)*max(fmu),xp,real(nz(fip)/3.5),'g',xp,Do,'c',xp,xmol,'k.')
  ylim([.001 50])
  pausak
  figure, plot(xp,real(nz(fip)),'g',xp,nref_mol,'k.')
  pausak

 FPerd=Fip.*Los;



 dz=[diff(xp); 0];
 smu=size(Los);
 F_dz=((Fip.*dz)*ones(1,smu(2)));
 FL_dz=((FPerd.*dz)*ones(1,smu(2)));
% Fnorm=Fiw.*dz;
% Fnorms=Fiw'*dz;
 dmic=(dto)/1000;
 fid=find(dmic==0);
 if length(fid)>0
  dmic(fic)=1.2e-7;
 end

 fista=find(abs(Dov)>0);
 kin=fista(1);
clear amedis 
 for k=1:length(puae)
  li=sum(dmic(1:puae(k)-1))-dmic(1);
  lu=li+dmic(puae(k));
  fi=find(xp>=li & xp<lu);
  R=sum(FL_dz(fi))/sum(F_dz(fi))
  amedis(k)=R;
%  [li lu],  k,   pausak
 end


 puc=puac(find(puac>0));

 Lo=zeros(size(dv));
 Lo(puc)=amedis; 
 k0cm2=k0cm*2;
 if ifp==-10
  figure, semilogy(cumsum(dv),-imag(nv(:,1))*k0cm2+.1,'o')
  hold on
  semilogy(cumsum(dv),Lo,'r.')
  pausak
 end 
 
 xmod=xm(:,1);
 fim=find(xmod>.01);
 nref_new=real(nAlGaAs(lambda*1e-6,xmod(fim)));
 nrdu=zeros(size(nv));
 nrdu(fim,1)=nref_new-1i*Lo(fim)/k0cm2;
 
 nr_sav=nr;
 

 nr.t(:,1)=nrdu(put,1);
 nr.b(:,1)=nrdu(pub,1);
 nr.a(:,1)=nrdu(pua,1);
 if ifp==-10
 'aggiornati indici e perdite', keyboard
 end