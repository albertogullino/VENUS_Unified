'qui EO entro', keyboard

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

 global EO
 z=EO.zT';
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

 
  axorien=1;  % 1: x=[110],  -1: x=[1-10],
  r411=-1.6e-6;  % in micron/V
  r412=-1.1e-6;
  core=polyfit([0 1],[r411 r412],1);
  r41vet=polyval(core,xfit);
  n2r41=r41vet.*real(perm).^2;
  epsr41=axorien*2*r41vet.*real(perm).^4;
  fi0=find(xfit<0);
  n2r41(fi0)=0;
  epsr41(fi0)=0;

 fista=find(abs(Dov)>0);
 finiz=1:fista(1)-1;


 dorigin=sum(dv(finiz));
 dorigin1=sum(dv(finiz(2:end)));
 x=(hz'-dorigin)*1e-3;
 Fiany=Fi.*n2r41';
 Fianyeps=Fi.*epsr41';
 %'dorig', keyboard
 
  fmu_0=mean(EO.Es,2);
  fmu=interp1(z,fmu_0,x);
  fiN=find(isnan(fmu)==1);
  fmu(fiN)=0;
  zo=z;
  z=x;
%  ' primapo fipe', keyboard

%z=x';
 Fiwaeps=spline(x,Fianyeps,z);
 Fiwa=spline(x,Fiany,z);
 Fiw=spline(x,Fi,z);
% pewa=spline(x,imag(relPerm),z);
%'qui', keyboard

clear amedis ame amz amzd ameps

 dz=[diff(z); 0];
 smu=size(fmu);
 F_dz=((Fiwa.*dz)*ones(1,smu(2))).*fmu;
 F_dzeps=((Fiwaeps.*dz)*ones(1,smu(2))).*fmu;
 F_dis=((Fiw.*dz)*ones(1,smu(2))).*fmu;
 Fnorm=Fiw.*dz;
 Fnorms=Fiw'*dz;
 dmic=(dto)/1000;
 fid=find(dmic==0);
 if length(fid)>0
  dmic(fic)=1.2e-7;
 end

 fista=find(abs(Dov)>0);
 kin=fista(1);

 for k=1:length(puae)
%  li=sum(dmic(1:puae(k)-1))-dto(1)/1000;
  li=sum(dmic(1:puae(k)-1))-dorigin/1000;
  lu=li+dmic(puae(k));
  fi=find(z>=li & z<lu);
  amedis(k,:)=sum(F_dis(fi,:))/sum(Fnorm(fi));
  du=sum(F_dzeps(fi,:))/sum(Fnorm(fi));
  if isnan(du)
  du=0;
  ' du = 0 ', keyboard
  end
  ameps(k,:)=du;
  ame(k,:)=sum(F_dz(fi,:))/sum(Fnorm(fi));
  amz(k,:)=li;
  amzd(k,:)=dmic(puae(k));
 end

 iaa=1;
 amedis=amedis(:,iaa)';
 ad=[amedis amedis];
 xd=[amz' (amz+amzd)'];
 xat=(reshape(xd',2*length(amz),1));
 %xd=[zeros(size(amz')) amz'];
 %xat=cumsum(reshape(xd',2*length(amz),1));
 at=(reshape(ad',2*length(amz),1));

 xts=(xt-dorigin)/1000;
if ifp>=-1 | ifp==-10

 figure, plot(z,fmu,x,real(sqrt(relPerm)),'r',z,Fiw*5,'w',...
          xts,nt,'c',xts,nt,'m.',xat,at,'r.'),
          pausak
 figure, plot(z,fmu(:,iaa),'r',z,Fiw*50,'b',xts,nt.^2,'g'),
 pausak
% figure, plot(z,fmu/100,x,sqrt(real(relPerm)),'y',z,Fiw*5,'w')
%          xts,nt*10,'c',xts,nt,'m.',xat,at,'r.'),
% if ifp>1, pausak, end
end


 anyv=zeros(length(dv),smu(2));

 puc=puac(find(puac>0));

 anyv(puc,:)=ameps;

% figure, plot(ameps), keyboard
 puaf=[kin:length(dto)];
 for k=1:length(puaf)
%  li=sum(dmic(1:puaf(k)-1))-dto(1)/1000;
  li=sum(dmic(1:puaf(k)-1))-dorigin/1000;
  lu=li+dmic(puaf(k));
  fi=find(z>=li & z<lu);
  if length(fi)==0
   Amedis(k)=0;
   Ame(k)=0;
  else
   sFn=sum(Fnorm(fi));
   if sFn~=0
    Amedis(k)=sum(F_dis(fi))/sFn;
    Ame(k)=sum(F_dz(fi))/sFn;
%    plot(z(fi),F_dis(fi),'.-',z(fi),Amedis(k)*ones(size(fi)),'w'), pausak
   else
    Amedis(k)=0;
    Ame(k)=0;
   end
  end
   Amz(k)=li;
   Amzd(k)=dmic(puaf(k));
 end
 fizd=find(Amzd==0);
 Amzd(fizd)=1e-7;

 Ad=[Amedis' Amedis'];
 Ade=[Ame' Ame'];
 Xd=[Amz' (Amz+Amzd)'-1e-6];
 Xat=(reshape(Xd',2*length(Amz),1));
 At=(reshape(Ad',2*length(Amz),1));
 Ate=(reshape(Ade',2*length(Amz),1));

 Atf=interp1(Xat,At,z,'linear');

 fiz=find(z<Xat(1));
 z(fiz)=Xat(1);
 Atfe=interp1(Xat,Ate,z,'linear');


if ifp>=-1 | ifp==-10
%if ifp>=-1
 figure, plot(z,fmu,x,3*real(sqrt(relPerm)),'c',z,Fiw*5,'w',...
         xts,3*nt,'b.',xat,at,'wo',Xat,At,'r.',z,Atf),
pausak
% figure, plot(z,fmu,xat,at,'wo',Xat,At,'r.',z,Atf),


  fip=find(fmu>0);
  fim=find(fmu<0);
  Fp=fmu(fip)/3;
  Fm=-fmu(fim)/3;
  figure, plot(z(fip),Fp,'y',z(fim),Fm,'r',z,Fiw*3,'w',xts,nt*2,'c')
  pausak
 if ifp>1, pausak, end
end


% return

 F_dz=Fiwa.*dz;
 ne=min(size(fmu));
 Fm=ones(ne,1)*F_dz';
 cSup=(Fm'.*fmu)/Fnorms;
 Sup=(F_dz'*fmu)/Fnorms;
 cobir=3e5/lambda;
 birt=Sup*cobir;

% fi0=find(pewa~=0);
 zcav_inf=sum(dv(1:fiC(1)-1).*abs(fst(1:fiC(1)-1,2)))-dorigin;
 zcav_sup=sum(dv(1:fiC(end)).*abs(fst(1:fiC(end),2)))-dorigin;
 z1000=z*1000;
 fi0=find(z1000<zcav_inf | z1000>zcav_sup);
 F_dz0=Fiwa.*dz;
 F_dz0(fi0)=0;
 Sup0=(F_dz0'*fmu)/Fnorms;
 bircav=Sup0*cobir;

% fi0=find(pewa==0);
 fi0=find(z1000>=zcav_inf & z1000<=zcav_sup);
 F_dz0=Fiwa.*dz;
 %F_dz0(fi0)=0;
 Sup0=(F_dz0'*fmu)/Fnorms;
 birmir=Sup0*cobir;


% verifica

 F_dz=Atfe.*Fiw.*dz;
 ne=min(size(Atfe));
 Fm=ones(ne,1)*F_dz';
 inonan=find(isnan(F_dz)==0);
 Supv=sum(F_dz(inonan))/Fnorms;
 cSupv=(F_dz)/Fnorms;
 cobir=3e5/lambda;
 birt_verifica=Supv*cobir;
 if ifp==-10
  disp(' ');
  disp(' verifica Bir. ');
  [birt       birt_verifica]
  disp(' ');
  keyboard
 end

 if ifp>=-1 | ifp==-10
  figure, plot(z,cumsum(cSup)*cobir,z,cumsum(cSupv)*cobir,'-'), title(' fy-fx')
  pausak
 if ifp>1, pausak, end
 end

 any.t=anyv(put,:);
 any.b=anyv(pub,:);
 any.a=anyv(pua,:);
 if ifp==-10
 'any', keyboard
 end