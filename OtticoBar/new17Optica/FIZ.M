clear
close all
load  fi
%load  fi0
for kstra=1:length(nv0)
% idu=find(xm(kstra,:)~=0);
 idu=1;
 [du,ip]=max(xm(kstra,idu));
 n1(kstra,1)=nv0(kstra,idu(ip));
 c1(kstra,1)=xm(kstra,idu(ip));
end

xd=[zeros(size(dv)) dv];
fl=abs(fst(:,2));
fim=find(fl>1);
fw=length(fim);


miup=fst(fim(1),2);
midw=fst(fim(end),2);
poup=miup-4;
podw=4;

dto=[];
nto=[];
ato=[];
cto=[];
puo=[];

piv=1;
ico=1;

fst=abs(fst);
fi=find(fst(:,1)==0);
fst(fi,1)=1;
figure, plot(fst(:,1),'.'), hold on, plot(fst(:,2),'ro'), pausak

%%%%

fstu=fst;
nilu=nv(:,1);
dilu=dv;
dsu=cumsum(dilu);
sfst=fstu(:,2);
itro=find(diff(sfst)~=0);
itro=[itro; length(fstu)];
ki=1;
xd=[];
yd=[];
nd=[];
sup=0;
for k=1:length(itro)

 ku=itro(k);
%' [k ki ku ] '
% [k ki ku ]
% pausak
 pu=ki:ku;
 ydu=dsu(pu);
 ndu=real(nilu(pu,1));
 RP=fstu(ki,2);
 for kk=1:RP
  if kk>1
   sup=sup+per0;
  else
   if RP>1
    per0=sum(dilu(pu));
    per0
 %   ' calcolo per0', keyboard
   end
  end
  yd=[yd; ydu+sup];
  nd=[nd; ndu];
 end
 ki=ku+1;
end
nds=nd;
figure, plot(yd,nd), pausak
% load sac, hold on, plot(yds,nds,'r'), grid
yds=yd;
nds=nd;
save sac yds nds
'salveto'
pausak

%%%


fsto=[];
while fw>0
 pd=fim(ico);
 nl=fst(pd,1);
 nlay=abs(fst(pd,2));
%` [pd piv nl nlay], pausak
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
%'ICI 0 bir', pausak
 pupe=pd:pd+nl-1;
 ay=anyf(pupe);
 for inp=1:nlay
  nto=[nto; n1(pupe)];
  cto=[cto; c1(pupe)];
  fsto=[fsto; fst(pupe,2)];
  dto=[dto; dv(pupe)];
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
 piv=pupe(length(pupe))+1;
 fw=fw-nl;
 ico=ico+nl;
end
 figure, plot(cumsum(real(dto)),real(nto),'.')
 cs=cumsum(real(dto)); rn=real(nto);
 save co cs rn
% load co, hold on, plot(cs,rn,'ro')
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

lmaxim=1e6;
%lmaxim=1e3;
fitro=find(dto>lmaxim);
if length(fitro)>0
 dto(fitro)=lmaxim;
end
dJ=cumsum([0; dto]);
fJ=([0; fsto]);
nJ=([nto(1); nto]);


dJ=[0; yd];
nJ=[nds(1); nds];

xJ=([cto(1); cto]);
lJ=length(dJ);
hz=[dJ(1):dJ(lJ)];
uFunc=hz*0;
uF0=hz*0;
perm=[];
xfit=[];
for il=1:lJ-1
 fiz=find(hz>=dJ(il) & hz<dJ(il+1));
 perm=[perm nJ(il+1)*ones(size(fiz))];
 xfit=[xfit xJ(il+1)*ones(size(fiz))];
 if fJ(il+1)==-1
  uFunc([fiz(1)-1 fiz])=1;
 end
 if fJ(il+1)==-1 & il==leqw
  uF0([fiz(1)-1 fiz])=1;
 end
end
xfit(length(hz))=xfit(length(hz)-1);
perm(length(hz))=perm(length(hz)-1);
uFunc(length(hz))=0;
uF0(length(hz))=0;
relPerm=conj(perm.^2);
%relPerm=real(relPerm);

if ifp==-10
rp=real(relPerm);
rpc=rp;
figure, plot(hz/1000,rp)
hzc=hz/1000;
pausak
save sa hzc  rpc
% load sa, hold on, plot(hzc,rpc,'r'), grid
'ICI 3 bir', keyboard
end
% lambda0=lambda/1e6

 [Ksi,lambdas,Fi,uLong,uLong0,Fa]=eiglmio(lambda0,uFunc,uF0,relPerm);
 uL(1)=uLong0;
 uL(2)=uLong;
 fatqw=uLong/(uLong0*nmir.a);
 lambda=lambdas*1e6;
 confzv=uL;


 dorigin=0;
 x=(hz'-dorigin)*1e-3;
 rperm=real(sqrt(relPerm));
 if ifp==-10
  figure, plot(x,rperm,'r',x,Fi*4,'w')
  title([' lambda_{res} = ',num2str(lambda)])
  pausak
 end

