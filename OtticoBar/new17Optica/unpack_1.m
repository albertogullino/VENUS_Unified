if ~exist('irepl_ind')
 irepl_ind=0;
end
%save fi
%save fi0
%' sopo dave', keyboard
for kstra=1:length(nv0)
% idu=find(xm(kstra,:)~=0);
 idu=1;
 [du,ip]=max(xm(kstra,idu));
% n1(kstra,1)=nv0(kstra,idu(ip));
 c1(kstra,1)=xm(kstra,idu(ip));
end

radi=radii.a;
finst=find(fst(:,2)==-1);
radi(finst,:)=0;


n1=fette_sub(0,dv,radi,nv,ifp);
%' dopo fette', keyboard
%ifp=-10

if irepl_ind==1
%' irepl=1', keyboard
 fi_re=find(shavet==4);
 if length(fi_re)>0
  nrepla=nv0(fi_re(1)-1,1);
  n1(fi_re)=nrepla;
 end
end 

xd=[zeros(size(dv)) dv];
fl=abs(fst(:,2));
flu=find(iauto(:,1)==3);
fim=find(fl>1);
%fim=[fim; flu];
fw=length(fim);

fime=find(real(n1)<0);
n1(fime)=-n1(fime);

if length(fim)>0
miup=fst(fim(1),2);
midw=fst(fim(end),2);
else
miup=0;
midw=0;
end
poup=miup-4;
podw=4;

dto=[];
nto=[];
ato=[];
shto=[];
cto=[];
puo=[];

%' qui nofield ', keyboard
%fif=find(fst(:,1)==0);
%fst(fif,1)=1;

piv=1;
ico=1;
fsto=[];
iauo=[];
while fw>0 
%while fw>0 | leqw==0
 pd=fim(ico);
 nl=fst(pd,1);
 nlay=fix(abs(fst(pd,2)));
% [pd piv nl nlay], pausak
 if piv-pd~=0
  pup0=piv:pd-1;
  ddu=dv(pup0);
  imz=find(ddu>0);
  pup0=pup0(imz);
  ay=anyf(pup0);
  nto=[nto; n1(pup0)];
  cto=[cto; c1(pup0)];
  fsto=[fsto; fst(pup0,2)];
  iauo=[iauo; iauto(pup0,:)];
  dto=[dto; dv(pup0)];
  ato=[ato; dv(pup0).*ay];
  shto=[shto; shavet(pup0,1)];  
  puo=[puo; pup0'.*ay];
 end
%'ICI 0 bir', keyboard
 pupe=pd:pd+nl-1;
 ay=anyf(pupe);
 for inp=1:nlay
  nto=[nto; n1(pupe)];
  cto=[cto; c1(pupe)];
  fsto=[fsto; fst(pupe,2)];
  iauo=[iauo; iauto(pupe,:)];  
  dto=[dto; dv(pupe)];
  shto=[shto; shavet(pupe,1)];    
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
% if length(find(dto>50000))>0
%  'dto', keyboard
% end
end

%'ICI 1 bir', keyboard
leqw=find(fst(:,1)==-1);
icav=find(iauto(:,2)==-4);
iat=find(iauto(:,1)==2);
iate=find(iauo(:,1)==2);

fiQ=leqw;
fiC=icav;
fiQWi=find(fst(:,2)==-1);

%' pack', keyboard

if piv<length(dv)+1
 pup0=piv:length(dv);
 nto=[nto; n1(pup0)];
 cto=[cto; c1(pup0)];
 fsto=[fsto; fst(pup0,2)];
 dto=[dto; dv(pup0)];
 shto=[shto; shavet(pup0,1)];     
 ay=anyf(pup0);
 ato=[ato; dv(pup0).*ay];
 puo=[puo; pup0'.*ay];
end

puac=puo;
puae=find(ato>0);
leqw=find(iauo==2);
%puae=find(ato>=0);


nd=[nto nto];
ccd=[cto cto];
%pudto=
xd=[zeros(size(dto)) dto];
xt=cumsum(reshape(xd',2*length(dto),1));
nt=real(reshape(nd',2*length(dto),1));
ct=(reshape(ccd',2*length(dto),1));
L_i=dto;
n_i=nd(:,1);
%'fi_z', keyboard
lmaxim=1e6;

nd=[nto nto];
ndplot=reshape(nd.',prod(size(nd)),1);
ddr=dto;
dd=[zeros(size(ddr))  ddr];
ddplot=cumsum(reshape(dd.',prod(size(nd)),1));


if ifp<=-10
figure, plot(ddplot,real(ndplot),'.-'), pausak
end
%lmaxim=1e3;

return
fitro=find(dto>lmaxim);
if length(fitro)>0
 dto(fitro)=lmaxim;
end

global Lmax
if length(Lmax)==0
 Lmax=150000;
end
fidt=find(dto>Lmax);
%'qui dto', keyboard
dto(fidt)=Lmax;
dJ=cumsum([0; dto]);
fJ=([0; fsto]);
nJ=([nto(1); nto]);
xJ=([cto(1); cto]);
lJ=length(dJ);
hn=STZ0;
hz=[dJ(1):hn:dJ(lJ)];
uFunc=hz*0;
uF0=hz*0;
perm=[];
xfit=[];
dJ=floor(dJ*1e4)*1e-4;
for il=1:lJ-1
 fiz=find(hz>dJ(il) & hz<=dJ(il+1));

% fiz=find(hz>=dJ(il) & hz<dJ(il+1));
 perm=[perm nJ(il+1)*ones(size(fiz))];
 xfit=[xfit xJ(il+1)*ones(size(fiz))];
 %il, 
 %fJ(il+1)
 %dJ(il)
 %pausak
 if fJ(il+1)==-1
%  uFunc([fiz(1)-1 fiz])=1;
  uFunc([fiz])=1;
 end
 if fJ(il+1)==-1 & il==leqw
  uF0([fiz])=1;
%  uF0([fiz(1)-1 fiz])=1;
% length(fiz)
 end
end
xfit(length(hz))=xfit(length(hz)-1);
perm(length(hz))=perm(length(hz)-1);
uFunc(length(hz))=0;
uF0(length(hz))=0;
relPerm=conj(perm.^2);
%relPerm=real(relPerm);
nd=[nto nto];
ndplot=reshape(nd.',prod(size(nd)),1);
ddr=dto;
dd=[zeros(size(ddr))  ddr];
ddplot=cumsum(reshape(dd.',prod(size(nd)),1));


if ifp<=-10
rp=real(sqrt(relPerm));
rpc=rp;
figure, plot(hz,rp,'r')
hzc=hz/1000;
hold on
plot(ddplot,real(ndplot),'.-')

%load sa,  hold on, plot(hzc,rpc,'r')
%keyboard
pausak
end
%'qui prima di mio dopo', keyboard
%'qui prima di mio dopo', keyboard
if ~exist('lambda0')
 lambda0=lambda/1e6
end

%'qui prima di mio dopo', keyboard
return
%lambda0=866*1e-9

% dorigin=0;
% x=(hz'-dorigin)*1e-3;
% rperm=real(sqrt(relPerm));
% load pri
% rpermp=real(sqrt(relPermp));
%  figure, plot(rperm,'r'), hold on, plot(rpermp,'g')
% figure, plot(
% lambda0p=lambda0; uFuncp=uFunc; uF0p=uF0; relPermp=relPerm;
% save pri lambda0p uFuncp uF0p relPermp
if nmir.a==1
 uF0=uFunc;
end
 [Ksi,lambdas,Fi,uLong,uLong0,Fa]=eiglmiu(lambda0,uFunc,uF0,relPerm,hn);
 uL(1)=uLong0;
 uL(2)=uLong;
 fatqw=uLong/(uLong0*nmir.a);
 lambda=lambdas*1e6;
 gpla=2e4*pi*rr/lambda*imag(Ksi)/uLong0/fatqw; 
 gpla=2e4*pi*rr/lambda*imag(Ksi)/uLong; 
 manv=max(uF0.*Fi.');
 confzv=uL;
%'ulong ', keyboard
 dorigin=0;
 x=(hz'-dorigin)*1e-3;
 rperm=real(sqrt(relPerm));
 faca=3;
 if ifp<=-10
  figure, plot(x,rperm,'r',x,faca*Fi/manv,'w')
%    title([' lambda_{res} = ',num2str(lambda),'  Ksi = ',num2str(Ksi)]), pausak
    title([' lambda_{res} = ',num2str(lambda),'  Gth = ',num2str(gpla)]), pausak
  pausak
 man=max(uF0.*Fi.')/3.5; 
 figure, plot(hz,real(perm),'w',hz,faca*Fi/man,'r',hz,-imag(perm)*6000,'c'),
 ax=axis;
 ax(4)=5;
 axis(ax)
  title([' lambda_{res} = ',num2str(lambda)]), pausak  
 end
