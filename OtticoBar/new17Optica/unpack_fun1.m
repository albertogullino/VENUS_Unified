function [shu,Lu,nu,fiQ,fiC,fiQWi,ato,nlast,iate]=unpack_fun1(dv,nv,xm,radii,fst,ifp,shavet,iauto,anyf)


dv0=dv;
nv0=nv;
if ~exist('irepl_ind')
 irepl_ind=0;
end
%save fi
%save fi0
%' sopo dave', keyboard

fic=find(shavet==6);
av=radii.a;
av(fic)=0;

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
n2=nv(:,2:end);

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
fim=find(fl>0);
%fim=[fim; flu];
fw=length(fim);
fw=flu;

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
ntoe=[];
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
 if nl==0
  nl=1;
 end
 nlay=fix(abs(fst(pd,2)));
% [pd piv nl nlay], pausak
 if piv-pd~=0
  pup0=piv:pd-1;
  ddu=dv(pup0);
  imz=find(ddu>0);
  pup0=pup0(imz);
  ay=anyf(pup0);
  nto=[nto; n1(pup0)];
  ntoe=[ntoe; n2(pup0,:)];
  cto=[cto; c1(pup0)];
  fsto=[fsto; fst(pup0,2)];
  iauo=[iauo; iauto(pup0,:)];
  dto=[dto; dv(pup0)];
  ato=[ato; av(pup0,:)];
  shto=[shto; shavet(pup0,1)];  
  puo=[puo; pup0'.*ay];
 end
%'ICI 0 bir', keyboard
 pupe=pd:pd+nl-1;
 ay=anyf(pupe);
 for inp=1:nlay
  nto=[nto; n1(pupe)];
  ntoe=[ntoe; n2(pupe,:)];
  cto=[cto; c1(pupe)];
  fsto=[fsto; fst(pupe,2)];
  iauo=[iauo; iauto(pupe,:)];  
  dto=[dto; dv(pupe)];
  shto=[shto; shavet(pupe,1)];    
  if nlay==miup | nlay==midw
   if (inp==poup & nlay==miup) | (inp==podw & nlay==midw)
    ato=[ato; av(pupe,:)];
    puo=[puo; pupe'.*ay];
   else
    ato=[ato; av(pupe,:)];
    puo=[puo; pupe'*0];
   end
  else
   if inp==fix(nlay/2)
    ato=[ato; av(pupe,:)];
    puo=[puo; pupe'.*ay];
   else
    ato=[ato; av(pupe,:)];
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

%fiQ=leqw;
fiC=icav;
%fiQWi=find(fst(:,2)==-1);

leqw=find(fsto==-1);
%icav=find(iauo(:,2)==-4);
%iat=find(iauo(:,1)==2);

fiQ=iat;
%fiC=icav;
fiQWi=leqw;
%' fitqagf', keyboard

if piv<length(dv)+1
 pup0=piv:length(dv);
 nto=[nto; n1(pup0)];
 ntoe=[ntoe; n2(pup0,:)];
 
 cto=[cto; c1(pup0)];
 fsto=[fsto; fst(pup0,2)];
 dto=[dto; dv(pup0)];
 shto=[shto; shavet(pup0,1)];     
 ay=anyf(pup0);
 ato=[ato; av(pup0,:)];
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
figure, plot(ddplot,real(ndplot),'.-'), 
  title(' Struttura scompattata ')
pausak
end
%lmaxim=1e3;

shu=shto;
Lu=L_i;
nu=n_i;
ntot=[nu ntoe];

for k=1:length(ntot)
 nr=ntot(k,:);
 fim=find(abs(nr)>0);
 fiz=find(abs(nr)==0);
 nm=nr;
 nm(fiz)=nr(fim(end));
 nlast(k,:)=nm;
end


%'ver in unpack', keyboard