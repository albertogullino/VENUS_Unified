function [ga,la,zet,Ez,nz,Gaqw,NQW_ef,az]=fiezCMMn(lam0,L_i,n_i1,iat,icav,fiQW,ifp,dlam,NPlam,rr,kt,ibast,par_grat,a_i);

%' in fiez vecchio ', keyboard

n_i=n_i1(:,1);
ifpc=-4;  %controllo soluzione
iplotCa=0;
%ifp=-10
if ifp<-10
ifpc=-10;  %controllo soluzione
end
if ifp<=-10
iplotCa=1;
end

if ~exist('rr')
 rr=3;
end 
if ~exist('kt')
 kt=0;
end 
if ~exist('ibast')
 ibast=[];
 par_grat=0;
else
 if ibast>0
  if ~exist('par_grat')
   'missing grating parameters ', keyboard
  end
 end
end 

PO=[909    83   672   711];
if ifp==-11
figure, plot(L_i,'o'), hold on, plot(iat,L_i(iat),'rp')
a=axis;
a(3)=0;
a(4)=max(L_i)+10;
axis(a)
pausak
end

lav=lam0+linspace(-1,1,NPlam)*dlam*lam0;
%lav=lam0+linspace(0,2,20)*dlam*lam0;
%lav=lam0;


for k=1:length(lav)
lai=lav(k);

[ei,g0]=Flam_gamLu(L_i,n_i,iat,icav,lai,rr,kt,ibast,par_grat);
%'dopo Flam vecchio', keyboard
gov(k)=g0;
eiv(k)=ei;
end
alf0d=imag(eiv);
gai0d=real(eiv);

fipo=find(gai0d>0);

alf0=alf0d(fipo);
gai0=gai0d(fipo);
lav=lav(fipo);

az=alf0(1:end-1).*alf0(2:end);
fiz=find(az<0);
if length(fiz)>0
 if fiz(end)==length(lav)
  fiz(end)=length(lav)-1;
 end
 gafiz=(gai0(fiz)+gai0(fiz+1))/2;
 [du,imi]=min(gafiz);
 fiz=fiz(imi);
end
der=diff(alf0);


ive=0;
if ive==1
 if fiz(end)==length(lav)
  fiz(end)=length(lav)-1;
 end
 gafiz=(gai0(fiz)+gai0(fiz+1))/2;
[du,imi]=min(gafiz);
der=alv0;
ip=1;
itro=0;
for ki=1:length(der)
 puv(ip)=ki;
 ip=ip+1;
  if ki>imi
   itro=1;
  end
  if der(ki)>0 
   if itro==0
    clear puv
    ip=1;
   else
    break
   end
  end
end 
%'trovati punti', keyboard
else
 if fiz>1
  puv=fiz+[-1 0 1 2];
 else
   puv=fiz+[0 1 2];
 end
end

if ifpc==-10
 figure
  subplot(211)
  plot(lav,alf0,lav(puv),alf0(puv),'r.'), grid
  subplot(212)
  semilogy(lav,gai0,lav(puv),gai0(puv),'r.'), grid, pausak
end

%'trovati punti', keyboard
 PU{1}=puv;

clear LV GV
kisol=0;
for kint=1:length(PU)
 pui=PU{kint};
 if length(pui)>2
 kisol=kisol+1;
 [du,imid] =min(abs(alf0(pui)));
 imi=pui(imid);
lav1=lav;
%imi=1;
fiz=[imi imi];
Npz=21;

da=diff(alf0);
dg=diff(diff(gai0));
dg=[dg(1) dg];
fipa=length(find(da>0));
fima=length(find(da<0));
fipg=length(find(dg>0));
fimg=length(find(dg<0));
if fipa*fima==0 & fipg*fimg==0
 npun=Npz;
else
 npun=50;
end
imi=1;
ic=0
lav1=linspace(lav1(puv(1)),lav1(puv(end)),npun);
%' primo lav1', keyboard
dii=2;
Npz=21;
dlau=diff(lav1(1:2));
ibre=0;
while ic<=2
 if ibre==1
  break
 end
 if npun==Npz
  ibre=1;
 end 
 ic=ic+1;
 if ic>1
  lav1=lac+linspace(-dii/2,dii*3/2,npun)*dlau;
 end 
% ' qui', keyboard
 clear gov eiv
 for k=1:length(lav1)
  lai=lav1(k);
  [ei,g0]=Flam_gamLu(L_i,n_i,iat,icav,lai,rr,kt,ibast,par_grat);

  gov(k)=g0;
%  if imag(ei)<0
%   'g0 < 0' , pausak
%  end
  eiv(k)=ei;
 end
 


%'qui CMMkk', keyboard

g00=g0;
alf=imag(eiv);
gav=real(eiv);

if npun==50
ad=alf(1:end-1).*alf(2:end);
fiz0=find(ad<0);
[du,fid]=min(abs(gav(fiz0)));
fiz=fiz0(fid);
else
if ic<2
[du,fiz]=min(abs(alf)+10*abs(gav));
else
[du,fiz]=min(abs(alf));
end
end


  npun=Npz;
  dii=1;
  dlau=diff(lav1(1:2));
  lac=lav1(fiz(imi));
%' lac', keyboard

  if ifpc==-10
  figure,  
  subplot(211)
  plot(lav1,alf,'.-',lav1(fiz),alf(fiz),'ro'), grid
    subplot(212)
  semilogy(lav1,gav,lav1(fiz),gav(fiz),'ro'), grid
  pausak
  end  
  
 %[du,imi]=min(gav(fiz)); 
 imi=1;
 
 
 %' fine while ', pausak
 end
 
 puz=fiz+[-1 0 1];
 lav1=lav1(puz);
 gav=gav(puz);
 alf=alf(puz);
 gov=gov(puz);

if length(lav1>7)
 [du,im]=min(abs(alf));
 pud=[-2:2]+im;
 puf=find(pud>0 & pud<=length(lav1));
 pu=pud(puf);
 lavf=lav1(pu);
 gvf=gav(pu);
 alvf=alf(pu);
 govf=gov(pu); 
else
 lavf=lav1;
 gvf=gav;
 alvf=alf;
 govf=gov;
end

coa=polyfit(lavf,alvf,1);
LAM=roots(coa);

if ifpc==-10
figure, plot(lavf,alvf,LAM,0,'wo'), pausak
end
cog=polyfit(lavf,gvf,2);
GTH=polyval(cog,LAM);
G0=GTH;
cog0=polyfit(lavf,govf,2);
Gga=polyval(cog0,LAM);

if ifpc==-10
h=figure, 
set(h,'pos',PO)
subplot(211)
plot(lav1,gav,LAM,GTH,'wo'), title('real part'), grid
subplot(212)
plot(lav1,alf,LAM,0,'wo'), title('imag part'), grid
pausak
end
 GV(kisol)=G0;
 LV(kisol)=LAM;
 end %length(pui)
end   %fine kint


if ifpc==-10
GV
[G0,imi]=min(GV);
G0
LAM=LV(imi)
'inizio iterazione', pausak
end

if ifpc==-10
 figure
  subplot(211)
  plot(lav,alf0,LAM,0,'ro'), grid
  subplot(212)
  semilogy(lav,gai0,LAM,G0,'ro'), grid, pausak
end

G00=G0;
Gprec=G0;
G0=0;
iter=1;
Gve(iter)=Gprec;
Lve(iter)=LAM;
while abs(G0/Gprec-1)>1e-3 & iter<2

lav=LAM*(1+linspace(-2,1,2)*10^-(3+iter));
lav=LAM;
g0=Gprec*1.1;
%g0=Gprec;
g0=Gprec*.9;

for k=1:length(lav)
lai=lav(k);
%[ei,g0]=Flam_gamLu(L_i,n_i,iat,icav,lai,rr,kt);
[ei]=Flam_gamLu(L_i,n_i,iat,icav,lai,rr,kt,ibast,par_grat,g0);
%[ei]=FlamU(L_i,n_i,iat,lai,g0);
eiv1(k)=ei;
end
eiv=eiv1;

G0=real(eiv);
iter=iter+1;
Gp=Gprec;
Gprec=G0;
Gve(iter)=Gprec;
Lve(iter)=LAM;

G0=Gp;
%Gprec=G0;
abs(G0/Gprec-1);
%pausak
end


Gx=[Gga g0];

Gx=g0*[1 1.2];
cg=polyfit(Gx,Gve,1);
x=linspace(min(Gx),max(Gx),10);
cg0=cg;
cg0(1)=cg0(1)-1;
gg0=roots(cg0);
Gstim=gg0;
ga=Gstim;
la=LAM;

if ifpc==-10
figure, plot(x,polyval(cg,x),gg0,gg0,'ro')
pausak
end

%' campi,' , pausak

%'qui nz', keyboard
Nx=50;
%Nx=2;

[fiez,nz,zet,Gaqw,NQW_ef,az]=Flam_field(L_i,n_i1,iat,fiQW,la,rr,kt,ibast,par_grat,ga,Nx,a_i);

Ez=sum(fiez);
Ez=Ez/max(Ez);
Hz=diff(fiez);
Hz=Hz/max(Hz);

fize=find(nz==0);
nz(fize,1)=NaN;

gatot=ga;
ga=gatot/NQW_ef;

if par_grat.itetm==1
 lab='TE';
else
 lab='TM';
end
if iplotCa==1
zqw=sum(L_i(1:fiQW(1)))/1000;
[du,imi]=min(abs(zet-zqw));
I0=abs(Ez).^2;
I=log10(I0/I0(imi))+3;
I(1,fize)=NaN;
h=figure, 
set(h,'pos',[150         300        1200         350])
subplot(121)
plot(zet,I, zet, nz(:,1),'r'), 
 title(['CMM: ',lab,'  lambda_{res} = ',num2str(la),' Gth= ',num2str(ga)]), 
  xlabel(' Long. coord. (um)')
 ylabel(' Intens (Log10) / Index')
 a=axis;
 a(3)=-1.5;
 a(4)=5;
 axis(a)
 subplot(122)
  plot(zet,I0*max(real(n_i)), zet, nz(:,1),'r'), 
  title(['CMM: ',lab,'  lambda_{res} = ',num2str(la),' Gth= ',num2str(ga)]), 
   xlabel(' Long. coord. (um)')
  ylabel(' Intens  / Index')
'fine campi', pausak
end
%'fine campi', pausak
%keyboard
