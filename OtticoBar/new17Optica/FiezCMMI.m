function [ga,la,zet,Ez,nz,Gaqw,NQW_ef,az]=fiezCMMI(lam0,L_i,n_i1,iate,dpac,npac,fpac,iat,fiQW,ifp,dlam,NPlam,rr,kt,ibast,ibaste,par_grat,a_i);


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
figure, plot(L_i,'o'), hold on, plot(iate,L_i(iate),'rp')
a=axis;
a(3)=0;
a(4)=max(L_i)+10;
axis(a)
pausak
end

lav=lam0+linspace(-1,1,NPlam)*dlam*lam0;
%lav=lam0+linspace(0,2,20)*dlam*lam0;
%lav=lam0-.01;






%lai=1.5928217
%keyboard
%[ei,g0]=Flam_gamLu(L_i,n_i,iat,icav,lai,rr,kt,ibast,par_grat);


for k=1:length(lav)
lai=lav(k);
%(L_in,n_i,iat,icav,la0,rr,kt,ibast,par_grat,g0);
%if k==1dpac,npac,fpac,
[ei,g0]=Flam_gamPack(dpac,npac,iat,fpac,lai,rr,kt,ibast,par_grat);
%'cont siy', keyboard
%else
%[ei,g0]=Flam_gamLu(L_i,n_i,iat,icav,lai,rr,kt,ibast,par_grat,real(ei));
%'passo'
%end
gov(k)=g0;
eiv(k)=ei;
end
alf0d=imag(eiv);
gai0d=real(eiv);

fipo=find(gai0d>0);

alf0=alf0d(fipo);
gai0=gai0d(fipo);
lav=lav(fipo);

if ifpc==-10
 figure
  subplot(211)
  plot(lav,alf0), grid
  subplot(212)
  semilogy(lav,gai0), grid, pausak
end


der=diff(alf0);

%ne=length(find(der<=0));
%po=length(find(der>0));
%if po==1 | ne==1
% if po==1
%  puf=(find(der>0)); 
% else
%  puf=(find(der<=0)); 
% end
% isw=1;
%else
% isw=0;
%end
%clear PU
% if isw==1
% PU{1}=1:puf;
% PU{2}=puf+2:length(lav);
%else 
% PU{1}=1:length(lav);
%end 
[du,imi]=min(abs(gai0));

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
 PU{1}=puv;

clear LV GV
kisol=0;
for kint=1:length(PU)
 pui=PU{kint};
 if length(pui)>2
 kisol=kisol+1;
 [du,imid] =min(abs(alf0(pui)));
 imi=pui(imid);
  if ifpc==-10
  figure,  
  subplot(211)
  plot(lav(pui),alf0(pui),lav(imi),0,'wo'), grid
  subplot(212)
  plot(lav(pui),gai0(pui),lav(imi),gai0(imi),'wo'), grid, pausak
  end 

 dlau=diff(lav(1:2));
 lav1=lav(imi)+linspace(-1,1,10)*dlau;
 clear gov eiv
 for k=1:length(lav1)
  lai=lav1(k);
  [ei,g0]=Flam_gamPack(dpac,npac,iat,fpac,lai,rr,kt,ibast,par_grat);  
  gov(k)=g0;
  eiv(k)=ei;
 end
 


%'qui CMMkk', keyboard

g00=g0;
alf=imag(eiv);
gv=real(eiv);
if length(lav1>7)
 [du,im]=min(abs(alf));
 pud=[-2:2]+im;
 puf=find(pud>0 & pud<=length(lav1));
 pu=pud(puf);
 lavf=lav1(pu);
 gvf=gv(pu);
 alvf=alf(pu);
 govf=gov(pu); 
else
 lavf=lav1;
 gvf=gv;
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
plot(lav1,eiv,LAM,GTH,'wo'), title('real part'), grid
subplot(212)
plot(lav1,imag(eiv),LAM,0,'wo'), title('imag part'), grid
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
%[ei]=Flam_gamLu(L_i,n_i,iat,icav,lai,rr,kt,ibast,par_grat,g0);
[ei]=Flam_gamPack(dpac,npac,iat,fpac,lai,rr,kt,ibast,par_grat,g0);
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

if ifpc==-10
' campi,' , pausak
end


Nx=50;
sNx=2;  %step in nm

%'qui nz', keyboard

[fiez,nz,zet,Gaqw,NQW_ef,az]=Flam_field1(L_i,n_i1,iate,fiQW,la,rr,kt,ibaste,par_grat,ga,sNx,a_i);



Ez=sum(fiez);
Ez=Ez/max(Ez);
Hz=diff(fiez);
Hz=Hz/max(Hz);

I=abs(Ez).^2;

figr=find(diff(zet)>.1);

if length(figr)==2
stm=2*mean(diff(zet));
ia=figr;
ib=ia+1;
fa=I(ia);
fap=diff(I(ia+[-2 0]))/diff(zet(ia+[-2 0]));
fb=I(ib);
fbp=diff(I(ib+[0 2]))/diff(zet(ib+[0 2]));
d=diff(zet([ia ib]));
za=zet(ia);
zb=zet(ib);
zed=[za:stm:zb];
zeh=linspace(za,zb,length(zed));
zeh=zeh(2:end-1);
df=fb-fa;
k=fbp/df;
fap=fa+df*exp(k*(zeh-zb));
if df*exp(k*(zeh(1)-zb))/fa>.1
 k=log(fb/fa)/(zb-za);
 fap=fa*exp(k*(zeh-za));
end
%figure, plot(zet,I,'r.',zeh,fap,'g.'),
' aggisto campo ', keyboard

In=[I(1:ia) fap*NaN I(ib:end)];
nz=[nz(1:ia,:); repmat(nz(ia,:),length(fap),1); nz(ib:end,:)];
zn=[zet(1:ia) zeh zet(ib:end)];

' cont', keyboard

Ez=sqrt(In);
zet=zn;
end


gatot=ga;
ga=gatot/NQW_ef;

%' cont NQE', keyboard

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
figure, plot(zet,I, zet, nz(:,1),'r'), 
 title(['CMM: ',lab,'  lambda_{res} = ',num2str(la),' Gth= ',num2str(ga)]), 
  xlabel(' Long. coord. (um)')
 ylabel(' Intens (Log10) / Index')
 a=axis;
 a(3)=I(end)-.5;
 a(3)=-2;
 axis(a)
'fine campi', pausak
end
%'fine campi', pausak
%keyboard
