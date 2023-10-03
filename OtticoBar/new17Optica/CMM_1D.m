function [ga,la,zet,nz,Ez,Hz,uL]=CMM_1D(lam0,L_i,n_i,iat,icav,fiQW,ifp,dlam,rr,kt,ibast,par_grat,i_campi,a_i);

i_campi=1;
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

lav=lam0+linspace(-1,1,20)*dlam*lam0;
%lav=lam0+linspace(0,2,20)*dlam*lam0;
%lav=lam0-.01;


for k=1:length(lav)
lai=lav(k);
%(L_in,n_i,iat,icav,la0,rr,kt,ibast,par_grat,g0);
%if k==1
[ei,g0]=Flam_gamLu(L_i,n_i,iat,icav,lai,rr,kt,ibast,par_grat);
%'cont 1', keyboard
%else
%[ei,g0]=Flam_gamLu(L_i,n_i,iat,icav,lai,rr,kt,ibast,par_grat,real(ei));
%'passo'
%end
gov(k)=g0;
eiv(k)=ei;
end
alf0=imag(eiv);
gai0=real(eiv);

der=diff(alf0);
ne=length(find(der<=0));
po=length(find(der>0));
if po>=1 | ne>=1
 if po>ne
  puf=(find(der>0)); 
 else
  puf=(find(der<=0)); 
 end
 isw=1
else
 isw=0
end
clear PU
 if isw==1
 fib=find(diff(puf)>1);
   
 puint=1:length(puf);  
  ipu=0;
  if  length(fib)>0
  puint=fib+1:length(puf);
  pud=puf(1:fib);
  adu=alf0(pud);
  se=adu(1:end-1).*adu(2:end);

  if length(find(se<0))==1
  ipu=ipu+1;
  PU{ipu}=pud;
  end
  end % fib
  pud=puf(puint);
  adu=alf0(pud);
  se=adu(1:end-1).*adu(2:end);
  if length(find(se<0))==1
   ipu=ipu+1;
   PU{ipu}=pud;
  end
%    ' PU', keyboard 
else 
 PU{1}=1:length(lav);
end 

%    ' PU', keyboard 

for kint=1:length(PU)
 pui=PU{kint};
 [du,imid] =min(abs(alf0(pui)));
 imi=pui(imid);
  if ifp==-10
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
  [ei,g0]=Flam_gamLu(L_i,n_i,iat,icav,lai,rr,kt,ibast,par_grat);
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

if ifp==-10
figure, plot(lavf,alvf,LAM,0,'wo'), pausak
end
cog=polyfit(lavf,gvf,2);
GTH=polyval(cog,LAM);
G0=GTH;
cog0=polyfit(lavf,govf,2);
Gga=polyval(cog0,LAM);

if ifp==-10
h=figure, 
set(h,'pos',PO)
subplot(211)
plot(lav1,eiv,LAM,GTH,'wo'), title('real part'), grid
subplot(212)
plot(lav1,imag(eiv),LAM,0,'wo'), title('imag part'), grid
pausak
end
 GV(kint)=G0;
 LV(kint)=LAM;

end   %fine kint

%' qui cont fine kint', keyboard

if ifp==-10
GV
[G0,imi]=min(GV);
G0
LAM=LV(imi)
'inizio iterazione', pausak
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
cg=polyfit(Gx,Gve,1);
x=linspace(min(Gx),max(Gx),10);
cg0=cg;
cg0(1)=cg0(1)-1;
gg0=roots(cg0);
Gstim=gg0;
ga=Gstim;
la=LAM;

if ifp==-10
figure, plot(x,polyval(cg,x),gg0,gg0,'ro')
pausak
end

if i_campi==1

if ifp==-10
' campi,' , pausak
end

%' campi,' , keyboard
Nx=50;
[fiez,nz,zet,Gaqw,NQW_ef,az,uL]=Flam_field(L_i,n_i,iat,fiQW,la,rr,kt,ibast,par_grat,ga,Nx);

%[fiez,nz,zet,Gaqw,NQW_ef,az]=Flam_field(L_i,n_i1,iat,fiQW,la,rr,kt,ibast,par_grat,ga,Nx,a_i);


Ez=sum(fiez);
Ez=Ez/max(Ez);
Hz=diff(fiez);
Hz=Hz/max(Hz);

lab='';
if isstr(par_grat)==1
if par_grat.itetm==1
 lab='TE';
else
 lab='TM';
end
end

if ifp==-10
figure, plot(zet,abs(Ez).^2*3, zet, nz,'r'), 
 title(['CMM: ',lab,'  lambda_{res} = ',num2str(la),' Gth pa= ',num2str(ga/NQW_ef)]), 
'fine campi', pausak
end

else
 zet=0;
 nz=0;
 Ez=0;
 Hz=0;
end

%'fine CMM', keyboard