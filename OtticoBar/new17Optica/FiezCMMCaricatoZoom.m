function [ga,la,zet,Ez,nz,Gaqw,NQW_ef,az]=fiezCMMnZ(lam0,L_i,n_i1,iat,icav,fiQW,ifp,dlam,NPlam,rr,kt,ibast,par_grat,a_i);

Gprec=1e-4

global Options
iAlb=Options.iAlb;
IC0=15;
if iAlb==1
 IC0=0;
end
%' in fiez vecchio ', keyboard

n_i=n_i1(:,1);
ifpc=-4;  %controllo soluzione
ifpc=-10;  %controllo soluzione
iplotCa=0;
%iplotCa=1;
ifp=-11
if ifp<-10
ifpc=-11;  %controllo soluzione
end
if ifp<=-10
iplotCa=1;
end
%'iplotaca', keyboard
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
if ifp==-12
figure, plot(L_i,'o'), hold on, plot(iat,L_i(iat),'rp')
a=axis;
a(3)=0;
a(4)=max(L_i)+10;
axis(a)
pausak
end

lav=lam0+linspace(-1,1,5*NPlam)*dlam*lam0;
%lav=lam0+linspace(0,2,20)*dlam*lam0;
%lav=lam0;

%' primo Caricato', keyboard
for k=1:length(lav)
lai=lav(k);

[ei,g0,T22]=Flam_gamCaricato(L_i,n_i,iat,icav,lai,rr,kt,ibast,par_grat);
tv(k)=T22;
%'dopo Flam Caricato', keyboard
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
 afiz=(abs(alf0(fiz))+abs(alf0(fiz+1)))/2;
 [du,imi]=min(gafiz+afiz);
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
der=alf0;
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
'trovati punti', keyboard
else
 if fiz>1
  if fiz<length(gai0)-2
   puv=fiz+[0 1 2];
   puv=fiz+[0 1];
  else
   puv=sort(length(gai0)-[0 1 2]);
%   puv=sort(length(gai0)-[0 1]);
  end
 else
  puv=fiz+[0 1 2];
  puv=fiz+[0 1];
 end
end

coa=polyfit(lav(puv),alf0(puv),1);
l0=roots(coa);
if ifpc==-10
 figure
  subplot(211)
  plot(lav,alf0,lav(puv),alf0(puv),'r.'), grid
  subplot(212)
  semilogy(lav,gai0,lav(puv),gai0(puv),'r.'), grid, pausak
  
end

%'trovati punti', keyboard
G0=mean(gai0(puv));
if G0>1e6
 ga=NaN;
 la=NaN;
 zet=NaN;
 Ez=NaN;
 nz=NaN;
 Gaqw=NaN;
 NQW_ef=0;
 az=0;
 uL=0;
 
 'NaN in FiezCMMCaricato', keyboard
 return
 %ancora
end

%'trovati punti', keyboard

 PU{1}=puv;

clear LV GV
kisol=0;
for kint=1:length(PU)
 pui=PU{kint};
 if length(pui)>1
 kisol=kisol+1;
 [du,imid] =min(abs(alf0(pui)));
 imi=pui(imid);
lav1=lav;
%imi=1;
fiz=[imi imi];
Npz=3;

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
 npun=20;
end
 npun=Npz;
imi=1;
ic=0
dlau=diff(lav(1:2));

dii=2;
Npz=3;
%dlau=diff(lav1(1:2));
ibre=0;


%G0=mean(gai0(puv));
%LAM=mean(lav(puv));

cG0f=polyfit(lav(puv),gai0(puv),1);
cLaf=polyfit(lav(puv),alf0(puv),1);
LAM=roots(cLaf);
G0=polyval(cG0f,LAM);
%'qui G0', keyboard

 icsog=1;
 lav1=mean(lav(puv))+linspace(-dlau/2,dlau/2,npun);
%' primo lav1', keyboard
lac=LAM;
dii=diff(lav1(1:2));
while ic<IC0
 if ibre==1
  break
 end
 if npun==Npz
  %ibre=1;
 end 
 ic=ic+1;
 dlau=diff(lav1(1:2));
 if ic>icsog
%  lav1=lac+linspace(-dii/2,dii/2,npun)*dlau;
  lav1=lac+linspace(-1/2,1/2,npun)*dlau;
  %lav1=lac+linspace(-dii/2,dii/2,npun)*dlau;
 end 
% ' qui', keyboard
 clear gov eiv
 for k=1:length(lav1)
  lai=lav1(k);
  [ei,g0,T22]=Flam_gamCaricato(L_i,n_i,iat,icav,lai,rr,kt,ibast,par_grat,G0);
  tov(k)=T22;
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
%if ic<2
%[du,fiz]=min(abs(alf)+10*abs(gav));
%else
[du,fiz]=min(abs(alf));
%end
end


  npun=Npz;
  dii=1;
  dlau=diff(lav1(1:2));
  lac=lav1(fiz(imi));

  
%' lac', keyboard


  puv=1:3;
  cG0f=polyfit(lav1(puv),gav(puv),2);
  cLaf=polyfit(lav1(puv),alf(puv),2);
  
  LAMd=roots(cLaf);
  [du,idu]=min(abs(LAMd-LAM));
  LAM=LAMd(idu);
  G0i=polyval(cG0f,LAM)
%  figure, 
%  subplot(121)
%  plot(lav1,alf,LAM,0,'ro'), grid
%  subplot(122)
%  plot(lav1,gav,LAM,G0i,'ro'), grid


  G0=G0+G0i;
  guess=[LAM G0];
  if ifp==-11
   T_zero(guess,L_i,n_i,iat,icav,rr,kt,ibast,par_grat)
  end 
  if abs(G0i)<Gprec
   ibre=1;
  end
  if ifpc==-11
  figure,  
  subplot(211)
  plot(lav1,alf,'.-',LAM,0,'ro'), grid
    subplot(212)
  plot(lav1,gav,LAM,G0i,'ro'), grid
  pausak
  end  
  
 %[du,imi]=min(gav(fiz)); 
 imi=1;
 
 

 end
 %' fine while ', pausak
 
       if iAlb==0											
	 GV(kisol)=G0;
	 LV(kisol)=LAM;
      end	 %iOLD
 end %length(pui)
end   %fine kint


if iAlb==0
 if ifpc==-10
	GV
	[G0,imi]=min(GV);
	G0
	LAM=LV(imi)
%	'inizio iterazione', pausak
 end
end 

if ifpc==-11
 figure
  subplot(211)
  plot(lav,alf0,LAM,0,'ro'), grid
  subplot(212)
  semilogy(lav,gai0,LAM,G0,'ro'), grid, pausak
end


if iAlb==1

igue=0

guess=[LAM G0];
if igue==1
 T_zeroFerma(guess,L_i,n_i,iat,icav,rr,kt,ibast,par_grat)
 guess1=[lav(puv(1)) G0];
 ' prima guess 1', keyboard
 T_zeroFerma(guess1,L_i,n_i,iat,icav,rr,kt,ibast,par_grat)
'fine guess1', pausak
 guess2=[lav(puv(2)) G0];
 T_zeroFerma(guess2,L_i,n_i,iat,icav,rr,kt,ibast,par_grat)
'fine guess2', pausak

lp=linspace(lav(puv(1)),lav(puv(2)),101);
lp=linspace(3.24585,3.2459,11);
lg=linspace(G0-30,G0,11);
clear tv
for kl=1:length(lp)
for kg=1:length(lg)
  guess1=[lp(kl) lg(kg)];
 tv(kl,kg)=T_zeroFerma(guess1,L_i,n_i,iat,icav,rr,kt,ibast,par_grat)
end 
end

end

OPTIONS = optimset('Display','Iter','TolFun',2,'MaxFunEvals',100);
OPTIONS = optimset('Display','Iter','TolFun',2);
OPTIONS = optimset('Display','Iter','TolX',1e-3,'MaxFunEvals',70,'TolFun',.5);
'vetop', keyboard
[X,FVAL,EXITFLAG,OUTPUT]=fminsearch(@(vet)T_zero(vet,L_i,n_i,iat,icav,rr,kt,ibast,par_grat),guess,OPTIONS)
ga=X(2);
la=X(1);
' dopo options', keyboard
 T_zeroFerma(X,L_i,n_i,iat,icav,rr,kt,ibast,par_grat)



else

%' qui T_zero', keyboard

X(1)=LAM;
X(2)=G0;
Res=T_zero(X,L_i,n_i,iat,icav,rr,kt,ibast,par_grat,ifp);
if Res>log10(Gprec*1000)
 'problemi ricerca soluzione ', %keyboard
end
%[Mtot,Mtot1]=T_zeroFerma(X,L_i,n_i,iat,icav,rr,kt,ibast,par_grat)

la=LAM;
ga=G0;

end  % iOLD

%' campi,' , pausak
%'qui nz', keyboard
Nx=50;
%Nx=2;

[fiez,nz,zet,Gaqw,NQW_ef,az]=Flam_fieldCaricato(L_i,n_i1,iat,fiQW,la,rr,kt,ibast,par_grat,ga,Nx,a_i);
%[fiez,nz,zet,Gaqw,NQW_ef,az]=Flam_field(L_i,n_i1,iat,fiQW,la,rr,kt,ibast,par_grat,ga,Nx,a_i);

Ez=sum(fiez);
Ez=Ez/max(Ez);
Hz=diff(fiez);
Hz=Hz/max(Hz);

fize=find(nz==0);
nz(fize,1)=NaN;

gatot=ga;
ga=gatot/NQW_ef;
if length(par_grat)>0
if par_grat.itetm==1
 lab='TE';
else
 lab='TM';
end
else
 lab='';
end 


if iplotCa==1
zqw=sum(L_i(1:fiQW(1)))/1000;
[du,imi]=min(abs(zet-zqw));
I0=abs(Ez).^2;
I=log10(I0/I0(imi))+3;
I(1,fize)=NaN;
h=figure, 
set(h,'pos',[150         127        1367         523])
subplot(121)
plot(zet,I, zet, nz(:,1),'r'), 
 title(['CMM: ',lab,'  lambda_{res} = ',num2str(la),' Gth= ',num2str(ga)]), 
  xlabel(' Long. coord. (um)')
 ylabel(' Intens (Log10) / Index')
 a=axis;
 a(3)=-1.5;
 a(4)=6;
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
