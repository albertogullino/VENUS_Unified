function [ga,la,zet,Ez,nz,Gaqw,NQW_ef,az,uL,E_Temp,z_Temp,velm1D]=fiezCMMn(lam0,L_i,n_i1,iat,icav,fiQW,ifp,dlam,NPlam,rr,kt,ibast,par_grat,a_i);

Gprec=1e-5;

global Options
iAlb=Options.iAlb;
IC0=15;
if iAlb==1
 IC0=0;
end

velm1D=0;

kNoll=2*pi/(lam0*1e-6);
n_i=n_i1(:,1);
ifpc=-4;  %controllo soluzione
% ifpc=-10;  %controllo soluzione
iplotCa=0;
%iplotCa=1;
%'qopt', keyboard

%ifp=-11        % da scommentare per vedere passo passo

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

%' in fiez vecchio ', keyboard

lav=lam0+linspace(-1.,1.,NPlam)*dlam*lam0;
%lav=lam0+linspace(-1,0,NPlam)*dlam*lam0;
%lav=lam0+linspace(0,2,20)*dlam*lam0;
%lav=lam0;

%' primo Caricato', keyboard
icav=icav([1 end]);
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
if ifpc<=-11
 figure, subplot(211)
 plot(lav,alf0,lav(fiz),alf0(fiz),'ro'), grid
 subplot(212)
 semilogy(lav,gai0,lav(fiz),gai0(fiz),'ro'), pausak
end
if length(fiz)>0
 if fiz(end)==length(lav)
  fiz(end)=length(lav)-1;
 end
 gafiz=(gai0(fiz)+gai0(fiz+1))/2;
 afizd=(abs(alf0(fiz))+abs(alf0(fiz+1)))/2;
 afiz=afizd/mean(afizd)*mean(gafiz);   %normalizzo afiz
 [du,imi]=sort(gafiz+afiz);
 if length(imi)>1;
  if length(imi)<=3
   fizu=fiz(imi(1:end-1));
  else
     fizu=fiz(imi(1:end-2));
  end
if ifpc<=-11
 figure, subplot(211)
 plot(lav,alf0,lav(fizu),alf0(fizu),'ro'), grid
 subplot(212)
 semilogy(lav,gai0,lav(fizu),gai0(fizu),'ro'), grid
 pausak
end  
  clear gafizu l0v
  kfok=0;
  for kf=1:length(fizu)
   puk=fizu(kf)+[0 1];
   izo=1;
   if izo==0
    coa=polyfit(lav(puk),alf0(puk),1);
    l0=roots(coa); 
    l0v(kf)=l0;
    cog=polyfit(lav(puk),gai0(puk),1);   
    gafizu(kf)=polyval(cog,l0);
   else 
    laiv=linspace(lav(puk(1)),lav(puk(2)),5);
    clear alf0i gai0i
    for kff=1:length(laiv)
     lai=laiv(kff);
     [ei]=Flam_gamCaricato(L_i,n_i,iat,icav,lai,rr,kt,ibast,par_grat);
     alf0i(kff)=imag(ei);
     gai0i(kff)=real(ei);     
    end 
    coa=polyfit(laiv,alf0i,2);
    coa0=polyfit(laiv,alf0i,1);
    ai0=polyval(coa0,laiv);
    %picco=sum(abs(ai0-alf0i))/mean(abs(ai0))
    picco=max(abs(ai0-alf0i))/median(abs(ai0));
    l0vi=roots(coa); 
    [dl0i,ila]=min(abs(l0vi-mean(laiv)));
    l0=l0vi(ila);
    piccov(kf)=picco;
    cog=polyfit(laiv,log10(gai0i),2);    
    g0f=10^polyval(cog,l0);    

    if g0f>0 & picco<.8
     kfok=kfok+1;
     gafizu(kfok)=g0f;    
     l0v(kfok)=l0;
     fizus(kfok)=fizu(kf);
    end 
if ifpc<=-11
 figure, subplot(211)
 plot(lav,alf0,laiv,alf0i,'.-',l0,0,'ro'), grid
 subplot(212)
 semilogy(lav,gai0,laiv,gai0i,'.-',l0,g0f,'ro'), grid
 pausak
end     
   end
  end
if ifpc<=-11
 figure, subplot(211)
 plot(lav,alf0,lav(fiz),alf0(fiz),'ro',l0v,zeros(size(l0v)),'gs'), grid
 subplot(212)
 semilogy(lav,gai0,lav(fiz),gai0(fiz),'ro',l0v,gafizu,'gs'), grid, pausak
end  
  [dus,imu]=sort(gafizu/min(gafizu));  
  fimi=find(dus<1.5);
  fiz=fizus(imu(fimi));
 else
%  [du,fiz]=min(gafiz);  
%  fiz=fiz(imi);
 end 
end
der=diff(alf0);



 for km=1:length(fiz)
  puv=fiz(km)+[0 1];
  PU{km}=puv;
 end

coa=polyfit(lav(puv),alf0(puv),1);
l0=roots(coa);
if ifpc<=-10
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


clear LV GV
kisol=0;
for kint=1:length(PU)
 pui=PU{kint};
 puv=pui;
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
ic=0;
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

if ifpc<=-11
 ' while ini ', keyboard
end
while ic<IC0
 if ibre==1
%  'prima di brek', keyboard
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
%   'vedi uL', keyboard  
  if abs(LAM-lav1(2))>2*diff(lav1([1 end]))
   [ga,la,zet,Ez,nz,Gaqw,NQW_ef,az]=FiezCMMCaricatoZoom(lam0,L_i,n_i1,iat,icav,fiQW,ifp,dlam,NPlam,rr,kt,ibast,par_grat,a_i);
%   'vedi uL', keyboard
   return
  end
  
  G0i=polyval(cG0f,LAM);
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
 if ifpc<=-11
  ' fine while int ', keyboard
 end
 
       if iAlb==0											
	 GV(kisol)=G0;
	 LV(kisol)=LAM;
      end	 %iOLD
 end %length(pui)
end   %fine kint

if ifpc<=-11
 ' fine while ', keyboard
end 

if iAlb==0
% if ifpc==-10
	[G0,imi]=min(GV);
	LAM=LV(imi);
%	'inizio iterazione', keyboard
% end
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
if Res>log10(Gprec*100)
 'problemi ricerca soluzione ', keyboard
end
%[Mtot,Mtot1]=T_zeroFerma(X,L_i,n_i,iat,icav,rr,kt,ibast,par_grat)

la=LAM;
ga=G0;

end  % iOLD


%' campi,' , pausak
%'qui nz', keyboard


%Nx=2;

iold=1;
if iold==1
 GV=ga;
 LV=la;
end

global Ps

if isfield(Ps,'i1D')
 i1D=Ps.i1D;
else
 i1D=0;
end 

%'qui caricato', keyboard

for kcam=1:length(GV)
 la=LV(kcam);
 ga=GV(kcam);

%kcam=1;
Nx=50;

if i1D==0
 [Ez,Hz,nz,zet,Gaqw,NQW_ef,az,uL]=Flam_fieldCaricato2016Mod(L_i,n_i1,iat,fiQW,la,rr,kt,ibast,par_grat,ga,Nx,a_i);
else
 [Ez,Hz,nz,zet,Gaqw,NQW_ef,az,uL,velm1D]=Flam_fieldCaricato2021gullino(L_i,n_i1,iat,fiQW,la,rr,kt,ibast,par_grat,ga,Nx,a_i);
end
%[fiez,nz,zet,Gaqw,NQW_ef,az]=Flam_field(L_i,n_i1,iat,fiQW,la,rr,kt,ibast,par_grat,ga,Nx,a_i);
%fie,nz,zet,Gaqw,NQW_ef,az,uL
 Ezv{kcam}=Ez;
 Hzv{kcam}=Hz;
 Gaqwv(kcam)=Gaqw;
 Gatot(kcam)=ga/NQW_ef;
 NQv(kcam)=NQW_ef;
 Uv{kcam}=uL;
 zv{kcam}=zet;
 Nzv{kcam}=nz;

end

ifitu=0;
if ifitu==1
figure, plot(zet,Ez,zet,Hz), 
title('parti reali'), pausak



figure, plot(zet,imag(Ez),zet,imag(Hz)), 
title('parti immaginarie'), pausak

%figure, plot(zet,Ez.*conj(Hz)/700+4,'r',zet,nz), pausak

Z0=377;
Zrn=1/rr;
Intensityz=abs(Ez).^2;
Pointingz=real(Ez.*conj(Hz))/Zrn;
figure, plot(zet,Intensityz), 
title('Intensità'), pausak
figure,  plot(zet,Pointingz), pausak
Teff_1D=abs(Pointingz(1))/max(Intensityz)

figure, plot(zet,real(Ez).*real(Hz),'g',zet,-imag(Ez).*imag(Hz),'r',zet,nz), pausak
%figure, plot(zet,real(Ez).*real(Hz),'g',zet,-imag(Ez).*imag(Hz),'r',zet,nz), pausak
%figure, plot(zet,real(Ez).*real(Hz),'r',zet,nz), pausak
end


% valori medi

Lsu=cumsum(L_i)/1000;

% fiM=find(real(n_i)>2.9);
fiM=1:length(L_i);
zte=linspace(Lsu(fiM(1)),Lsu(fiM(end)-1),201);

E2=abs(Ez.^2);
dE=[diff(E2) 1];
dS=[diff(dE) 1];
pp=[dE(1:end-1).*dE(2:end) 1];
fiM=find(pp<0 & dS<0 & zet>Lsu(1) & zet<Lsu(end-1));
fi1=find(zet<Lsu(1));
fi2=find(zet>Lsu(end-1));
fiM=sort([fiM fi1(end) fi2(1)]);

%nzM=nz(fiM);
%fz=find(nzM<1.7);
%fiMs=fiM;
%fiM=fiM(fz(end)+1:end);

%coe=polyfit(zet(fiM),log10(E2(fiM)),17);
%Em=10.^polyval(coe,zte);
 zf=zet(fiM);
 Ef=E2(fiM);

EmL=interp1(zf,log(Ef),zte,'linear',0);
Em=exp(EmL);

if ifp==-10
figure, semilogy(zet,E2,zte,Em,'r',zet(fiM),E2(fiM),'r.'), 
ylim([1e-2 max(E2)*1.2])
pausak
%figure, plot(zet,E2,zte,Em,'r',zet(fiM),E2(fiM),'r.'), pausak
end
E_Temp=Em/max(Em);
% z_Temp=zte-Lsu(fiM(1));
z_Temp=zte-zte(1);
%coe=polyfit(zet,(E2),15);
%Efl=polyval(coe,zte);
%figure, plot(zet,E2,zte,(Efl),'r')

 %Ezv{kcam}=Ez;
 %Hzv{kcam}=Hz;
 %Gaqwv(kcam)=Gaqw;
 %Gatot(kcam)=ga/NQW_ef;
 %NQv(kcam)=NQW_ef;
 %Uv{kcam}=uL;

[gas,fis]=sort(Gatot);

isce=1;
fimo=fis(1);
if length(fis)>1
 fimo=fis(isce);
end 


la=LV(fimo);
NQW_ef=NQv(fimo);
Ez=Ezv{fimo};
nz=Nzv{fimo};
zet=zv{fimo};
Hz=Hzv{fimo};
uL=Uv{fimo};

fize=find(nz==0);
nz(fize,1)=NaN;
if isfield(par_grat,'itetm')==1 
if par_grat.itetm==1
 lab='TE';
else
 lab='TM';
end
else
 lab='';
end 


%Po=real(Ez.*conj(Hz));
%figure, plot(zet,Po)
% iplotCa=1;
if iplotCa==1 
    zqw=sum(L_i(1:fiQW(1)))/1000;
    [du,imi]=min(abs(zet-zqw));
    I0=abs(Ez).^2;
    I=log10(I0/I0(imi))+max(real(nz(:,1)));
    I(1,fize)=NaN;
    
    figure
    vg=uL(2)/uL(1);
    set(gcf,'pos',[150         127        1367         523])
    subplot(121)
    plot(zet,I, zet, nz(:,1),'r'),
    title(['CMM: ',lab,'  lambda_{res} = ',num2str(la),' Gth= ',num2str(ga/vg)]),
    xlabel(' Long. coord. (um)')
    ylabel(' Intensity (Log_{10}) / Index')
    a=axis;
    a(3)=-1.5;
    a(4)=max(real(nz(:,1)))+.5;
    axis(a)
    
    subplot(122)
    plot(zet,I0/max(I0)*max(real(n_i)), zet, nz(:,1),'r'),
    title(['CMM: ',lab,'  lambda_{res} = ',num2str(la),' Gth= ',num2str(ga/vg)]),
    xlabel(' Long. coord. (um)')
    ylabel(' Intensity  / Index')
    
    'fine campi', pausak
end
%'fine campi', pausak
%keyboard
