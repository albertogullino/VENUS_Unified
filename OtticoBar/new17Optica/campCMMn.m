clear all
clear global
close all

icloo=1;
icloo=0;
iret_BW=0;

teta=4;
%teta=0;

ifp=-11;
ifp=-4;

itetm=3;  % 3 tutte, 1-TE, 2-TM
%itetm=1;  % 3 tutte, 1-TE, 2-TM

nomeF='mims_okn.str';
nomeF='HCG_VDAYflat.str';
nomeF='hcg_Ber0.str';   
%nomeF='mims_strain_simp.str';
%nomeF='bottouCor.str';
%nomeF='preflex_botto1.str';
%nomeF='preflex_botto2.str';
%nomeF='preflexp.str';    % no reticolo
%nomeF='giel_paper3D_ret0.str';
%' prima', pausak
[nref,aref,nv0,dv,xm,Dov,radii,fst,perd,anyf,iauto,dvw,xmw,fsw,dov,shavet,...
   ipar,ifield,lambda0,igrad_tr]=Lay_tran1D(nomeF);
%   pua=find(iauto(:,1)==2);

%' dopo', pausak

MOSTRO_SOLUZIONE=0;
rr=nref;
iLP=0;
iref=1;  
ks=1;



lambda=lambda0*1e6;

R0=0;
nitot=nv0;   
nv=nv0;
%'dop', keyboard 

fs=.98;
fs=1;

period0=.4*fs;
DC=.6;
dret=260*fs;

%sol0_sub
fish=find(shavet==6);
if length(fish)>0  
  ngrating=nv(fish,:);
  par_grat.r1=ngrating(1);
  par_grat.r2=ngrating(2);
  par=radii.array{fish};
  n1g=ngrating(1);
  n2g=ngrating(2);
  period=period0;
  par_grat.itetm=itetm;
  par_grat.per=period;
  par_grat.DC=DC;
  NModi=11;
  par_grat.NModi=NModi;
  if ~exist('itetm')
  itetm=3;
  end
  icarico=0;   % calcola T reticolo caricato sui G_i
%  dret=dv(fish);
  %' cont grat', keyboard
  ibast=fish;
  ivap=ibast+1;
par_grat.iret_BW=iret_BW;
par_grat.itetm=itetm;
par_grat.r_in=(nv(ibast-1,1));
par_grat.r_out=(nv(ibast+1,1));
else 
 par_grat1=[];
 ibast=[];
 dret=0;
end 

dlam=5e-2;
NPlam=20;



kref=2*pi/lambda*rr;
kout=2*pi/lambda*par_grat.r_out;
kt=kout*sin(teta/180*pi)/kref;
%' prima di kt', keyboard
%kt=0;

% scompatto@!
unpack_0

fiCav=icav;
fiQW=iat;

ibast=find(shto==6);

fiNp=find(fst(:,2)==4);

ivar=ibast+1;

ivar=1;

par=700+[0:.05:2.5]*1000;
par=700+[0:.05:.9]*1000;
%par=700+[0:50:900];
pari=linspace(250,270,11);
parip=linspace(110,510,51);
Np=[1:2:9];
Np=[1];
%par=[1200];
ipari=1;   %dret
%ipari=2;   %DC
%ipare=3;   %period
%ipari=4;   %last layer
%ipari=5;   %paia uscita
%ipari=0;   %last layer


ipare=1;   %dret
ipare=2;   %DC
ipare=3;   %period
%ipare=4;   %last layer
%ipare=0;   %last layer
%ipare=5;   %last layer

pare=0;
pret=[220:5:250]*fs;
pdc=[6 7]/10;
pper=[.38 .40 .42];
pper=(linspace(350,500,51))/1000;   %periodo reticolo
playi=parip;

preti=[250:10:330]*fs;
preti=[220:5:260]*fs;
preti=linspace(230,280,15);  
%preti=linspace(225,265,11)*fs;
pdci=linspace(.65,.75,11);
pperi=linspace(.38,.42,11)*fs;
play=[60 110 200];


  par_grat.per=period;
  par_grat.DC=DC;
  
 if ipare==1
  pare=pret;
 elseif ipare==2
 pare=pdc;
 elseif ipare==3
  pare=pper;
 elseif ipare==4
  pare=play;  
 elseif ipare==5
  pare=Np;    
 end

 if ipari==1
  pari=preti;
  lab='Thickness';
 elseif ipari==2
 pari=pdci;
  lab='DC';
 elseif ipari==3
  pari=pperi;
  lab='Period';
 elseif ipari==4
  pari=playi; 
  lab='Oxide thick';
 elseif ipari==5
  pari=Np;    
  lab='Pairs';
 end
' cont', keyboard
tic

global ie iloop
for ie=1:length(pare)
%for ie=4
  parei=pare(ie);
 if ipare==1
  L_i(ibast)=parei;
 elseif ipare==4
  L_i(ivap)=parei;
 elseif ipare==2
  par_grat.DC=parei;
 elseif ipare==3
  par_grat.per=parei;
 elseif ipare==5
   fst(fiNp,2)=parei;
   unpack_0
   fiCav=icav;
   fiQW=iat; 
 end

for iloop=1:length(pari)
%for iloop=13

parii=pari(iloop)
 if ipari==1
  L_i(ibast)=parii;
 elseif ipari==4
  L_i(ivap)=parii;
 elseif ipari==2
  par_grat.DC=parii;
 elseif ipari==3
  par_grat.per=parii;
 elseif ipari==5
   fst(fiNp,2)=parii;
   unpack_0
   fiCav=icav;
   fiQW=iat; 
 end 
 
if icloo==1 
' controllo loop',
[ie iloop], 
[parei parii], pausak
else
 if itetm==2 | itetm==3
  'TM'
  par_grat.itetm=2;
  [gam,lam]=fiezCMMn(lambda,L_i,n_i,fiQW,fiCav,fiQ,ifp,dlam,NPlam,rr,kt,ibast,par_grat,ato);
  gmv(iloop,ie)=gam;
  lmv(iloop,ie)=lam;
 end 


%[gam,lam]=fiezCMM(lambda,L_i,n_i,fiQ,fiCav,ifp,dlam,rr,kt,ibast,par_grat1);

 if itetm==1 | itetm==3
   'TE'
 par_grat.itetm=1;
  [gae,lae]=fiezCMMn(lambda,L_i,n_i,fiQW,fiCav,fiQ,ifp,dlam,NPlam,rr,kt,ibast,par_grat,ato);
  gev(iloop,ie)=gae;
  lev(iloop,ie)=lae;
 end
if itetm==1
  gmv(iloop,ie)=NaN;
  lmv(iloop,ie)=NaN;
elseif itetm==2
  gev(iloop,ie)=NaN;
  lev(iloop,ie)=NaN;
end
end %if icloo
% pausak
end %ii
end %ie

toc
par=pari;
h=figure, 
set(h,'pos',[150    400   1000   400])
subplot(131)
plot(par,lmv), hold on, plot(par,lev,'--')
a=axis;
a(1:2)=[min(par) max(par)];
axis(a)
subplot(132)
plot(par,gmv), hold on, semilogy(par,gev,'--')
%semilogy(par,gmv), hold on, semilogy(par,gev,'--')
a=axis;
a(3:4)=[100 300];
a(1:2)=[min(par) max(par)];
axis(a)
subplot(133)
plot(par,(1-gmv./gev)*100)
a=axis;
a(1:2)=[min(par) max(par)];
axis(a)
title(' dichroism ')
fsr=(max(lmv)-min(lmv))*1000

pausak

colordef white
if itetm==3 | itetm==2
h=figure, 
set(h,'pos',[150    300   1000   400])
subplot(121) 
contourf(pari,pare,-log10(gmv'))
 xlabel(' HCG thickness (nm)')
 ylabel(' HCG period (um)')
 title(' -log10(Gth) (1/cm)')
 colorbar('Southoutside')
 subplot(122) 
 contourf(pari,pare,1000*(lmv'))
 xlabel(' HCG thickness (nm)')
 ylabel(' HCG period (um)')
 colorbar('Southoutside')
 title(' wavelenght (nm)')
pausak
end
if itetm==3 | itetm==1
h=figure, 
set(h,'pos',[150    200   1000   400])
subplot(121) 
contourf(pari,pare,-log10(gev'))
 xlabel(' HCG thickness (nm)')
 ylabel(' HCG period (um)')
 title(' -log10(Gth) (1/cm)')
 colorbar('Southoutside')
 subplot(122) 
 contourf(pari,pare,1000*(lev'))
 xlabel(' HCG thickness (nm)')
 ylabel(' HCG period (um)')
 colorbar('Southoutside')
 title(' wavelenght (nm)')
pausak
end



h=figure;
set(h,'pos',100*[3 2 9 5])
subplot(121)
plot(par,lmv*1000,'linewidth',2),
ylabel('wavelength (nm)')
xlabel(lab)
a=axis;
a(1:2)=[min(par) max(par)];
axis(a)
subplot(122)
semilogy(par,gmv,'linewidth',2),
xlabel(lab)
ylabel('Gain/QW (1/cm)')

a=axis;
a(3:4)=[100 1000];
a(1:2)=[min(par) max(par)];
axis(a)
pausak
figure
semilogy(par,gmv)
a=axis;
%a(3:4)=[100 400];
%a(1:2)=[min(par) max(par)];
%axis(a)