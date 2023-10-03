
n_sav=n_i;
n_i=nto;
if dret>0 & imet<2
n_i(fiirv)=nreticolo;
end
%' in guessf', keyboard
iscattu=1;
if ifp==-10 & icheck_guess==1
 iplot=1
else
 iplot=0;
end
isub_old=1;

Ldum=L_i;
 cs=cumsum(Ldum);
% z_cav=cs(fiCav);
% z_cav(1)=z_cav(1)-100;
% z_cav(end)=z_cav(end)+100;


iproga=0;
iff=0;
freq=0;
if ks==1
 lambdau=lambda_cen;
 nl0=61;
 famul= NP_dec;
 str=str0;
else
 str=str1;
 nl0=61;
 famul=1;
end
nl=famul*nl0;

lambda_cenmod=lambdau;
k02=4e4*pi/lambdau;


      

dlav=linspace(-str*7/5,str*7/5,nl);
lav=lambdau+dlav;


%nv0=nv(:,1);
fiQ=find(fsto==-1);
NQW=length(fiQ);
iqw=round(NQW/2);
fiQ01=fiQ(iqw);
if ~exist('kv')
 kv=0;
else
 kv
end 

 %'qui guessf: KV', keyboard

clear gae gam Lee Lem Gev Gmv s11e s11m s22e s22m fve fvm env
if ifp==-10
% 'qui guessf', keyboard
end 
for k=1:length(lav)

 lambdai=lav(k);
 [gazk,enan,fak,Tu,Tb,Gue,Gum,Gbe,Gbm,Lf,Lcav,ztot,Ez,Hz,indz,nmean,Perd,Ge,Gm,du,du,KKe,KKm]=...
         th_scattu(fiQ,fiCav,L_i,n_i,rr,rfd,rfu,lambdai,freq,kv,iLP,ifp,iff,ibast,par_grat);

 gae(k)=gazk(1);
 gam(k)=gazk(2);
 Lee(k)=Lf(1);
 Lem(k)=Lf(2);
 Gev(k)=Ge;
 Gmv(k)=Gm;
 s11e(k)=Gue;
 s11m(k)=Gum;
 
 Ke(k,:)=KKe(1,:); 
 Km(k,:)=KKm(1,:); 
 s22e(k)=Gbe;
 s22m(k)=Gbm;
 fve(k)=fak(1);

 if iscattu==1
  fvm(k)=fak(2);
 else
  fvm(k)=fak(1);
 end 
 env(k,:)=enan;
 %pausak
end

%'fine loop ', keyboard

%f11e=unwrap(angle(s11e));
%f11m=unwrap(angle(s11m));
%f22e=unwrap(angle(s22e));
%f22m=unwrap(angle(s22m));
f11e=angwr(s11e);
f11m=angwr(s11m);
f22e=angwr(s22e);
f22m=angwr(s22m);
NQW=length(fiQ);
iqw=round(NQW/2);
fi_qw0=fiQ(iqw);
dqw=L_i(fi_qw0);
nqw=n_i(fi_qw0);

gthe=-log(abs(s11e.*s22e))/dqw*1e7;  % in 1/cm
gthm=-log(abs(s11m.*s22m))/dqw*1e7;  % in 1/cm

%gthe=-log(abs(s11e.*s22e))./Lee*1e4./gae;  % in 1/cm
%gthm=-log(abs(s11m.*s22m))./Lem*1e4./gam;  % in 1/cm
%'cont gthw', keyboard
%  Gthp = 2e-2*pi*rr/lambdap*imag(Ksip)/Gamp;
fas_qw=4*pi./lav*real(nqw)*dqw*1e-3;
fa_tot0e=f11e+f22e-fas_qw;
fa_tot0m=f11m+f22m-fas_qw;
dfas=max(fa_tot0e)-min(fa_tot0e);
dfase= floor(mean(fa_tot0e)/2/pi)*2*pi;
nsol=ceil(dfas/2/pi)+4;
gth=gthe;
f11=f11e;
s11=s11e;
f22=f22e;
s22=s22e;
fa_tot0=fa_tot0e;
ene=env(:,1)';
eneu=ene.*fve;
%' ene ', keyboard

%TE
iso=0;
dfase=0;
clear dlv gtv gamv dlv0 enmv
for isol=1:nsol
fa_tot=fa_tot0e-(isol-1)*2*pi-dfase;
%[du,imi]=min(abs(fa_tot));
fa1=fa_tot(1:end-1);
fa2=fa_tot(2:end);
fap=fa1.*fa2;
%' stop figure', keyboard
imi=find(fap<=0)+1;
%imi
%keyboard
if length(imi)>0
imi=imi(1);
if imi==length(fap)
 imi=imi-1;
end
if imi>1  & imi<nl

fi=imi+[-1 0 1];
cof=polyfit(dlav(fi),fa_tot(fi),2);
dlavd=roots(cof);
[du,iro]=min(abs(dlavd-dlav(imi)));
dlav0=dlavd(iro);
cog=polyfit(dlav(fi),log10(gth(fi)),2);
gth0=10^polyval(cog,dlav0);
if isnan(gth0)==1
cog=polyfit(dlav(fi),(gth(fi)),2);
gth0=polyval(cog,dlav0);
end

coga=polyfit(dlav(fi),log10(fve(fi)),2);
gam0=10^polyval(coga,dlav0);
coen=polyfit(dlav(fi),log10(eneu(fi)),2);
en0=10^polyval(coen,dlav0);
if gth0>0
 iso=iso+1;
 enmv(iso)=real(en0); 
 dlv(iso)=real(dlav0);
 gtv(iso)=real(gth0);
 gamv(iso)=real(gam0); 
 dlv0(iso)=(isol-1)*2*pi;
% 'isol ', isol, pausak
end
end %if
end %length
end
%' fine loop gamma ', keyboard
%' fine loop gamma ', keyboard
if ~exist('gtv')
 ierrla=1;
' fine loop gamma, problemi ', keyboard 
 str0=str0*5;
 'set str', keyboard
 guessf_mod 
else
 ierrla=0;
end 


if itetmt==1 | dret==0 | itetmt==3
gtve=gtv;
dlve=dlv;
[GT0e,igm]=min(real(gtv));
L0e=real(dlv(igm));
Ga0e=real(gamv(igm));
En0e=real(enmv(igm));

if ks==1

global la1Ds ga1Ds 
 if length(gtv)>2
  if igm==1
   igm=igm+1;
  end
  if igm==length(dlv)
   igm=igm-1;
  end     
  igms=igm+[-1 0 1];
  la1Ds=lambda_cen+dlv(igms);
  ga1Ds=gtv(igms);
  [GT0e1,igm2]=min(real(gtv(igms)));
  L0e1=real(dlv(igms(igm2)));
  Ga0e1=real(gamv(igms(igm2)));
  gth1=GT0e1/(NQW*Ga0e1)/2;
  lambda1=lambda_cen+L0e1;
 else
  la1Ds=lambda_cen+dlv(igm);
  ga1Ds=gtv(igm);
  [GT0e1,igm2]=min(real(gtv));
  L0e1=real(dlv(igm2));
  Ga0e1=real(gamv(igm2));
  gth1=GT0e1/(NQW*Ga0e1)/2;
  lambda1=lambda_cen+L0e1;
 end
end

if ks==1
str1=mean(diff(dlv));
end
if iplot==1
figure, plot(dlav,f11,dlav,f22), title(' fase TE' ), pausak
figure, plot(dlav,abs(s11),dlav,abs(s22)),
title(' moduli S_ii TE'), pausak
figure, plot(dlav,fa_tot0,dlv,dlv0,'ro'), title(' cond. fase  TE'), pausak
%figure, semilogy(dlav,fve,dlv,gamv,'wo',L0e,Ga0e,'gx'), pausak
%figure, semilogy(dlav,eneu,dlve,enmv,'wo',L0e,En0e,'gx'), pausak
figure, semilogy(dlav,gth,dlv,gtv,'ro',L0e,GT0e,'gx'), title(' cond modulo TE'), 
pausak
figure, semilogy(1000*lav,gth,1000*(lambda_cen+dlv),gtv,'ro',...
                    1000*(lambda_cen+L0e),GT0e,'gx'), title(' cond modulo TE'), 
grid, pausak
end


GteQW=GT0e/(NQW*Ga0e)/2




la_ver=lambdau+L0e;
nimm=la_ver*GteQW/(4*pi)*1e-4;
%n_i(fiQ)=n_i(fiQ)+j*nim;
lambda=la_ver;





%'cont COE', keyboard
  Gthp=real(GT0e)
  Lvpi=real(la_ver)
  nime=Gthp/(4e4*pi)*Lvpi;
  G0e=Gthp;
  lime=L0e;
  
  gthe=gth;
  gtve=gtv; 
  dlve=dlv;
  gamev=gamv;
  enmev=enmv;

end  %itetmt

if dret>0  & itetmt>=2
%TM
dfas=max(fa_tot0m)-min(fa_tot0m);
dfasm= floor(mean(fa_tot0m)/2/pi)*2*pi;
dfasm=0;
nsol=ceil(dfas/2/pi)+1;
gth=gthm;
f11=f11m;
s11=s11m;
f22=f22m;
s22=s22m;
ene=env(:,2)';
eneu=ene.*fvm;
fa_tot0=fa_tot0m;
iso=0;
clear dlv gtv gamv dlv0 enmv
%' wui ', keyboard
for isol=1:nsol
fa_tot=fa_tot0-(isol-1)*2*pi-dfasm;
fa1=fa_tot(1:end-1);
fa2=fa_tot(2:end);
fap=fa1.*fa2;
imi=find(fap<0)+1;

if length(imi)>0
imi=imi(1);
if imi==length(fap)
 imi=imi-1;
end
if imi>1  & imi<nl
iso=iso+1;
fi=imi+[-1 0 1];
cof=polyfit(dlav(fi),fa_tot(fi),2);
dlavd=roots(cof);
[du,iro]=min(abs(dlavd-dlav(imi)));
dlav0=dlavd(iro);
cog=polyfit(dlav(fi),gth(fi),2);
gth0=polyval(cog,dlav0);
coga=polyfit(dlav(fi),log10(fvm(fi)),2);
gam0=10^polyval(coga,dlav0);
coen=polyfit(dlav(fi),log10(eneu(fi)),2);
en0=10^polyval(coen,dlav0);
 enmv(iso)=real(en0); 
 dlv(iso)=real(dlav0);
 gtv(iso)=real(gth0);
 gamv(iso)=real(gam0);
 dlv0(iso)=(isol-1)*2*pi;
end %if
end %length(imi)
end

if ~exist('gtv')
 ierrla=1;
' fine loop gamma, problemi ', keyboard 
 str0=str0*5;
 'set str', keyboard
 guessf_mod 
else
 ierrla=0;
end 

gtvm=gtv;
dlvm=dlv;
gammv=gamv;
enmmv=enmv;

[GT0e,igm]=min(real(gtve));
[GT0m,igm]=min(real(gtvm));
if GT0e<GT0m
 L0e=real(dlv(igm));
 Ga0e=real(gamev(igm));
 En0e=real(enmev(igm));
 [du,filam]=min(abs(dlvm-L0e));
 GT0m=gtvm(filam);
 L0m=dlvm(filam);
 Ga0m=gammv(filam);
 En0m=enmmv(filam);
else
 L0m=dlvm(igm);
 Ga0m=gammv(igm);
 En0m=enmmv(igm);
 [du,filam]=min(abs(dlve-L0m));
 GT0e=gtve(filam);
 L0e=dlve(filam);
 Ga0e=gamev(filam);
 En0e=enmev(filam);
end

%[GT0e,igm]=min(real(gtv));
%L0e=real(dlv(igm));
%Ga0e=real(gamv(igm));

if iplot==1
figure, plot(dlav,f11,dlav,f22), title(' fase TM' ), pausak
figure, plot(dlav,abs(s11),dlav,abs(s22)),
title(' modulo TM'), pausak
figure, plot(dlav,fa_tot0,dlv,dlv0,'ro'), title(' cond. fase  TM'), pausak
figure, semilogy(dlav,fvm,dlv,gamv,'wo',L0m,Ga0m,'gx'), pausak
figure, semilogy(dlav,eneu,dlvm,enmv,'wo',L0m,En0m,'gx'), pausak
figure, semilogy(dlav,gth,dlv,gtv,'ro',L0m,GT0m,'gx'),
hold on, semilogy(dlav,gthe,'--',dlve,gtve,'mo',L0e,GT0e,'cx'), title(' cond modulo TM (cont) -TE (- - - -)'), pausak
figure, semilogy(lav,gth,lambda_cen+dlv,gtv,'ro'),
a=axis;
%'key', keyboard
a(1)=min(lav);
a(2)=max(lav);
a(3)=min(gth);
a(4)=min(abs(gth))*100;
axis(a)
grid, pausak

end



  Gthv=real(GT0m)
  G0m=Gthv;
  la_ver=lambdau+L0m;
  GTM=Gthv;
  LTM=real(la_ver);
  nimm=GTM/(4e4*pi)*LTM;
  limm=L0m;

  GTM=GT0m;
  LTM=lambdau+L0m;  
  nimm=GTM/(4e4*pi)*LTM;
  limm=L0m;
  
  GTE=GT0e;
  LTE=lambdau+L0e;  
  nime=GTE/(4e4*pi)*LTE;
  lime=L0e;

end  % dret

n_i=n_sav;
%' fine guessmot', keyboard