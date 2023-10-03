STZ0=.1;
STZ0=1;
iLOC=1;
iwo=1;
%NP_dec=5;
NP_dec=1;
%NP_dec=15;
%ifceck_zero=-10;   

MOSTRO_SOLUZIONE=0;
MOSTRO_SOLUZIONE=1;

if MOSTRO_SOLUZIONE==0
 ifceck_zero=0;  % plotta zero
 icheck_guess=0;   % 1 plotta guess
else
 ifceck_zero=1;  %1  plotta zero
 icheck_guess=1;   % 1 plotta guess
end

iplotta=0;
iraffina=iref;

%'sol0_sub ', keyboard

if ifp==-10
 iplotta=1;
end 

fi_znofield

izetrasm=2;
iscan=0;
ifollow=0;
irepiano=0;
iben=1;  %attivo in res_anis
ibasto_sing=1;    % applico Bastonero ai singoli strati, caricandolo di mano in mano
ibasto_sing=0;    % BW  ai singoli strati,
iret_sep=0;   %=1 tratta i 2 reticoli separati, 0 mette tutto insieme, anche strati piani.
icarico=1;   % calcola T reticolo caricato sui G_i

ires_scatt=0;        % risolve con risonanza trasversale
ineff=1;   %0: BW, 1:Bastonero isolato, 2: Bast caricato
irefine=0;   % ricalcola n_eff perche varia lambda

itetm=-1;
nvero=1;



%isub_old=1;   %aggiusta guadagno su una solo standing wave
isub_old=1;   %aggiusta guadagno su una solo standing wave
%isub_old=0;   %metodo L_ef (chiama sub_transunuu in sub_transu)



fish=find(shavet==6);
fishu=find(shu==6);

if length(fish)>0  
%' iBast', keyboard
  par=radii.array{fish};
  igr_NEW=0;
  
  if length(par)>=11
   par_gratNEW=par{11};
    if isstruct(par_gratNEW)==1
     igr_NEW=1;
    else
%     'errore reticolo', keyboard
     'reticolo VECCHIO', keyboard
    end
   end 
  NModi=11; 
  if igr_NEW==0 
   ngrating=nv(fish,:);
   par_grat.n1=ngrating(1);
   par_grat.n2=ngrating(2);
   n1g=ngrating(1);
   n2g=ngrating(2);
   period=par{5};
   t1=par{6};
   t2=period-t1;
   DC=t1/period;
 %  par_grat.itetm=3;
   par_grat.px=period;
   par_grat.DC=DC;
   par_grat.NModi=NModi;
%  itetm=3;
   icarico=0;   % calcola T reticolo caricato sui G_i
   dret=L_i(fishu);
  else
   par_grat.n1=par_gratNEW.nr1;
   par_grat.n2=par_gratNEW.nr2;
   par_grat.r1=par_gratNEW.nr1;
   par_grat.r2=par_gratNEW.nr2;
   par_grat.px=par_gratNEW.Pe;
   par_grat.per=par_gratNEW.Pe;
   par_grat.DC=par_gratNEW.dcv(1);
   par_grat.th=par_gratNEW.th;
   par_grat.NModi=NModi;
 
  
  %'spn qwio', keyboard
  end
   
  %' cont grat', keyboard
  ibast=fishu+1;
else 
 ibast=[];
 dret=0;
 par_grat=0;
end  

imet=2;
ioldeig=0;
global isempl ibast  iorta icrescita icarico
icrescita='top';
step=1;
if imet==2
iorta=1;   % in th_scattu, f_cam usa orta
%iorta=0;   % in th_scattu, f_cam usa orta
%iorta=input(' iorta = ');   % in th_scattu, f_cam usa orta
end
isempl=1;  
iord_long=[1];
% STZ -> 
% step discretizzazione in z: 
%  dJ -> hz
%------------------------------------------------%


if ioldeig==1
 iord_long=[1];
end

lambda_cen=lambda;  

str=lambda/150;  
str=lambda/15;  
str=lambda/20;  
str0=str;  
rfu=n1(1);
rfd=n1(end);
dJ=cumsum([0; dto]);
nJ_ve=nJ;
nJ_pa=nJ;
lJ=length(dJ);


if imet==0
 lambda0=lambda0p;
 STZ=STZ1;
 set_zeta_mult
 fizcav=find(hz>z_cav(1) & hz<z_cav(end));
 zc=hz(fizcav);
 dzc=diff(zc);
 dzc=[dzc dzc(end)];
 dzc=dzc/sum(dzc);

% mEqw=sum(dzc.*Em(fizcav))*2;    
 for kmm=1:length(lambda0) 
  lambda0i=lambda0(kmm); %lambda0i: singola lambda guess 
  [Ksip,lambdap,Fipi,uLongp,uLong0p,Fap]=eiglmio(lambda0i,uFunc,uF0,rel_pa,STZ);
  uLp(1)=uLong0p;
  uLp(2)=uLongp;
 
  fatqwp=uLongp/(uLong0p*nmir.a);
%  manp=max(uF0.*Fipi.');
  manp=sum(dzc.*Fipi(fizcav)')*2;    
  rmed=sum(dzc.*sqrt(real(rel_pa(fizcav))).*Fipi(fizcav)')*2/manp; 
  rm_pa=rmed;
  Gamp=uLongp;
  Gthp = 2e-2*pi*rr/lambdap*imag(Ksip)/Gamp;
  Fipi=Fipi/manp*rmed;
  Lvpi=lambdap*1e6;
  ztot=hz;
 end
else
  %ver_sub
  if imet==1
   itetmt=1;
  else
   itetmt=3;
   if dret==0
    itetmt=1;
   end
  end
  nreticolo=1.5;

%  'prima sub_transu 0', keyboard  
  %fiQdu=fiQ;
%  par_grat.itetm=itetmt;
  fiQ=find(fsto==-1);
  fiQW=fiQ(ceil(end/2));
  %fiQ=fiQdu;

global NPlam dlam

if length(NPlam)==0
  dlam=lambda*.03;
  NPlam=71;
end  

%' NP lam E2', keyboard  
  
  kt=0.12;
  kt=0.;
  global Options
  Options.iAlb=0;
  if ifp==-10
%   'prima di Caricato', keyboard  
  end
%   keyboard
%   'prima di Caricato', keyboard  
global E2 ztot nz

if length(E2)==1
 ifpSA=ifp;
 if Ps.ifpstop==1
 ifp=-10;
 end
 
 
  [gam,lam,ztot,Ez,nz,Gaqw,NQW_ef,az,uL,E_Temp,z_Temp,velm1D]=fiezCMMCaricato(lambda,L_i,n_i,fiQW,fiCav,fiQ,ifp,dlam,NPlam,rr,kt,ibast,par_grat,ato);
  ifp=ifpSA;
else

 Ez=E2;

end  
%ga,la,zet,Ez,nz,Gaqw,NQW_ef,az,uL
%  sub_transu
if iany==1
% global EO ABS
% close all
% save abs
' anisotropia', keyboard
 anisotropiaEO
end
if ifp==-100
' Losses', keyboard
end


global EO ABS

if ABS.e+ABS.h>0
 Sub_Losses_zeta
 nitn=nr.t.'; 
 nib=nr.b.'; 
 niat=nr.a.'; 
 nitot=[nitn.'; niat.'; nib.']; 
% 'dp nt\', keyboard
end


  ztote=ztot*1000;  
  if imet==2 & itetmt==30
   ztotm=ztotm*1000;
   Fivi=FiTM;
   Gthv=GTM;
   Lvvi=LTM;   
  end 
end
  Fip0=abs(Ez);
  Fipi=Fip0;
  Gvp0=gam;
  Lvp0=lam;

if imet==0
 zetae=hz;
else 
 zetae=ztote;
end
imetplo=imet;
imetplo=0;
lambda=Lvp0;
gth=Gvp0;

if ~exist('St_wave')
 St_wave(:,1)=ztot';
 St_wave(:,2)=Ez';
end

return

if ks==1 
%'fermo per cont', keyboard
iraffz=1;
[Ksi,lambdas,Fi,uLong,uLong0,Fa]=eiglmiu(Lvp0*1e-6,uFunc,uF0,relPerm,hn,d,iraffz);
%'fermo per cont', keyboard
 uLv(1)=uLong0;
 uLv(2)=uLong;
%'fermo per cont', keyboard
% uuF
uL=uLv;
end

%'fermo per cont iplotta', keyboard
if iplotta==1 & ifp==-10
 h=figure, 
 set(h,'position',[ 129    75   667   645])
 
 if itetm==1
  subplot(211)
  if iwo==0
   if imetplo<2
    plot(hz,real(perm_pa),'r',zetae,Fip0.^2),
  else
    plot(hz,real(perm_pa),'r',zetae,Fip0.^2,'.'),
  end
 else
      plot(ddplot,real(ndplot),'r',zetae,(Fipi/max(Fipi)).^2*max(ndplot)),
 end 
 ax=axis;
 ax(4)=5;
 axis(ax)
 title(['TMM:   lambda_{res} = ',num2str(Lvp0),' Gth pa= ',num2str(Gvp0)]), 
 

 fatqw=uLong/(uLong0*nmir.a);
 lambda_num=lambdas*1e6;

 %'fermo lambdas', keyboard
 gpla=2e4*pi*rqw/lambda*imag(Ksi)/uLong0/fatqw; 
 gpla=2e4*pi*rqw/lambda*imag(Ksi)/uLong; 
 manv=max(uF0.*Fi.');
 confzv=uL;
%'ulong ', keyboard
 dorigin=0;
 x=(hz'-dorigin)*1e-3;
 rperm=real(sqrt(relPerm));
 faca=3;
subplot(212)
 
  plot(x*1000,rperm,'r',zetae,(Fipi/max(Fipi)).^2*faca,x*1000,faca*Fi/max(Fi),'g')
%    title([' lambda_{res} = ',num2str(lambda),'  Ksi = ',num2str(Ksi)]), pausak
    title([' Numerical: lambda_{res} = ',num2str(lambda_num),'  Gth = ',num2str(gpla), '  --- GIALLO TMM']), 
    keyboard
 man=max(uF0.*Fi.')/3.5; 
 
else  %itetm

 subplot(211)
 faca=3;
  if iwo==0
   if imetplo<2
    plot(hz,real(perm_pa),'r',zetae,Fip0.^2),
  else
    plot(hz,real(perm_pa),'r',zetae,Fip0.^2,'.'),
  end
 else
      plot(ddplot,real(ndplot),'r',zetae,(Fipi/max(Fipi)).^2*faca),
 end 
 ax=axis;
 ax(4)=5;
 axis(ax)
 title(['TMM:  TE  lambda_{res} = ',num2str(Lvp0),' Gth pa= ',num2str(Gvp0)]), 
 

 fatqw=uLong/(uLong0*nmir.a);
 lambda_num=lambdas*1e6;

 %'fermo lambdas', keyboard
 gpla=2e4*pi*rqw/lambda*imag(Ksi)/uLong0/fatqw; 
 gpla=2e4*pi*rqw/lambda*imag(Ksi)/uLong; 
 manv=max(uF0.*Fi.');
 confzv=uL;
%'ulong ', keyboard
 dorigin=0;
 x=(hz'-dorigin)*1e-3;
 rperm=real(sqrt(relPerm));

subplot(212)
 
  plot(x*1000,rperm,'r',ztotm,(FiTM/max(FiTM)).^2*faca,'g')
%    title([' lambda_{res} = ',num2str(lambda),'  Ksi = ',num2str(Ksi)]), pausak
 title(['TMM:   TM  lambda_{res} = ',num2str(LTM),' Gth pa= ',num2str(GTM)]), 
    keyboard
 man=max(uF0.*Fi.')/3.5; 


end

end

lambda=Lvp0;
gth=Gvp0;

if imet==2 & dret>0
 lambda=LTE;
 gth=GTE; 
 if GTE>GTM
  lambda=LTM;
  gth=GTM; 
  ' attenzione ai GAMMA confinamento con reticolo hcg!!!!!!!'  , %keyboard
 end
end

 %'fermo lambdas', keyboard

%' fine ', keyboard