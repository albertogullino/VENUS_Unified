
global icomp_grat 

%if ifp==-10
% iplot=1;
%else 
% iplot=0;
%end 
%iplo=iplot;

%'chain', keyboard


nl=101;
if exist('div_lam')
 str=lambda_cen/div_lam;  
else
 str=lambda_cen/15;  
end
str0=str;
if exist('nl0')
 nl=nl0;
else
 nl=101;
end

iproga=0;
iff=0;
freq=0;
lambdau=lambda_cen;
k02=4e4*pi/lambdau;


dlav=linspace(-str,str,nl);
lav=lambdau+dlav;

%lav=lambdau;

if exist('ibast')
else
 ibast=-2;
 par_grat=0;
end
global fsto L_i n_i rr rfd rfu  iLP ifp iff Lam0 ifun


n_i=nto;
if dret>0 & imet==1
n_i(fiirv)=nreticolo;
end
L_i=dto;
lambda_cenmod=lambda_cen;
%'fermo trans', keyboard
freq=0;
ierrla=0;

%guessf
'guess_mod !!!!!!!!! '
guessf_mod
%' dopo guessf_mod', keyboard


nimes=nime;
limes=lime;

if itetmt>=2
nimms=nimm;
limms=limm;
end

%'itetm', keyboard
global itetmt

itetmt_sav=itetmt;

if itetmt==3
 itetmtv=[1 2];
else
 itetmtv=itetmt;
end
%'chain', keyboard

for itetmt=itetmtv


 if itetmt==2
  nime=nimms;
  lime=limms;
 else 
  nime=nimes;
  lime=limes;
 end
%n_i=n_i_sav;

%z0r=



 if itetmt==1
  gav=gamev;
  gtv=gtve;
  dlv=dlve;
  env=enmev;
 else
  gav=gammv;
  gtv=gtvm;
  dlv=dlvm;
  env=enmmv; 
 end 
  
% GtQW=real(GT0/(NQW*Ga0)/2);
% nim=la_ver*GtQW/(4*pi)*1e-4; 
% z0=0+j*nim;
% z0=0;
% icomp_grat=1;
% [Ez,ztot,indz]=f_cam(z0);
% iff=1;
%[gazk,enan,fak,Tu,Tb,Gue,Gum,Gbe,Gbm,Lf,Lcav,ztot,Ez,Hz,indz,nmean,Perd,Ge,Gm,lqw]=...
%         th_scattu(fiQ,fiCav,L_i,n_i,rr,rfd,rfu,lambda,freq,0,iLP,ifp,iff,ibast,par_grat,NPZ);
clear GtQW GtQWo

[gad,iso]=sort(gtv);
dlv_sav=dlv;
gtv_sav=gtv;

if length(iso)>1
iscelta=iso(1:2);
else
iscelta=1;
end
%'gav', keyboard
gtv=gtv(iscelta);
dlv=dlv(iscelta);
for lmodo=1:length(dlv)
 lmodoi=lmodo;
 Lam0=lambda_cenmod+dlv(lmodo);
% la_ver=lambdau+L0e;
 la_ver=Lam0;
 lambda=la_ver;
 icomp_grat=1;
 if izetrasm==0
  NPZ=45;
  iff=1;
  [gazk,enan,fak,Tue,Tb,Gue,Gum,Gbe,Gbm,Lf,Lcav,ztot,Ez,Hz,indz,nmean,Perd,Ge,Gm,lqw,Ezm,KKie,KKim]=...
         th_scattu3(fiQ,fiCav,L_i,n_i,rr,rfd,rfu,lambda,freq,0,iLP,ifp,iff,ibast,par_grat,NPZ);    

  Em=abs(Ez').^2;
  lacav=real(lambda/rqw)/4;
  zce=mean(ztot(lqw));
  fice1=find(ztot>zce-lacav*1.5 & ztot<zce);
  [du,fii]=min(Em(fice1));
  zin=ztot(fice1(fii));
  fice2=find(ztot>zin*1.01 & ztot<zin+lacav*3);
  fice=[fice1; fice2];
  [du,fii]=min(Em(fice1));
  [du,fiu]=min(Em(fice2));
  pe=[fice1(fii) fice2(fiu)];
 
  if iplot==1
   figure, plot(ztot,Em*3,'.',ztot,abs(Hz).^2*3,ztot,real(indz),ztot(fice),Em(fice)*3,'r.',ztot(pe),Em(pe)*3,'wo')
   pausak
  end
 
  pev=pe(1):pe(2);
  pevd=pe(1):pe(2)+1;
  lqwd=[lqw lqw(end)+1];
  dnm=d*1e6;
  ze=diff(ztot(pev([1 end])));
  fizcav=find(ztot>z_cav(1)/1000 & ztot<z_cav(end)/1000);
  zc=ztot(fizcav);
  dzc=diff(zc);
  dzc=[dzc; dzc(end)];
  dzc=dzc/sum(dzc);
  mEqw=sum(dzc.*Em(fizcav));  
  rmed=real(sum(dzc.*indz(fizcav)'));   
  rmed=real(sum(dzc.*indz(fizcav)'.*Em(fizcav)))/mEqw;   
  epsi=abs(indz.^2).';
  eat=sum(Em(lqw).*epsi(lqw).*diff(ztot(lqwd)));
  eaca=sum(Em(pev).*epsi(pev).*diff(ztot(pevd)));
  fat=eat/eaca;
  rd=dnm/ze;
  enan=fat/rd;
 
  GtQWo(lmodoi)=gtv(lmodo)/(NQW*gav(lmodo))/enan;
  GtQW(lmodoi)=gtv(lmodo)/(NQW*env(lmodo));
  Fidud{lmodoi}=real(Em/mEqw*rmed)/2;
  zidud{lmodoi}=ztot;
  nidud{lmodoi}=indz;
 else
  fimag=find(gtv>0);
  if gtv(lmodo)<5*min(gtv(fimag))
   GtQW(lmodoi)=3000;
  else
     GtQW(lmodoi)=1e9;
  end
 end
  if izetrasm==2  & GtQW(lmodoi)<1e6              %metodo risonanza trasv
   iff=0;  %per calcolare il campo in z
   imap=1;
   lime=dlv(lmodo);
   gth=GtQW(lmodo);
   Lam0=lambda_cenmod+lime;
   nime=Lam0*gth/(4*pi)*1e-4; 


%  ' stop prima di nuova sub ', keyboard
   fzero_comp
   
   [fzever,fzmver]=f_mulut(z0,itetmt);
   if itetmt==1
    ver_zer=abs(fzever)
   else
    ver_zer=abs(fzmver)
   end
   if ifp==-10
    'verifica zero' 
    keyboard   
   end 
   if ver_zer>10000
    'entro fzero_raff', keyboard
    iscan_sav=iscan;
    iscan=1;
    fzero_raff
    iscan=iscan_sav;
   end

   

   if iplo==1
   %' stop prima di campo ', keyboard
   icomp_grat=1;
   [Ez,ztot,indz]=f_cam(z0,itetmt);
   Em=abs(Ez(1,:).^2);
   fizcav=find(ztot>z_cav(1)/1000 & ztot<z_cav(end)/1000);
   zc=ztot(fizcav);
   dzc=diff(zc);
   dzc=[dzc; dzc(end)];
   dzc=dzc/sum(dzc);
   mEqw=sum(dzc.*Em(fizcav)');  
   rmed=real(sum(dzc.*indz(fizcav)'));   
   rmed=real(sum(dzc.*indz(fizcav)'.*Em(fizcav)'))/mEqw;   
   Fidu=real(Em/mEqw*rmed)/2;
   else
   Em=0;
   ztot=0;
   Fidu=0;
   Fipi=0;
   end
   la_ver=Lam0+real(z0);
   k0=2*pi/la_ver;
   gth=2*k0*imag(z0)*1e4
   if gth<0
    iscan_sav=iscan;
    iscan=1;
      fzero_comp
      la_ver=Lam0+real(z0);
      k0=2*pi/la_ver;
      gth=2*k0*imag(z0)*1e4      
     iscan=iscan_sav; 
%    gth=1e6;
%    Fidu=0;
%    ztot=0;
%    indz=0;
   end 
     if ver_zer>1
      GtQW(lmodo)=1e9;
     else 
      GtQW(lmodo)=gth;
     end   
    Fidud{lmodoi}=Fidu;
    zidud{lmodoi}=ztot;
    nidud{lmodoi}=indz;
%    ' ver lmodo ', keyboard
  end  % fine metodo determinazione modo
  
  
 end %lmodo
 ga_spet=1;
% ga_spet=1-(dlv/.08).^2;
 
  [dgth,is]=min(GtQW./ga_spet);
  gth=GtQW(is);
%  ' GtQW ver', keyboard
  la_ver=lambda_cenmod+dlv(is);
  lime=dlv(is);
  Fidu=Fidud{is};
  Lvdu=la_ver;
  Gthdu=gth;
  indz=nidud{is};
  ztot=zidud{is};
  if ifp==-10
   ' fine controllo lmodo', keyboard
  end
%  ' fine controllo lmodo', keyboard  
if izetrasm==1                %metodo risonanza trasv
 iff=0;  %per calcolare il campo in z
 imap=1;
 nime=la_ver*gth/(4*pi)*1e-4; 
 Lam0=lambda_cenmod+lime;

%' stop prima di nuova sub ', keyboard
 fzero_comp

 if iplo==1
 %' stop prima di campo ', keyboard
 icomp_grat=1;
 [Ez,ztot,indz]=f_cam(z0);
 Em=abs(Ez(1,:).^2);
 fizcav=find(ztot>z_cav(1)/1000 & ztot<z_cav(end)/1000);
 zc=ztot(fizcav);
 dzc=diff(zc);
 dzc=[dzc; dzc(end)];
 dzc=dzc/sum(dzc);
 mEqw=sum(dzc.*Em(fizcav)');  
 rmed=real(sum(dzc.*indz(fizcav)'));   
 rmed=real(sum(dzc.*indz(fizcav)'.*Em(fizcav)'))/mEqw;   
 Fidu=real(Em/mEqw*rmed)/2;
 else
 Em=0;
 ztot=0;
 Fidu=0;
 Fipi=0;
 end
 la_ver=Lam0+real(z0);
 k0=2*pi/la_ver;
 gth=2*k0*imag(z0)*1e4
end  % fine metodo determinazione modo



  Gthdu=gth;
  Lvdu=real(la_ver);
  

if ifp==-10 & iplo==1
 figure, semilogy(ztot,Fidu,ztot,real(indz)), 
 title([' lambda_{res} = ',num2str(Lvdu),' Gth = ',num2str(gth)]), pausak
end  

  
if itetmt==1 
 ztote=ztot;
 Fipi=Fidu;
 Lvpi=Lvdu;
 FiTE=Fipi;
 Gthp=Gthdu;
 GTE=Gthp;
 LTE=Lvpi;
% 'assegno valori TE', keyboard
end


if itetmt==2 
 ztotm=ztot;
 Fivi=Fidu;
 Lvvi=Lvdu;
 FiTM=Fivi;
 Gthv=Gthdu;
 GTM=Gthv;
 LTM=Lvvi;
% 'assegno valori TM', keyboard

end  


end  %itetmtv

itetmt=itetmt_sav;
