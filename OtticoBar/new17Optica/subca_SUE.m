
global icomp_grat 
if ifp==-10
 iplot=1;
 iplot=ifp;
else 
 iplot=0;
end 
iplo=1;
%iplo=0;



nl=31;

str0=str;
%str=lambda_cen/10;   
% 'SUE dentro', keyboard

iproga=0;
iff=0;
freq=0;
lambdau=lambda_cen;
k02=4e4*pi/lambdau;


dlav=linspace(-str,str,nl);
lav=lambdau+dlav;

%lav=lambdau;
%' ibast', keyboard
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
%'guess_mod !!!!!!!!! ', keyboard
guessf_mod

%' dopo guess ', keyboard

if iraffina==0
 Fipi=0;
 Gthp=GteQW;
 Lvpi=lambda_cen+lime;
 ztote=0;
% 'return no raffina ', keyboard
 return
end
 %'raffina ', keyboard

nimes=nime;
limes=lime;

global itetmt

itetmt_sav=itetmt;

if itetmt==3
 itetmtv=[1 2];
 nimms=nimm;
 limms=limm; 
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
%' prima di ricerca', keyboard

clear GtQW GtQWo GtQWb dlvbuo
lmodobuo=0;

[du,imi]=min(gtv);
if length(gtv)>=400
lsce=imi+[-1 0 1];
 if lsce(end)>length(gtv)
  lsce=lsce-1;
 end
 if lsce(1)<=0
  lsce=lsce+1-lsce(1);
 end 
else
lsce=imi;
end

%for lmodo=1:length(gtv)
for lmodo=lsce
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
  Em=abs(Ez');
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
 
  GtQWo(lmodo)=gtv(lmodo)/(NQW*gav(lmodo))/enan;
  GtQW(lmodo)=gtv(lmodo)/(NQW*env(lmodo));
  Fidud{lmodo}=real(Em/mEqw*rmed)/2;
  zidud{lmodo}=ztot;
  nidud{lmodo}=indz;
 else
  fimag=find(gtv>0);
%  if gtv(lmodo)<2*min(gtv(fimag))
   GtQW(lmodo)=3000;
%  else
%     GtQW(lmodo)=1e9;
%  end
 end

  if izetrasm==2  & GtQW(lmodo)<1e6              %metodo risonanza trasv
%    'calcolo', keyboard
   iff=0;  %per calcolare il campo in z
   imap=1;
   lime=dlv(lmodo);
   gth=GtQW(lmodo);
   Lam0=lambda_cenmod+lime;
   nime=Lam0*gth/(4*pi)*1e-4; 


%  ' stop prima di nuova sub ', keyboard
   fzero_comp
   
   [fzever,fzmver]=f_mulut(z0);
   if itetmt==1
    ver_zer=abs(fzever)
   else
    ver_zer=abs(fzmver)
   end
   if ifp==-10
%    'verifica zero' ,    keyboard   
   end 
   if ver_zer/fp_conf>.02
    'entro fzero_raff'
    iscan_sav=iscan;
    iscan=1;
    fzero_raff
    iscan=iscan_sav;
   end
   
%' qui fie', keyboard   

   if iplo==1
%   ' stop prima di campo ', keyboard
   icomp_grat=1;
   [Ez,ztot,indz]=f_cam(z0);
   Em=abs(Ez(1,:).^2);
   Em=abs(Ez(1,:));
   if iLOC==0
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
      Fidu1=real(Em);
      Fidu=Fidu1/max(Fidu1)*rr;
   end
   else
   Em=0;
   ztot=0;
   Fidu=0;
   Fipi=0;
   end
   la_ver=Lam0+real(z0);
   k0=2*pi/la_ver;
   gth=2*k0*imag(z0)*1e4
   
   if gth>0
    lmodobuo=lmodobuo+1;
    dlvbuo(lmodobuo)=dlv(lmodo);
    GtQWb(lmodobuo)=gth;
    Fidud{lmodobuo}=Fidu;
    zidud{lmodobuo}=ztot;
    nidud{lmodobuo}=indz;
   end 
  end  % fine metodo determinazione modo
  
  
 end %lmodo
 
if ~exist('dlvbuo')
    lmodobuo=lmodobuo+1;
    dlvbuo(lmodobuo)=dlv(lmodo);
    GtQWb(lmodobuo)=GteQW;
    Fidud{lmodobuo}=Fidu;
    zidud{lmodobuo}=ztot;
    nidud{lmodobuo}=indz;
end

 if exist('dlvbuo')
  dlv=dlvbuo;
  ga_spet=1-(dlv/.08).^2;
  [dgth,is]=min(GtQWb./ga_spet);
  gth=GtQWb(is);
  lime=dlv(is);

  lime=dlv(is);
  Fidu=Fidud{is};
  Lvdu=la_ver;
  Gthdu=gth;
  indz=nidud{is};
  ztot=zidud{is};
 else 
  [gth,imi]=min(gtv);
  lime=dlv(imi);  
 end
  la_ver=lambda_cenmod+lime;
 if ifp==-10
%   ' fine controllo lmodo', keyboard
  end
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
% Em=abs(Ez(1,:));
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
  

%if ifp==-10 & iplo==1
% figure, semilogy(ztot,Fidu,ztot,real(indz)), 
% title([' lambda_{res} = ',num2str(Lvdu),' Gth = ',num2str(gth)]), pausak
%end  

  
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
global la1Dr ga1Dr
la1Dr=lambda_cenmod+dlvbuo;
ga1Dr=GtQWb;
%' fine hcgw ', keyboard