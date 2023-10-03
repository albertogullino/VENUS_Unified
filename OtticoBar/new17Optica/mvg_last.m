%' entro mvg_last', keyboard


if ~exist('sT')==1
 sT=0;
end
if ~exist('sN')==1
 sN=0;
end
if ~exist('ivacc')==1
 ivacc=0;
end
if ~exist('ilossk')==1
 ilossk=[0 0];
end
xroJ=xroI;
deltany=0;
imod_err=0;
%'mvt', keyboard

immax=length(mmvet);
global igamveb igamveu GGbext GGuext

 if exist('igamveb')==0 | length(igamveb)==0
  igamveb=1;
 end
 if exist('igamveu')==0 | length(igamveu)==0
  igamveu=1;
 end


 if exist('iplan')==0 | length(iplan)==0
  iplan=0;
 end

 if exist('ifnm')==0 | length(ifnm)==0
  ifnm=0;
 end

 if exist('isopro')==0
  isopro=0;
 end
 %' set_stru ', pausak
 set_stru

if i1D==1
  return
end
%'setstru', keyboard


 nk1max=sum(nK_dis);
 nk1=nK_dis(1);


 if length(alim_in)>1
   alimi=alim_in(1);
   alim=alim_in(2);
   if length(alim_in)==3
    alimu=alim_in(3);
    if length(nK_dis)==1
     'nK_dis must have 2 elements'
     keyboard
    else
     nk2=nK_dis(2);
    end
   end

 else
  alimi=0;
  if alim_in==0
   fiam=find(aitot~=0);
   a_min=min(aitot(fiam));
   afit=[2.1  3    4   6  10  20 ];
   kfit=[.25  .22 .15 .1 .05 .03];
   [du,fiam]=min(abs(a_min-afit));
   alim=kfit(fiam);
  else
   alim=alim_in;
  end
 end


%'mvg Dlam', keyboard


 Dlam_mo=Dlam_mod(1:2);
 if length(Dlam_mod)==2
  Ndisp=3;
 else
  Ndisp=Dlam_mod(3);
 end

 Dso=sort(Dlam_mo);
 Frisi=0;
 Frisi=Dso(1)*1e-3/lambda;
 Frisu=Dso(2)*1e-3/lambda;
 if length(Dlam_mod)>=4
  DelFr=Dlam_mod(4)*1e-3/lambda;
 else
  DelFr=0;
 end
 if isnan(DelFr)==1
 DelFr=0;
 end

%'contr Dle', keyboard

nv_sa=nv;

 avero=a0ref;
 %ra=niat(1);

 k0=2*pi/lambda;

 kcav=2*pi/lambda*rr;
 kcav0=kcav;
 kcav00=kcav*1e6;
 a=kcav0*avero;     %adimensionale:normalizzato rispetto a k0 della cavita'

%' kcav', keyboard
 igint=1;   % coupl. coeff. integration
 igint=0;   % coupl. coeff. integration



         if igau==0 | igau==1
          npx1=100;
          npx2=100;
          xi0=1;
          yi0=1;
          lxi0=1;
          nst=1;
          yiv=1;
         end

          npx1=50;
          npx2=20;
        if sT>=2
         yiT=[1:100]; %solo per renderlo un vettore
         if sT==2
          yiT=[];
          if igau>=4 & length(T)>1
           xiT1=rodisT(2:length(rodisT));
           xiT2=rodisT(1:length(rodisT)-1);
           xiT=(xiT1+xiT2)/2;
%           xiT=[(xiT1+xiT2)/2; rodisT(end)];
           s=size(Tdis);
           for iz=1:s(2)
           yiT(:,iz)=-diff(Tdis(:,iz));
%            yiT(:,iz)=[-diff(Tdis(:,iz)); Tdis(end,iz)];
           end
          end
%          zdis=(zdis(1:end-1)+zdis(2:end))/2;
          zdis=zdis(1:end-1);
          if ifp~=-4
           'in mvg_last '
           hT=figure;
           set(hT,'pos',[98   423   1200  525])
           subplot(121)
           for ky=1:length(zdis)
            yiTric(:,ky)=flipud(cumsum(flipud(yiT(:,ky)),1))+Tdis(end,ky);
           end 
           plot(zedisp,yiTric,'.'), 
           title('Temp ricostruita')
           subplot(122)
%           xiTp=[xiT; ro_in(end)];
%           plot(xiT,Tdis,'.'), pausak
           plot(zdis,Tdis,'.'), 
           title('Temp originale')           
           pausak           
          end
         end
        end
        if sN==2
         yiN=[1:100]; %solo per renderlo un vettore
        end

iLPr=iLP;
imm=0;
imod=0;
itetmv=ivacc;
if ipolar==0
 ivacc=-1;
end
%' cont itetmv mvg_last', keyboard
global nucompv nucomp
nmasces=nmasce;
%' cont nucont', keyboard

if length(nucompv)==0
  nucompv=nucomp;
end 
sovet=[3 5];
sovet=mmvet*2+3;
nucomp=nucompv;
icomo=0;
numodiacc_sav=numodiacc;
for mm=mmvet
clear Trcrit
icomo=icomo+1;
so=sovet(icomo);
mmsav=mm;
if length(nmasces)>1
 if mm==1
  nmasce=nmasces(1);
 else
  nmasce=nmasces(2);
 end
end
if length(nucompv)>1
 if mm==1
  nucomp=nucompv(1);
 else
  nucomp=nucompv(2);
 end
end
%'kazzo', keyboard
global igr_app
if length(numodiacc_sav)>1
 numodiacc=numodiacc_sav(icomo);
end
ipolar_sav=ipolar;
 if iLP==0 & mm==0 & igr_app==1 & numodiacc==0
  numodiacc=1
  'FORZO numodiacc !!!!!!!! ', pausak
%  ipolar=1;
%  ipolar=2;
 end

 
 if iLP==10 & mm==0 & ivacc(1)~=-1
  if ivacc(1)==0
   itetmv=[1 2];
  else
   itetmv=ivacc;
  end
  iLP=1;
  iLPr=iLP;
  mm=mm+1;
  ipolar=1;
 end
 %' prima sosta', keyboard
 
 for itet=1:length(itetmv)
   itetm=itetmv(itet);
   iLP=iLPr;

 %' seconda sosta', keyboard
   Frisi=Dso(1)*1e-3/lambda;
   Frisu=Dso(2)*1e-3/lambda;
%   if (is_even(mm)==1 & iLP==0) | (mm==1 & iLP==1)
   if (mm==0 & iLP==0) 
    mmul=1;
   elseif (mm==1 & iLP==0) 
    mmul=0;    
   else
    mmul=mm;
   end 

    Frisi=Frisi+DelFr*mmul;
    Frisu=Frisu+DelFr*mmul;
%    if exist('Fsov') & isnan(Dlam_mod(4))==1
%     Frisi=min(Fsov);
%     Frisu=Frisi+diff(Dso)*1e-3/lambda;
%    end
%   end
   mm_ver=mm;
   imm=imm+1;
  clear KAp KA icousav
  icalcola=1;
  igainshape=0;

%' qui Fris', keyboard

     if imod_err==1
      return
     end
 if ~exist('lambdavet')
  lambdavet=lambda;    
 end
 for ilaord=1:length(lambdavet)
  ipolar=polvet(ilaord);
  uL=uLv(ilaord,:);
  igouadd=1;
  lambda=lambdavet(ilaord)
%  'carico altra lambda', keyboard
%  load sa  
%  if ilaord==1
%   dla_fre=0;
%  else
%   dla_fre=-(lambdavet(ilaord)/lambdavet(1)-1);
%  end
  
  k0=2*pi/lambda;
  kcav=2*pi/lambda*rr;
  kcav0=kcav;
  kcav00=kcav*1e6;
  a=kcav0*avero;      
%'prima sub', keyboard
     sub
     global lp1 isavetutto lp2
     if length(isavetutto)==0
      isavetutto=0;
     end
%   'dopo sub', keyboard
isavetutto=-10;
    if isavetutto>0
      noms=['rad',num2str(isavetutto),'_',num2str(ilaord),num2str(lp1),num2str(lp2)];
      if iLP==0
       save BMat MKoszmd MKoszpd  MKospd MKosmd  Koszp Koszm  Kosp Kosm Kos Kosz  ary ary_s Fi
       clear Tstort Tstortp  MKoszmd MKoszpd MKoszd  MKospd MKosmd  Koszp Koszm  Kosp Kosm Kos Kosz Tstor Tstorb var
       clear DR ary ary_s Fi Tmirro Tstof Tstorbp Tdum Tcm Tc Pow Oo Odu Mns Mis
      else
       save BMat  MKospd MKosmd  Kosp Kosm Kos  ary ary_s Fi
       clear Tstort Tstortp  MKospd MKosmd   Kosp Kosm Kos  Tstor Tstorb var
       clear DR ary ary_s Fi Tmirro Tstof Tstorbp Tdum Tcm Tc Pow Oo Odu Mns Mis
      end
      eval(['save ',noms])
      load BMat
      
     end     
   end    
lambda=lambdavet(1);
     
     numodiacc=numodiacc_sav;
     ipolar=ipolar_sav;
   Ppo{imm}=Ppol;
   Ppt{imm}=Plot;
%   'dopo sub 1', keyboard
%   'dopo sub', keyboard
 end % itetm
%   'dopo sub 2', keyboard

end   %mm
     numodiacc=numodiacc_sav;
% prendo solo i modi voluti da nmasce
%save modo

isalto=0;
if isalto==0
%' scleta modi', keyboard
%' scleta modi', keyboard
%' scleta modi', keyboard
Plots=Plot;
clear Plot;
%iso=find(Gsov>0);
% 'controllo Plot', keyboard
% 'controllo Plot', keyboard
 %'controllo Plot', keyboard


  inval=find(Gsov<=0);
%  Gdu=Gsov.*M2v;
  Gdu=Gsov;
  Fdu=Fsov;
  Gdu(inval)=1e100;
%  isoga=0;
  if isoga==1
   [du,idu]=sort(Gdu);
   %ival=find(du<1e100);
   %iso=idu(ival);
   iso=idu;
  else
   [du,iso]=sort(Fdu);
%   iso=find(Gdu<1e90);
  end
  %'1ui', keyboard
 io=0;
 if io==1
  InR=IndRad(iso);
  if nmasce<0
   if length(iso)>abs(ipolar)*abs(nmasce)
   isu=[];
    for kmod=mmvet
     fij=find(InR==kmod);
     isu=[isu iso(fij(1:abs(nmasce)))];
    end
    iso=isu;
   else
    isu=iso;
   end 
  end
 else
  niso=abs(ipolar)*abs(nmasce);
  if length(iso)>niso
   iso=iso(1:niso);
  end
 end
%  ' nmasce', keyboard
%  iso=isu;
kco=0;
for kis=iso                                             
  kco=kco+1;					
  Plot.Ap{kco}=Plots.Ap{kis};		
  Plot.ApQ{kco}=Plots.ApQ{kis};		
  Plot.parmod{kco}=Plots.parmod{kis};	
  Plot.XP{kco}=Plots.XP{kis};		
  Plot.YP{kco}=Plots.YP{kis};		
  Plot.X{kco}=Plots.X{kis};			
  Plot.Y{kco}=Plots.Y{kis};			
  Plot.E2xo{kco}=Plots.E2xo{kis};		
  Plot.E2xp{kco}=Plots.E2xp{kis};		
  %Plot.FF{kco}=Plots.FF{kis};			
  if iLP==0
   Plot.E2yo{kco}=Plots.E2yo{kis};		
   Plot.E2yp{kco}=Plots.E2yp{kis};		
   Plot.E2zo{kco}=Plots.E2zo{kis};		
   Plot.E2zp{kco}=Plots.E2zp{kis};		
  end 
  Plot.Ef{kco}=Plots.Ef{kis};			
  Plot.Efx{kco}=Plots.Efx{kis};			
  Plot.Efy{kco}=Plots.Efy{kis};			
  Plot.Cug{kco}=Plots.Cug{kis};		
  Plot.gou{kco}=Plots.gou{kis};		
  Plot.aou{kco}=Plots.aou{kis};		
  Plot.fou{kco}=Plots.fou{kis};		
  Plot.ze{kco}=Plots.ze{kis};			
  Plot.gg0{kco}=Plots.gg0{kis};		
%  Plot.FF{kco}    =Plots.FF{kis};
  Plot.parmod{kco}=Plots.parmod{kis};

 end 					

 
% 'controllo Plot', keyboard 
   ADc=           ADc(:,iso);			
   ADout=        ADout(:,iso);	
   M2v=	M2v(iso);
   if exist('Purspet')
    Purspet=	Purspet(iso);
   else
     Purspet=	1;
   end
   Gsov=	         Gsov(iso);
   Fsov=	          Fsov(iso);
   Nazim=	         Nazim(iso);
   Nrad=	         Nrad(iso);
%   tyE=		   tyE(iso);
   EDc.x=      EDc.x(:,:,iso);
   EDout.x=   EDout.x(:,:,iso);
   if iLP==0
    EDc.y=      EDc.y(:,:,iso);
    EDout.y=   EDout.y(:,:,iso);
   end      
%   'qui Lsov', keyboard
 if idyn>=1      
   Lsov=         Lsov(iso);
   Gam_v=	    Gam_v(iso);
   Wv=	      Wv(iso);
   czv=	       czv(iso);
   fqwv=	      fqwv(iso);
   fa4v=	       fa4v(iso);
   Lefv=	       Lefv(iso);
  Tefv=	      Tefv(iso);
  Tefiv=	      Tefiv(iso);
  taut=	      taut(iso);
  tauu=	      tauu(iso);
  taub=	      taub(iso);
  teff=	       teff(iso);
  alca=	      alca(iso);
  pvol=	      pvol(iso);
  pmet=	      pmet(iso);
end
     
end     
%if ifp==-10
%' qui mvg_last ', keyboard
%end

% iLP=iLP1;
if ~exist('Gsov')
 Eqw=[];
 Eout=[];
 xro=[];
 fian=[];
 lambda=[];
 delf=[];
 gain=[];
 ord=[];
 nrAzim=[];
 Cu=[];
 PaDy=[];
 ADom=[];
 ADcm=[];
 Plot=[];
 Ppol=[];
 Dla_out=[];
 imod_err=[];
 PPlot=[];
 return
end

%'qio ve', keyboard
tyPmo=tyE;
nrmo=nazim;
namo=nrad;
Wvmo=PD.a;
czmo=PD.b;
fqwmo=PD.c;
Lfmo=PD.d;
Tfmo=PD.e;
Teff_all=PD.eu;
Tpemo=PD.e1;
Tmemo=PD.e2;
Tfimo=PD.f;
gtmo=PD.g;
famo=PD.t;

gmo=gsov;
fmo=fsov;
%  fi=find(gmo>0);
%  if isoga==1
%  [gdu,ig]=sort(gmo(fi));
%  else
%    [gdu,ig]=sort(fmo(fi));
%  end
%  if nmasce<0
%   if length(ig)>abs(nmasce)
%    ig=ig(1:abs(nmasce));
%   end 
%  end  
if length(iso)>nmasce
 ig=iso(1:nmasce);
else
 ig=iso;
end 
fig=ig;  
%  fig=fi(ig)';
  gmod=gmo(fig)';
  fmod=fmo(fig)';
  tyPmod=tyPmo(fig);
  nrmod=nrmo(fig);
  namod=namo(fig);
  Wvmod=Wvmo(fig);
  famod=famo(fig);
  czmod=czmo(fig);
  fqwmod=fqwmo(fig);
  Lfmod=Lfmo(fig);
  Tfmod=Tfmo(fig);
  Teff=Teff_all(fig);
  Tpemod=Tpemo(fig);
  Tmemod=Tmemo(fig);
  Tfimod=Tfimo(fig);
  gtmod=gtmo(fig);

  nord=[namod; nrmod];
%  'ord', keyboard
%  ADom=ADo;
%  ADcm=ADc;
%   ADc=           ADc(:,iso);			
%   ADout=        ADout(:,iso);	  
  EDom=M_EDo;
  if iLP==0
   EDcm=sqrt(M_EDc.x.^2+M_EDc.y.^2);
  else
   EDcm=M_EDc.x;
   EDom.y=zeros(size(EDom.x));
  end

ig=fig;
sF=size(EDcm);
  if length(sF)==3
   Eqw=EDcm(:,:,ig);
   Eout.x=EDom.x(:,:,ig);
   Eout.y=EDom.y(:,:,ig);
  elseif length(s)==2
   if i2D==3
    Eqw=EDcm(:,:);
    Eout.x=EDom.x(:,:);
    Eout.y=EDom.y(:,:);
   else
    Eqw=EDcm(:,ig);
    Eout.x=EDom.x(:,ig);
    Eout.y=EDom.y(:,ig);
   end
  end
  ADom=ADo(:,ig);
  ADcm=ADo(:,ig);
%  ADom=ADo;
%  ADcm=ADc;  
  xro=xvero;
%'qui controllo', keyboard

   PaDy.a=Wvmod;
   PaDy.b=czmod;
   PaDy.c=fqwmod;
   PaDy.d=Lfmod;
   PaDy.e=Tfmod;
   PaDy.e1=Tpemod;
   PaDy.e2=Tmemod;
   PaDy.f=Tfimod;
   PaDy.g=gtmod;
   PaDy.M2=Purspet;

   d_at=d*1e6;
   PaDy.h=fatqw;
   PaDy.i=d_at;
   PaDy.l=tyPmod;
   PaDy.m=NQW;
   PaDy.n=rr;
   PaDy.o=avero;
   PaDy.q=confztot;
   PaDy.t=famod;
   PaDy.nef=PD.nefgra;

   Pa.taut=Wvmod;   % il nome e` dummy!
   Pa.tauu=Tfmod;
   Pa.taub=Tfimod;
   Pa.taup=Tpemod;
   Pa.taum=Tmemod;
   Pa.gt=gtmod;
   Pa.dat=d_at;
   Pa.type=tyPmod;
   Pa.rmed=rg;
   Pa.rr=rr;
   
   uLtot=[];
   for kmodi=1:length(Gsov)
    if size(uLv,1)==1
       pu=1;
    else   
     if nord(1,kmodi)==-1
      pu=1;
     else
      pu=2;
     end
    end 
    uLtot=[uLtot; uLv(pu,:)];
   end
   
   Pa.confztot=uLtot;
   Pa.NQW=NQW;
   Pa.losm=gtmod;
   Pa.trasm=famod;
   Pa.Lf=Lfmod;
   Pa.Teff=Teff;
   Pa.fatqw=fatqw;
   Pa.pvol=pvol;
   Pa.E_Temp=E_Temp;
   Pa.z_Temp=z_Temp;

  delf=fmod;
  gain=gmod';
  fian=fian0;
  if ~exist('ierrMod')
   ierrMod=0;
  end
  if ierrMod==0 & exist('Cug')
  Cu=Cug;
  end
  ord=nord;
  nrAzim=immax;



%'fine mvg_last', keyboard
%'fine mvg_last', keyboard
