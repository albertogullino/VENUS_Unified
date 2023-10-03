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
 'controllo Plot', keyboard


  inval=find(Gsov<=0);
%  Gdu=Gsov.*M2v;
  Gdu=Gsov;
  Fdu=Fsov;
  Gdu(inval)=1e100;
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
  InR=IndRad(iso);
  if nmasce<0
   if length(iso)>abs(nmasce)
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
  ' nmasce', keyboard
  %iso=isu;
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
  Plot.FF{kco}=Plots.FF{kis};			
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
  Plot.FF{kco}    =Plots.FF{kis};
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
   EDc.y=      EDc.y(:,:,iso);
   EDout.y=   EDout.y(:,:,iso);
      
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
if ifp==-10
' qui mvg_last ', keyboard
end

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


tyPmo=tyE;
nrmo=nazim;
namo=nrad;
Wvmo=PD.a;
czmo=PD.b;
fqwmo=PD.c;
Lfmo=PD.d;
Tfmo=PD.e;
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
ig=iso;
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
  if length(s)==3
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
%  ADom=ADom(:,ig);
%  ADcm=ADcm(:,ig);
  ADom=ADo;
  ADcm=ADc;  
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
   Pa.rr=rg;
   Pa.confztot=confztot;
   Pa.NQW=NQW;
   Pa.losm=gtmod;
   Pa.trasm=famod;
   Pa.Lf=Lfmod;


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
%'fine mvg_last', keyboard
