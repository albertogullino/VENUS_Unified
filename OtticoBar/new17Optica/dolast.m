if ~exist('Plots')
Plots=Plot;
end
  inval=find(Gsov<=0);
%  Gdu=Gsov.*M2v;
  Gdu=Gsov;
  Gdu(inval)=1e100;
  [du,idu]=sort(Gdu);
  ival=find(du<1e100);
  iso=idu(ival);

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
  Plot.E2xp{kco}=Plots.E2xp{kis};		
  if iLP==0
   Plot.E2yo{kco}=Plots.E2yo{kis};		
   Plot.E2yp{kco}=Plots.E2yp{kis};		
   Plot.E2zo{kco}=Plots.E2zo{kis};		
   Plot.E2zp{kco}=Plots.E2zp{kis};		
  end 
  Plot.Ef{kco}=Plots.Ef{kis};			
  Plot.Cug{kco}=Plots.Cug{kis};		
  Plot.gou{kco}=Plots.gou{kis};		
  Plot.aou{kco}=Plots.aou{kis};		
  Plot.fou{kco}=Plots.fou{kis};		
  Plot.ze{kco}=Plots.ze{kis};			
  Plot.gg0{kco}=Plots.gg0{kis};		
  Plot.FF{kco}    =Plots.FF{kis};
  Plot.parmod{kco}=Plots.parmod{kis};

 end 					

 
 if ifp==-100
 'controllo Plot', keyboard
 'controllo Plot', keyboard
 'controllo Plot', keyboard
 end
 ADc=           ADc(:,iso);			
   ADout=        ADout(:,iso);			
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
end
     
if ifp==-100
' qui mvg_last ', keyboard
' qui mvg_last ', keyboard
' qui mvg_last ', keyboard
end
% iLP=iLP1;
if ~exist('gsov')
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

gmo=gsov;
fmo=fsov;

%'qui', keyboard
tyPmo=tyE;
nrmo=nazim;
namo=nrad;
Wvmo=PD.a;
czmo=PD.b;
fqwmo=PD.c;
Lfmo=PD.d;
Tfmo=PD.e;
Tfimo=PD.f;
gtmo=PD.g;
famo=PD.t;

  fi=find(gmo>0);
  gmo1=gmo(fi);
  [gmod,fig]=sort(gmo1);
  [gdu,ig]=sort(gmo(fi));
  fig=fi(ig)';
'qui', keyboard

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
  Tfimod=Tfimo(fig);
  gtmod=gtmo(fig);

  nord=[namod; nrmod];
  ADom=ADo;
  ADcm=ADc;
  EDom=M_EDo;
  if iLP==0
   EDcm=sqrt(M_EDc.x.^2+M_EDc.y.^2);
  else
   EDcm=M_EDc.x;
   EDom.y=zeros(size(EDom.x));
  end


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
  ADom=ADom(:,ig);
  ADcm=ADcm(:,ig);
  xro=xvero;


   PaDy.a=Wvmod;
   PaDy.b=czmod;
   PaDy.c=fqwmod;
   PaDy.d=Lfmod;
   PaDy.e=Tfmod;
   PaDy.f=Tfimod;
   PaDy.g=gtmod;

   d_at=d*1e6;
   PaDy.h=fatqw;
   PaDy.i=d_at;
   PaDy.l=tyPmod;
   PaDy.m=NQW;
   PaDy.n=rr;
   PaDy.o=avero;
   PaDy.q=confztot;
   PaDy.t=famod;


   Pa.taut=Wvmod;   % il nome e` dummy!
   Pa.tauu=Tfmod;
   Pa.taub=Tfimod;
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
  Cu=Cug;
  ord=nord;
  nrAzim=immax;
 figm=fig;
inews=1;
if inews==0
  Plotme    =Plot;
 for ig=1:length(figm)
  fig=figm(ig);
  Plot.FF{ig}    =Plotme.FF{fig};
  Plot.gou{ig}   =Plotme.gou{fig};
  Plot.aou{ig}   =Plotme.aou{fig};
  Plot.fou{ig}   =Plotme.fou{fig};
  Plot.ze{ig}    =Plotme.ze{fig};
  Plot.gg0{ig}   =Plotme.gg0{fig};
  Plot.Ap{ig}    =Plotme.Ap{fig};
  Plot.ApQ{ig}   =Plotme.ApQ{fig};
  Plot.parmod{ig}=Plotme.parmod{fig};
  Plot.XP{ig}    =Plotme.XP{fig};
  Plot.YP{ig}    =Plotme.YP{fig};
  Plot.X{ig}     =Plotme.X{fig};
  Plot.Y{ig}     =Plotme.Y{fig};
  Plot.E2xo{ig}  =Plotme.E2xo{fig};
  Plot.E2xp{ig}  =Plotme.E2xp{fig};
  if iLP==0
   Plot.E2yo{ig} = Plotme.E2yo{fig};
   Plot.E2yp{ig} = Plotme.E2yp{fig};
   Plot.E2zo{ig} = Plotme.E2zo{fig};
   Plot.E2zp{ig} = Plotme.E2zp{fig};
  end
  Plot.Ef{ig}    =Plotme.Ef{fig};
  Plot.Cug{ig}   =Plotme.Cug{fig};
 end
end  %inews
% save sa Plotme Plot
'fine mvg_last'
return
if ifp==-10
 keyboard
 keyboard
end