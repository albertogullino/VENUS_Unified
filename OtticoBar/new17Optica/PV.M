gmo=gsov;
fmo=fsov;
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
  [gmod,ig]=sort(gmo(fi));
  fig=fi(ig)';
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
  Plotme    =Plot;
 for ig=1:length(figm)
  fig=figm(ig);
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
  end
  Plot.Ef{ig}    =Plotme.Ef{fig};
  Plot.Cug{ig}   =Plotme.Cug{fig};
 end
