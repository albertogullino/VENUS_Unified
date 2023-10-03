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
% ' set_stru ', pausak
 set_stru




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


%'mvg', keyboard


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
 if length(Dlam_mod)==4
  DelFr=Dlam_mod(4)*1e-3/lambda;
 else
  DelFr=0;
 end




 avero=a0ref;
 %ra=niat(1);

 k0=2*pi/lambda;

 kcav=2*pi/lambda*rr;
 kcav0=kcav;
 kcav00=kcav*1e6;
 a=kcav0*avero;     %adimensionale:normalizzato rispetto a k0 della cavita'


 igint=1;   % coupl. coeff. integration



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
           s=size(Tdis);
           for iz=1:s(2)
            yiT(:,iz)=-diff(Tdis(:,iz));
           end
          end
          zdis=(zdis(1:end-1)+zdis(2:end))/2;
          if ifp~=-4
           'in mvg_chaos '
           figure, plot(abs(cumsum(yiT',2)),'.'), pausak
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
%' cont', keyboard
for mm=mmvet
mmsav=mm;
%'kazzo', keyboard
 if iLP==0 & mm==0 & ivacc(1)~=-1
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
% ' prima sosta', keyboard
 for itet=1:length(itetmv)
   itetm=itetmv(itet);
   iLP=iLPr;

% ' seconda sosta', keyboard
   Frisi=Dso(1)*1e-3/lambda;
   Frisu=Dso(2)*1e-3/lambda;
   if (mm==0 & iLP==0) | (mm==1 & iLP==1)
    Frisi=Frisi+DelFr;
    Frisu=Frisu+DelFr;
   end
   mm_ver=mm;
   imm=imm+1;
  clear KAp KA icousav
  icalcola=1;
  igainshape=0;



     if imod_err==1
      return
     end

     sub
   Ppo{imm}=Ppol;
   Ppt{imm}=Plot;
 end % itetm
end   %mm

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
 save sa Plotme Plot
