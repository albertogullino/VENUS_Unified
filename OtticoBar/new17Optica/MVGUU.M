%' passo mvguu ', keyboard

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

%'ipla', keyboard
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
 Pa.glotp=dglo;




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
          if igau==4 & length(T)>1
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
           'in mvguu '
           figure, plot(abs(cumsum(yiT',2)),'.'), pausak
          end
         end
        end
        if sN==2
         yiN=[1:100]; %solo per renderlo un vettore
        end

imm=0;
imod=0;
for mm=mmvet
 Frisi=Dso(1)*1e-3/lambda;
 Frisu=Dso(2)*1e-3/lambda;
 if (is_even(mm)==1 & iLP==0) | (is_even(mm)==0 & iLP==1)
  ' Frisi ', keyboard
  Frisi=Frisi+DelFr;
  Frisu=Frisu+DelFr;
 end
 mm_ver=mm;
 imm=imm+1;
  clear KAp KA icousav
  icalcola=1;
  igainshape=0;
%  '[icalcola mm]'
%  [icalcola mm]
%  pausak
% if length(mmvet)>1
% imod=0;
% end



     if imod_err==1
      return
     end
%     ' sub ', keyboard

     sub
   Ppo{imm}=Ppol;
   Ppt{imm}=Plot;

end %mm

     if icampi>0
      if iLP==1
       EDo.x=M_EDo.x;
       EDo.y=M_EDo.x*0;
       EDc=M_EDc.x;
      else
       EDo.x=M_EDo.x;
       EDo.y=M_EDo.y;
       EDc=M_EDc.x+M_EDc.y;
      end


       ske=size(EDc);
       if i2D==3
        s1e=1:ske(1);
        s2e=1:ske(2);
        if imod==1
         s3e=1;
        else
         s3e=1:ske(3);
        end
         EDom.x(s1e,s2e,s3e)=EDo.x;
         EDom.y(s1e,s2e,s3e)=EDo.y;
         EDcm(s1e,s2e,s3e)=EDc;
       else
%        disp(' per ora non fatto '), keyboard
         s1e=1:ske(1);
         s2e=1:ske(2);
         EDom.x(s1e,s2e)=EDo.x;
         EDom.y(s1e,s2e)=EDo.y;
         EDcm(s1e,s2e)=EDc;
       end
%       KKm(imm,1:length(KK))=KK';
       sk=size(ADo);
       s1k=1:sk(1);
       s2k=1:sk(2);
       ADom(s1k,s2k)=ADo;
       ADcm(s1k,s2k)=ADc;

       s=1:length(gsov);

%       nmodm(imm,s)=ones(size(gsov))*mm;
       tyPmo=tyE;
       namo=nazim;
       nrmo=nrad;
       Wvmo=PD.a;
       czmo=PD.b;
       fqwmo=PD.c;
       Lfmo=PD.d;
       Tfmo=PD.e;
       Teffmo=PD.eu;
       Tfimo=PD.f;
       Tfmo1=PD.e1;
       Tfimo1=PD.f1;
       gtmo=PD.g;
       famo=PD.t;

     end  %icampi
       gmo=gsov;
       fmo=fsov;



if icampi>0
  clear Eqw Eout
  nubev=[1 0 2:20];

  fi=find(gmo>0);
  [gmod,ig]=sort(gmo(fi));
  fig=fi(ig)';
  fmod=fmo(fig);
  tyPmod=tyPmo(fig);
  nrmod=nrmo(fig);
  namod=namo(fig);
  Wvmod=Wvmo(fig);
  famod=famo(fig);
  czmod=czmo(fig);
  fqwmod=fqwmo(fig);
  Lfmod=Lfmo(fig);
  Tfmod=Tfmo(fig);
  Teffmod=Teffmo(fig);
  Tfimod=Tfimo(fig);
  Tfmod1=Tfmo1(fig);
  Tfimod1=Tfimo1(fig);
  gtmod=gtmo(fig);

  nord=[namod; nrmod];
  s=size(EDcm);
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

%  if idyn==1

   PaDy.a=Wvmod;
   PaDy.b=czmod;
   PaDy.c=fqwmod;
   PaDy.d=Lfmod;
   PaDy.e=Tfmod;
   PaDy.ee=Teffmod;
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
   Pa.pvol=Tfmod1;
   Pa.alca=Tfimod1;
   Pa.gt=gtmod;
   Pa.dat=d_at;
   Pa.type=tyPmod;
   Pa.rr=rg;
   Pa.aiat=Ppol.aiat;
%   ' confztot', keyboard
   Pa.confztot=confzv;
   Pa.NQW=NQW;
   Pa.losm=gtmod;
   Pa.trasm=famod;
   Pa.Lf=Lfmod;
   Pa.Teff=Teffmod;
   Pa.fatqw=fatqw;


  delf=fmod';
  gain=gmod';
  fian=fian0;
  Cu=Cug;
  ord=nord;
  nrAzim=immax;


  disp(' fine mvg ')
  keyboard
  if ifp>-4, keyboard, end
else

  fi=find(gmo>0);
  [gmod,ig]=sort(gmo(fi));
  fig=fi(ig)';
  fmod=fmo(fig)';

  tyPmod=zeros(size(fmod));
  nrmod=zeros(size(fmod));
  namod=zeros(size(fmod));
  Wvmod=zeros(size(fmod));
  famod=zeros(size(fmod));
  czmod=zeros(size(fmod));
  fqwmod=zeros(size(fmod));
  Lfmod=zeros(size(fmod));
  Tfmod=zeros(size(fmod));
  Tfimod=zeros(size(fmod));
  gtmod=zeros(size(fmod));
  glumod=zeros(size(fmod));
  glbmod=zeros(size(fmod));

  nord=[namod; nrmod];


  delf=fmod;
  gain=gmod';



end
