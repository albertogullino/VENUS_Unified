xroJ=xroI;
deltany=0;



immax=length(mmvet);

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

         yiT=[];
         if igau==4 & T~=0
          xiT=rodisT(2:length(rodisT));
          s=size(Tdis);
          for iz=1:s(2)
           yiT(:,iz)=-diff(Tdis(:,iz));
          end
         end




imm=0;
imod=0;
for mm=mmvet
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





     sub
%' DOPO sub' ,keyboard, keyboard

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

       siE=0;
%      if imm==1
%       siE=0;
%      else
%       sE=size(EDcm);
%       if length(sE)==3
%        siE=sE(3);
%       else
%        siE=sE(2);
%       end
%      end


%     if sha>1
       ske=size(EDc);
       if i2D==3
        s1e=1:ske(1);
        s2e=1:ske(2);
        if imod==1
         s3e=1;
        else
         s3e=1:ske(3);
        end
         EDom.x(s1e,s2e,s3e+siE)=EDo.x;
         EDom.y(s1e,s2e,s3e+siE)=EDo.y;
         EDcm(s1e,s2e,s3e+siE)=EDc;
       else
%        disp(' per ora non fatto '), keyboard
         s1e=1:ske(1);
         s2e=1:ske(2);
         EDom.x(s1e,s2e+siE)=EDo.x;
         EDom.y(s1e,s2e+siE)=EDo.y;
         EDcm(s1e,s2e+siE)=EDc;
       end
       KKm(imm,1:length(KK))=KK';
       sk=size(ADo);
       s1k=1:sk(1);
       s2k=1:sk(2);
%       ADom(imm,s1k,s2k)=ADo;
%       ADcm(imm,s1k,s2k)=ADc;
       ADom(s1k,s2k+siE)=ADo;
       ADcm(s1k,s2k+siE)=ADc;
       s=1:length(gsov);
       nmodm(imm,s)=ones(size(gsov))*mm;
       tyPm(imm,s)=tyE;
       nazm(imm,s)=nazim;
       nram(imm,s)=nrad;
       Wvm(imm,s)=PD.a;
       czm(imm,s)=PD.b;
       fqwm(imm,s)=PD.c;
       Lfm(imm,s)=PD.d;
       Tfm(imm,s)=PD.e;
       Tfim(imm,s)=PD.f;
       gtm(imm,s)=PD.g;
%       glum(imm,s)=PD.r;
%       glbm(imm,s)=PD.s;
       fam(imm,s)=PD.t;

     end  %icampi
       gsovm(imm,s)=gsov;
       fsovm(imm,s)=fsov;


end   %mm


if icampi>0
  clear Eqw Eout
  nubev=[1 0 2:20];
  Ep=[];
  gmo=[];
  fmo=[];
  nrmo=[];
  namo=[];
  Wvmo=[];
  famo=[];
  czmo=[];
  fqwmo=[];
  Lfmo=[];
  Tfmo=[];
  Tfimo=[];
  gtmo=[];
%  glumo=[];
%  glbmo=[];
  tyPmo=[];
  ord=[];
  imm=0;
  s=size(EDom.x);
  sp=size(gsovm);
  for mm=mmvet(1:immax)
   mmd=mm+1;
   imm=imm+1;
   if iLP==0
    nube=nubev(mmd);
   else
    nube=imm;
   end
   gmo=[gmo reshape(gsovm(imm,:),1,sp(2))];
   tyPmo=[tyPmo reshape(tyPm(imm,:),1,sp(2))];
   fmo=[fmo reshape(fsovm(imm,:),1,sp(2))];
   nrmo=[nrmo reshape(nram(imm,:),1,sp(2))];
   namo=[namo reshape(nazm(imm,:),1,sp(2))];
   Wvmo=[Wvmo reshape(Wvm(imm,:),1,sp(2))];
   famo=[famo reshape(fam(imm,:),1,sp(2))];
   czmo=[czmo reshape(czm(imm,:),1,sp(2))];
   fqwmo=[fqwmo reshape(fqwm(imm,:),1,sp(2))];
   Lfmo=[Lfmo reshape(Lfm(imm,:),1,sp(2))];
   Tfmo=[Tfmo reshape(Tfm(imm,:),1,sp(2))];
   Tfimo=[Tfimo reshape(Tfim(imm,:),1,sp(2))];
   gtmo=[gtmo reshape(gtm(imm,:),1,sp(2))];
%   glumo=[glumo reshape(glum(imm,:),1,sp(2))];
%   glbmo=[glbmo reshape(glbm(imm,:),1,sp(2))];
  end

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
%  glumod=glumo(fig);
%  glbmod=glbmo(fig);

  nord=[namod; nrmod];
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
%   PaDy.r=glumod;
%   PaDy.s=glbmod;
   PaDy.t=famod;

%   Pa.W=Wvmod;
%   Pa.cz=czmod;
%   Pa.fqw=fqwmod;
%   Pa.Lf=Lfmod;
%   Pa.Tf=Tfmod;
%   Pa.Tfinf=Tfimod;
%   Pa.gt=gtmod;
%   Pa.fatqw=fatqw;
%   Pa.dat=d_at;
%   Pa.type=tyPmod;
%   Pa.NQW=NQW;
%   Pa.rr=rg;
%   Pa.avero=avero;
%   Pa.conztot=confztot;
%   Pa.glu=glumod;
%   Pa.glb=glbmod;
%   Pa.famod=famod;

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

%   Pa.Gam_u=abs((rfu-nitn(1,1))/(rfu+nitn(1,1)));


  delf=fmod;
  gain=gmod';
  fian=fian0;
  Cu=Cug;
  ord=nord;
  nrAzim=immax;

  if ifp~=-4
  disp(' fine mvg ')
  end
%  keyboard
  if ifp>-4, keyboard, end
else
  gmo=[];
  fmo=[];
  sp=size(gsovm);
  imm=0;
  for mm=mmvet(1:immax)
   mmd=mm+1;
   imm=imm+1;
   if iLP==0
    nube=nubev(mmd);
   else
    nube=imm;
   end
   gmo=[gmo reshape(gsovm(imm,:),1,sp(2))];
   fmo=[fmo reshape(fsovm(imm,:),1,sp(2))];
  end

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