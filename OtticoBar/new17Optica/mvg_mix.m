function [Eqw,Eout,xro,fmod,gmod,nord,immax,lambda,PaDy]=...
          mvg_mix(nomepa,Ndis,Tdis,rodisN,rodisT,zedis0,xroJ0);
format compact
format short e

%% Costanti universali

global  h hbar kB mi c eps0 Z0 mass0 q
mass0=9.1e-31;
h=6.626e-34;
hbar=h/(2.*pi);
kB=1.38e-23;
j=sqrt(-1);
mi=pi*4e-7;
eps0=8.8541e-12;
c=1/sqrt(mi*eps0);
Z0=120*pi;
q=1.6e-19;

global ilo
global CP W Lemi Lef fatqw confz fTras fTrasinf fPsup fPinf confztot
global r_pil r_met
global par_in fimaxi lfi_inp fian0 Cug
global tauN Bric Cric Pumv
global r1r r2r r0r nstrau nstrad L_cav fconv ra  L_i n_i a_i fiQ numoradv
global confz0 confz0t
global PD
global imm

global shape desi i2D

global yiN yiT xiN xiT isuT zedis anirid isu
global xroJ
xroJ=xroJ0;
zedis=zedis0;
global Luvb Luv Lbb Lbt Lvbr Lvtr nuvb nuv nbb nbt nvbr nvtr  Lvt nvt Lvb nvb
global Fiu Xiu Fid Xid
global Ndt dgra d_cap n_ext  emme icasa iloelet
global deltanyb deltanyn iany iinf iallin r41 numodiacc fattany aniat
global aniats deltanyns deltanybs Strain aniatsv deltanynsv deltanybsv ianys
global igamveb igamveu GGbext GGuext
global Frisims Frisums alimms mmvets
global Lv rv Lb rb Litot nitot igrad iperat fcrit
global ndisg Lf Gaf2 L1g L2g L3g L4g nfe10 rfe r11 r12 nfe1 ifem ifeed iaut
global deltany iliplan iafr
global isopro ikiaut ifpdum
global nextvet istrutpil istrutpilcav istrutpilme
global gapla iztm icR
global pia piaext igrac Nstu
global nvert dap daq dae daa fc rc
global nvar nk1max alim alimi ilat ipolar nubesi KKv nubesu iprvol
global niat aiat nmetv posoib istrmet
global r1su r2su r1sb r2sb
global iemmeg ilin
global ibes mm  nmod iperd nl perd
global pardiff ia iaext
global xi yi nst lxi yiv xs ldiff x1 y1 ifdiff kcav kcav0 lambda
global xvn0 yvn0
global iload  ifp iLP iprimo  iriga
global Temp k0 lambda NQW d nstratiu nstratid
global aplan avero  Lcav aos2 amet a0ref
global bvero ras0
global Lcav1 Lcav11 Lcav12 Lcav13 Lcav2 Lcav21 Lcav22 Lcav23 ross ross1
global dos L21 L22 L2t L1t ceoss
global m0 m1 m2 ma rr rqw0 r0 r1 r2 rfu rfd rforo rfum
global aniatv deltanynv deltanybv Fist zsta
global flagthrougt Bias pucavi rqw isu ipuos
global mbv pasnu KK npk nup1 nube
global DYo
global pis10 pis20 pic0
global igau ite isim
global xi yi lxi ldiff  x1 y1 ifdiff Ldiff fatdiff icorrente
global ioss iposs idiffos ral ifn
global igint numodi ivar itetm ivarext
global xm1 xm2 xmc alphap N00 dr10 dr20 dr00 drforo0 r10 r20 r00 rr0 rforo0
global Lic Lis Liu nic nis niu aic ais aiu atot igen ivfre nirid istt ameti
global Litn aitn nitn fmlst Lib aib nib fmlsb
global Litn_0 aitn_0 nitn_0 nib_0 aib_0 niat_0 aiat_0
global ifpil ipilat ipsup istox
global coeffgv coeffnv sgal falf g0 ra epsatt me mh
global Frisi Frisu ifit Nmax Ndisp
global x1 y1 xs ys
global nmet next istrutt ipsup iplan
global next1 next2 next1c next1s
global isi itemperat apila Tvet fatd apilat
global alphapf alphap0 alphap1 alphap2 alphapc
global Dla ios iosr
global ifnm
global a_Lu a_LHu a_HLu a_Hu n_Lcu n_LHu n_HLu n_Hcu L_Lcu L_LHu L_HLu L_Hcu Fiu Xiu
global a_Ld a_LHd a_HLd a_Hd n_Lcd n_LHd n_HLd n_Hcd L_Lcd L_LHd L_HLd L_Hcd Fid Xid


eval(['load ' nomepa]);
nvar=nvar1;
Frisi=Frisi0;
Frisu=Frisu0;
Ndisp=Ndisp0;

immax=length(mmvet);

 if exist('igamveb')==0 | length(igamveb)==0
  igamveb=1;
 end
 if exist('igamveu')==0 | length(igamveu)==0
  igamveu=1;
 end

 if exist('fcrit')==0 | length(fcrit)==0
  fcrit=1;
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

 eval(nvar)
 avero=a0ref;
 ra=niat(1);

 k0=2*pi/lambda0*1e-6;
 igen=1;

 kcav=2*pi/lambda*rr;
 kcav0=kcav;
 kcav00=kcav*1e6;
 a=kcav0*avero;     %adimensionale:normalizzato rispetto a k0 della cavita'


igint=1;


global nomes1
nomes1=nomes;
imetrice=0;

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
         if igau==4
          xiT=rodisT(2:length(rodisT));
          s=size(Tdis);
          for iz=1:s(2)
           yiT(:,iz)=-diff(Tdis(:,iz));
          end
         end

         yiN=[];
         if length(Ndis)>1
          xiN=rodisN(2:length(rodisN));
          NN=Ndis/max(Ndis);
          yiN=-diff(NN);
         end




imm=0;
for mm=mmvet
 imm=imm+1;
  icalcola=1;




     [KK,ADc,ADo,gsov,fsov,nazim,nrad,xvero,M_EDc,M_EDo,imod,tyE]=...
     sub_mix(icalcola);
%     disp(' dopo sub_mix '), keyboard

     if iLP==1
      EDo.x=M_EDo.x;
      EDo.y=M_EDo.x*0;
      EDc=M_EDc.x;
     else
      EDo.x=M_EDo.x;
      EDo.y=M_EDo.y;
      EDc=M_EDc.x+M_EDc.y;
     end

     if imm==1
      siE=0;
     else
      sE=size(EDcm);
      if length(sE)==3
       siE=sE(3);
      else
       siE=sE(2);
      end
     end
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
       ADom(imm,s1k,s2k)=ADo;
       ADcm(imm,s1k,s2k)=ADc;
       s=1:length(gsov);
       nmodm(imm,s)=ones(size(gsov))*mm;
       gsovm(imm,s)=gsov;
       tyPm(imm,s)=tyE;
       fsovm(imm,s)=fsov;
       nazm(imm,s)=nazim;
       nram(imm,s)=nrad;
%       disp(' m.sub'), keyboard
       Wvm(imm,s)=PD.a;
       czm(imm,s)=PD.b;
       fqwm(imm,s)=PD.c;
       Lfm(imm,s)=PD.d;
       Tfm(imm,s)=PD.e;
       Tfim(imm,s)=PD.f;
       gtm(imm,s)=PD.g;
       glum(imm,s)=PD.r;
       glbm(imm,s)=PD.s;
       fam(imm,s)=PD.t;


end   %mm

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
glumo=[];
glbmo=[];
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
 glumo=[glumo reshape(glum(imm,:),1,sp(2))];
 glbmo=[glbmo reshape(glbm(imm,:),1,sp(2))];
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
glumod=glumo(fig);
glbmod=glbmo(fig);

nord=[namod; nrmod];
if length(s)==3
 Eqw=EDcm(:,:,ig);
 Eout.x=EDom.x(:,:,ig);
 Eout.y=EDom.y(:,:,ig);
elseif length(s)==2
 Eqw=EDcm(:,ig);
 Eout.x=EDom.x(:,ig);
 Eout.y=EDom.y(:,ig);
end
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
PaDy.r=glumod;
PaDy.s=glbmod;
PaDy.t=famod;

disp(' fine mvgsub ')
%keyboard
if ifp>-4, keyboard, end
