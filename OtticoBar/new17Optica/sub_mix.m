function   [KK,ADc,ADout,Gsov,Fsov,Nazim,Nrad,xvero,EDc,EDout,imod,tyE]=...
          sub_mix(icalcola);

%pausak

%in entrata
global tauN Bric Cric Pumv
global r1r r2r r0r nstrau nstrad L_cav fconv ra  L_i n_i a_i fiQ numoradv
global confz0 confz0t
global CP W Lemi Lef fatqw confz fTras fTrasinf fPsup fPinf confztot

global ilo
global r_pil r_met
global par_in fimaxi lfi_inp fian0 Cug
global iraff idyn imm
global shape desi i2D
global xroJ
global yiN yiT xiN xiT  isuT zedis anirid isu
global Luvb Luv Lbb Lbt Lvbr Lvtr nuvb nuv nbb nbt nvbr nvtr  Lvt nvt Lvb nvb
global Fiu Xiu Fid Xid
global Ndt dgra d_cap n_ext  emme icasa iloelet
global a_Lu a_LHu a_HLu a_Hu n_Lcu n_LHu n_HLu n_Hcu L_Lcu L_LHu L_HLu L_Hcu Fiu Xiu
global a_Ld a_LHd a_HLd a_Hd n_Lcd n_LHd n_HLd n_Hcd L_Lcd L_LHd L_HLd L_Hcd Fid Xid
global deltanyb deltanyn iany iinf iallin r41 numodiacc fattany aniat
global aniats deltanyns deltanybs Strain aniatsv deltanynsv deltanybsv ianys
global igamveb igamveu GGbext GGuext
global Frisims Frisums alimms mmvets
global Lv rv Lb rb Litot nitot igrad
global ifnm igrac Nstu iperat fcrit
global ndisg Lf Gaf2 L1g L2g L3g L4g nfe10 rfe r11 r12 nfe1 ifem ifeed iaut
global deltany iliplan iafr
global isopro ikiaut ifpdum
global igau iosr ios
global bvero ras0
global nextvet istrutpil istrutpilcav istrutpilme
global gapla iztm icR
global ameti amete apila
global nvert dap daq dae daa rc fc
global r1su r2su r1sb r2sb posoib
global nvar nk1max alim alimi ilat nmetv rfum
global ifpil ipilat ipsup istox ipolar nubesi KKv nubesu itetm
global iemmeg ilin iprvol
global Litot aitot nitot istrmet
global isi itemperat apila Tvet fatd apilat
global nomes1
global x1 y1 xs ys iplan xvn0 yvn0
global Frisi Frisu ifit Nmax Ndisp
global ibes mm iload iprimo
global pardiff ia iaext
global ifp iLP iriga
global iperd perd
global nstratiu nstratid k0
global coeffgv coeffnv sgal falf g0 ra m0 lav port
global rr rqw0 r0 r1 r2 rfu rfd rforo d NQW avero Lcav h hbar kB Temp c aplan kcav kcav0
global flagthrougt Bias pucavi rqw isu ipuos
global aniatv deltanynv deltanybv Fist zsta
global epsatt eps0 Z0 lambda lambdac mass0
global me mh
global xi yi lxi nst yiv xs ldiff kcav x1 y1 ifdiff Ldiff fatdiff icorrente
global ioss iposs idiffos ral
global nmod nl
global dos L21 L22 L2t L1t ceoss
global mbv pasnu KK npk nup1 nube
%global Lib Lit nib nit aib ait igen nirid ipsup
global Lic Lis Liu nic nis niu aic ais aiu atot igen ivfre nirid ipsup istt ameti
global Litn aitn nitn fmlst Lib aib nib fmlsb
global Litn_0 aitn_0 nitn_0 nib_0 aib_0 niat_0 aiat_0


%in uscita
global pis pic
global pis10 pis20 pic0

global Lcav1 Lcav11 Lcav12 Lcav13 Lcav2 Lcav21 Lcav22 Lcav23 ross

global igint numodi ivar aos2 amet ivarext
global xm1 xm2 xmc alphap N00 dr10 dr20 dr00 drforo0 r10 r20 r00 rr0 rforo0
global m0 m1 m2 ma rr r0 r1 r2 rfu rfd rforo
global nme next istrutt ipsup rforo ross1
global alphapf alphap0 alphap1 alphap2 alphapc
global Dla
global next1 next2 next1c next1s
global niat aiat

%icasa=0;
T0=Temp;

iprvol=0;
KKsav=0;
npstop=3;
ibw=0;
ifpil0=2;
vgv=50*c/rr;
%keyboard

global ifn


   clear Anuvet Acavet Amevet Gvet alvet Fint Kvet
   clear Anuvet1 Acavet1 Amevet1 Gvet1 alvet1

%%%%%%%%%%%%%%%%  PARTE 2: Cerco la frequenza di risonanza KoL, pic, fris:
%                         I)  nel caso planare avero>aplan la fisso direttamente
%                             a KoL= pic=pi, fris=0 (quella di costruzione)
%                         II) nel caso non planare avero<aplan, ripeto la
%                             soluzione dell'intero problema npfr volte e poi
%                             trovo la fris con l'autoconsistenza in autoc.m
%
%%%%%%%%%%%%%% II]

' in submix', keyboard
' in submix', keyboard
 nube=mm;

nupd=nube/2-fix(nube/2);

     pasnu=2;
     nubesi=nube-pasnu*numodiacc;
     if nubesi<0
      if abs(nube/2-fix(nube/2))==0
       nubesi=0;
      else
       nubesi=1;
      end
     end
     nubesu=nube+pasnu*numodiacc;

     numodi=(nubesu-nubesi)/pasnu+1;


ibr=1;


  fre_camp=linspace(Frisi,Frisu,Ndisp);




  nk1=nk1max;
  eint=alim;


  if iLP==1
   nni=numodi;
  else
   nni=2*numodi;
  end

  pes0=1;


   eval(nvar);

for ifr=1:length(fre_camp)

    freq=fre_camp(ifr);

    [ifr freq]

    kcav=kcav0*(1+freq);

%        if ivfre==1
%         frequer=freq;
%        end

%      disp(' prima acc0 ')

      tic
      acc0
      toc

% N1i=-inv(N1);
% Mtot=N1i*N2;
% mapab(Mtot)

%      disp(' dopo acc0  in sub_mix')
%      keyboard


   Fint(ifr)=freq;
   sAz=size(Anz);

    Azvet(1:sAz(1),1:sAz(2),1:sAz(3),ifr)=Anz;
    Azvetf(1:sAz(1),1:sAz(2),1:sAz(3),ifr)=Anzf;


    Gvet(1:sAz(2),ifr)=Gvav;
    sk=size(KK);
    Kvet(1:max(sk),ifr)=KK;
    alvet(1:sAz(2),ifr)=alphavv;

    if ipolar==2

     sAz=size(Anz1);

     Azvet1(1:sAz(1),1:sAz(2),1:sAz(3),ifr)=Anz1;
     Azvetf1(1:sAz(1),1:sAz(2),1:sAz(3),ifr)=Anzf1;

     Gvet1(1:sAz(2),ifr)=Gvav1;
     alvet1(1:sAz(2),ifr)=alphavv1;
    else
     Azvet1=0;
     Azvetf1=0;
     Gvet1=0;
     alvet1=0;
    end

end % while ifr

%  disp(' fuori dal while ')
%  keyboard
%  if iriga==0 & ifp>=-3
   ifps=ifp;
   if ifp>=-2
    ifp=-1;
   end

  if ipolar==2
   pola=1;
  else
   pola=ipolar;
  end


  Azvet0=Azvet;
  Gvet0=Gvet;
  alvet0=alvet;

  save O

  imod=0;
  camv_new
  if ifp>=-3 | ifp==-10
   figure,
   puf=1:length(fso);
   subplot(211), plot(fou(puf,:)',aou(puf,:)',fso,fso*0,'wo'), grid;
   dista=(max(gou(puf,:))-min(gou(puf,:)))/min(gou(puf,:));
   if dista>10
    subplot(212), semilogy(fou(puf,:)',gou(puf,:)'/vg,fso,gso/vg,'wo');
   else
    subplot(212), plot(fou(puf,:)',gou(puf,:)'/vg,fso,gso/vg,'wo');
   end

   keyboard
  end

% icasa=1;

 if ipolar==2
  Azvet=Azvet1;
  Gvet=Gvet1;
  alvet=alvet1;
  pola=-1;
  camv_new
 end
 ifp=ifps;

%  end

 disp(' autdiff1 ')
% keyboard
 if icampi>2
  if ifp>-4
   keyboard
  end
 end

global PD

 if idyn==1
  PD.a=Wv;
  PD.b=czv;
  PD.c=fqwv;
  PD.d=Lefv;
  PD.e=Tefv;
  PD.f=Tefiv;
  PD.g=gtot;
  PD.r=Perdu;
  PD.s=Perdb;
  PD.t=fa4v;
 else
  PD.a=Gsov*0;
  PD.b=Gsov*0;
  PD.c=Gsov*0;
  PD.d=Gsov*0;
  PD.e=Gsov*0;
  PD.f=Gsov*0;
  PD.g=Gsov*0;
  PD.r=Gsov*0;
  PD.s=Gsov*0;
  PD.t=Gsov*0;
 end

ierr=0;
if iriga==1
 if npfr~=0
  modaut1
 else
  fris=fr;
 end
 if exist('gg0')==0
  ierr=1;
 end
% disp(' dopo modaut')
% pausak
 if ierr==0
  if ifnm==0
   riga
  else
   save mop
   mop1
   disp(' verifica qui ')
   keyboard
%   autdiff1
  end
   if ifeed==1
    subrife
   else
    Pfi=0;
    Pfi1=0;
    pi0=0;
   end
 else
  dnup=0;
  ae=0;
  dnupv=0;
  aev=0;
  gsoav=0;
  fsoav=0;
  asoav=0;
  Nth=0;
  fTras=0;
  Pfi=0;
  Pfi1=0;
  pi0=0;
  Ww=0;
  Is=0;
  Nno=0;
  Psi=0;
  Alfp=0;
  fTras=0;
  nsp=0;
 end
else
 dnup=0;
 ae=0;
 dnupv=0;
 aev=0;
end
if ifp>=-2
 pausak
end
if ierr==1
 disp(' ierr=1 ')
 pausak
end
% figure, subplot(121), plot(fr,alvet), subplot(122), semilogy(fr,Gvet)
% figure, subplot(121), plot(fr,alvet'), subplot(122), semilogy(fr,Gvet')
%
% Ab=abs(EDc).^2;
% Abc=Ab/max(Ab);
% Ab=abs(EDo).^2;
% Abo=Ab/max(Ab);
% if exist('EDox')
%  Ab=abs(EDox).^2;
%  Abox=Ab/max(Ab);
% else
%  Abox=Ab*0;
% end
% figure, plot(xvero,Abc,xvero,Abox,xvero,Abo)
% AD=abs(ADo);
% ADo=AD/max(AD);
% AD=abs(ADc);
% ADc=AD/max(AD);
% figure, plot(KK,abs(ADc),KK,abs(ADo))
%
% figure,
% Ab=abs(EDc).^2;
% Ab=Ab/max(Ab);
% plot(xvero,Dvn,Ab,xvero,,xvero,Dvn.*Ab,'--')
% title(' giallo: Dvn,  magenta: campo, ciano: Dvn*campo')
%
% figure,
% Ab=abs(Evi).^2;
% Ab=Ab/max(Ab);
% plot(xvero,Dvn,xvero,Ab,xvero,Dvn.*Ab,'--')
% title(' giallo: Dvn,  magenta: campo, ciano: Dvn*campo')
