function gen_berlin(Pi,iint)

clear global
global tolk kretvero farm k_ret ikr iredmat nucomp ispeed ireset_int

global ifmonitor 

global lamm uL Dlam0 demm i1D
global la1Ds ga1Ds 
global la1Dr ga1Dr 

global isa_file nfile_sa 

%global if_only_out
%if_only_out=1;

ifmonitor=1;
ipvo=0;
global irid_bas
irid_bas=1;
if nargin==1
 iint=1;
end

global del_n_ag ianti_gui Nref
del_n_ag=-.1;

ifp=       Pi.ifp;

if iint==1
 irese0=1;
else
 irese0=-1;
end

fi_rad=     Pi.fin0;
fi_sa=     Pi.fi_sa;


%' fsop', keyboard
nfile_sa=fi_sa;

ee=exist([fi_sa,'.mat']);
ees=ee;
if ee==2
 Pss=Pi;
 save sa Pss
 fi_sas=fi_sa;
 eval(['load ',fi_sa]);
 clear P Pss
 load sa
 uL=Ppol.uL;
 fi_sa=fi_sas;
 D1ee=D1;
% if ifp==-10
%  restart
% else 
  iRestart=0;
% end
 if iRestart~=0
  if iRestart<0
   D1ee(abs(iRestart))=NaN;
  else
   D1ee(iRestart:end)=NaN;
  end
 else 
  if ifp==-10
%  ' quit ',
% keyboard
% keyboard
  end
 end 
 ee=2;
 ees=2;
 P=Pi;
D1ee
 if ifp==-10
% ' riparto'
% keyboard
 end
else
 P=Pi;
end
%' ifp', keyboard
doxv= P.dox;
dox= doxv(1);
%load data
if isfield(P,'pasnu')==1 
 global pasnu
 pasnu=P.pasnu;
end
%np_mir

Ps=P.Ps;
Pf=P.Pf;
if isfield(Pf,'salvasempre')==1
    isalvasempre=Pf.salvasempre;
else    
    isalvasempre=1;
end    

if isfield(Pf,'ianti_gui')==0
 ianti_gui=1; 
else
 ianti_gui=Pf.ianti_gui;
end
Ps.fi_sa=fi_sa;
if isfield(Pf,'isalvafile')==1
 isa_file=Pf.isalvafile;
else 
 isa_file=0;
end

iraff=0;
if isfield(Ps,'iraff')==1
 iraff=Ps.iraff;
end

global iKexact
iKexact=0;
if isfield(Ps,'iKexact')==1
 iKexact=Ps.iKexact;
end

istop_fie=0;
if isfield(Pf,'istop_fie')==1
 istop_fie=Pf.istop_fie;
end

if isfield(Ps,'ired_ret')==1
 ired_ret=Ps.ired_ret; 
else
 ired_ret=0;
end

global isoga GMA fiload iBEL lamm Dlam_mod istop_fie Pf   ired_ret



isoga=Ps.isoga;
i1D=Ps.i1D;
%'GNA', keyboard
if isfield(Pf,'GMA')==1
 GMA=Pf.GMA;
else
 GMA=1e6;
end
if isfield(Ps,'fiload')==1
 fiload=Ps.fiload;
else
 fiload='sol0';
end
if isfield(Ps,'iBEL')==1
 iBEL=Ps.iBEL;
else
 iBEL=0;
end
global isavetutto
if isfield(Ps,'isavetutto')==1
 isavetutto=Ps.isavetutto
else
 isavetutto=0;
end

global igr_app
if isfield(Ps,'igr_app')==1
 igr_app=Ps.igr_app;
else
 igr_app=0;
end

if isfield(Ps,'igr_app')==1
 igr_app=Ps.igr_app;
else
 igr_app=0;
end

%'Ps grapp', keyboard

if isfield(Ps,'veac')==1
 veac=Ps.veac;
else
 veac=0;
end

if isfield(Ps,'Rc_acc')==1
 Rc_acc=Ps.Rc_acc;
else
 Rc_acc=0;
end


Par=P.Ps.Par;


%' Par', keyboard
close all
 co(1:2,1)='w-';
 co(1:2,2)='m-';
 co(1:2,3)='c-';
 co(1:2,4)='g-';
 co(1:2,5)='r-';
 co(1:2,6)='b-';
 co(1:3,7)='y--';
 co(1:3,8)='m--';
 co(1:3,9)='c--';
 co(1:3,10)='g--';
 co(1:3,11)='r--';
 co(1:3,12)='b--';

 co1(1:2,1)='w.';
 co1(1:2,2)='m.';
 co1(1:2,3)='c.';
 co1(1:2,4)='g.';
 co1(1:2,5)='r.';
 co1(1:2,6)='b.';
%
r_pil=3;
ivfre0=0;


i2D=3;
mmvet=[0 1];


nomes='dp';





if isfield(Ps,'Npa1')==1
 Npa1=Ps.Npa1;
else
 Npa1=23;
end


if isfield(P,'farm')==0
 farm=1;
else
 farm=      P.farm;
end
if isfield(P,'att')==0
 del_att=0;
else
 del_att=      P.att;
end
ipro0=      P.ipro;
ipro=      P.ipro;
iLP=       P.iLP;
fi_sa=     P.fi_sa;
ifp=       P.ifp;

%alim0=      P.alim;
%alim=      alim0;
%nK_dis=    P.nK_dis;

alimv=      P.alim;
alim0=      alimv(1);
dkv=        P.dk;
dk=    dkv(1);
nK_dis=    P.nK_dis;
if length(alimv)>1
 alim1=      alimv(2);
 if length(dkv)>1
  dk1=    dkv(2);
 else
  dk1=    dkv(1);
 end
end
ikr=       P.ikr;
Ps.KKmax=alim1;
[lam0,rr]=lay_par(P.fil_str);
%'cont', keyboard

%lam0=.8836

numodiacc= P.numodiacc;
tolk=      P.tolk;
%dk=        P.dk;
ipro=      P.ipro;

%'quin contro', keyboard
%Dlam_mod=  P.Dlam_mod;
if isfield(P,'FRint')==1
 FRint=     P.FRint;
else
 FRint=1;
end
if isfield(P,'shfom')==1
 shfom=     P.shfom;
else
 shfom=0.5;
end
if isfield(P,'FRint0')==1
 FRint0=     P.FRint0;
else
 FRint0=0;
end
if isfield(P,'Del_su')==1
 Del_su = P.Del_su;
else
 Del_su = 1.5;
end
NPfreq=    P.NPfreq;
Ev_or_Od=  P.Ev_or_Od;
Dreti=     P.Dreti;
Dreli=     P.Dreli;
nucompi=   P.nucompi;
ikr=       P.ikr;
iany=      P.iany;
if isfield(P,'ianys')==1
 ianys=      P.ianys;
else
 ianys=0;
end 
dsub=      P.dsub;
tistru=    P.tistru;
DC=        P.DC;
period=    P.period;
dpm=       P.dpm     ;
dpm1=      P.dpm1    ;
fil_str=   P.fil_str ;
np_mir=    P.np_mir ;
if isfield(P,'np_miru')==1
 np_miru = P.np_miru;
end
if isfield(P,'H')==1
 H = P.H;
end
if isfield(P,'r')==1
 r_lens = P.r;
end
if isfield(P,'Lbuf')==1
 Lbuf = P.Lbuf*1000;
end
if isfield(Par,'ri')==1
 if Par.ri<0
  Par.ri=abs(Par.ri);
  iBase = 1;
  Lbuf = P.Lbuf;
 else
  iBase=0;
 end
end


Lm=mean(Lbuf);
if isfield(P,'dmet')==1
 Dmet = P.dmet;
else
 Dmet=0;
end
[lam0,rr]=lay_par(fil_str);
%'rr', keyboard
%lam0=.8836

Dlam_mod(1)=FRint0(1);   %in nm
Dlam_mod(2)=FRint(1);   %in nm
Dlam_mod(3)=NPfreq;
Dlam_mod(4)=shfom;   %in nm
if isfield(Pf,'stepF')==1
 Dlam_mod(5)=Pf.stepF;   %in nm
else
 Dlam_mod(5)=0;   %in nm
end
Dlam0=Dlam_mod;
if ee==2
% fi=find(delFre~=0);
% fi=find(delFre~=0 & delFre~=-1 & D1ee;
%'rifo fuori inter.', keyboard

% if length(fi)>0
%  D1ee(fi)=NaN;
%  iBEL=-1;
%  Dnew=[-.1 .1 3];
% end 
%'rifo fuori inter.', keyboard

end 
if exist('Dnew')
 Dlam0(1:3)=Dnew;
 Dlam_mod(1:3)=Dnew;
end
%'Dlam', keyboard

%iany=1;  % 1: anisotropia planare, 2: confinata (effetto elettro-ottico )
%iany=0;  % 1: anisotropia planare, 2: confinata (effetto elettro-ottico )
%ianys=0;   % strain
global fattany; fattany=1;
if iLP==1
 iany=0;
 ianys=0;
end

idyn=0;
nmasce=[30 2];
%nmasce=[1 1];

if pasnu==1
 nmasce=[4];
end 

if isfield(Pf,'nmasce')==1
 nmasce=Pf.nmasce;
end
%'cont Nmasce', keyboard
% nmasce=[6];



nucomp=nucompi;

%fakri=tolk;
%fakrs=tolk/2;
%fakrs=tolk;

fakri=tolk(1);
fakrs=tolk/2;
fakrs=tolk(1);
if length(tolk)>1
 fashif=tolk(2);
else
 fashif=0;
end


ispeed=1;
%'ispeed', keyboard
ispeed=0;
if isfield(Ps,'fast')==1
 ispeed=Ps.fast;
end
iredmat=1;
if igr_app==1
 iredmat=0;
end

global iNOsav
iNOsav=0;
if isfield(Ps,'iNOsav')==1 
 iNOsav=Ps.iNOsav; 
end

ipolar=P.ipolar;




%% ifp:  -10 alcuni stop campi
%% ifp:  -5  no display campi e no stop
%% ifp:  -4  display campi e no stop

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% calcolo e plot modi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%xroI1=linspace(0,3,30);
%xroI2=linspace(3,6,10);
%xroI=[xroI1 xroI2(2:length(xroI2))];

ploma=10*r_pil;
ploma=r_lens;

if isfield(Ps,'rplot')==1
 ploma=Ps.rplot;
end
Np=100;
npfi0=25;
Np=500;
Np=400;
%npfi0=45;
r=linspace(0,ploma,Np);
xroI=r;
%xroI=Np;
xro=r;

%' ploma ', keyboard

global igamveb igamveu GGbext GGuext
global ilo

igamveb=1;   %=1, gamma vero,  =0, riflessione GGext
igamveu=1;
GGbext=1;
GGuext=1;


lfi_inp=2*npfi0;
fimaxi=pi;

lfi_inp=4*npfi0;
fimaxi=2*pi;

if isfield(Pf,'adapt')==1
 FRadapt=Pf.adapt;
else
 FRadapt=0;
end
if isfield(Pf,'iplan')==1
 iplan=Pf.iplan;
else
 iplan=0;
end


fileeps='Gres_jo2';
Dla=-.005;


disp(' call caloptg '),


% variation management (2 nested loops)
% if centers are varying, put a negative number in the corresponding place:
% its absolute value will be the column of the matrix mat_cen where the centers
% must be put.

%Grating D=1|P1 d=.3|P2 shape=square|P5 Ry=3|P3 Rx=3|P4 orien.=0 shift=0 circle=1|P6

scas(1,:)='Case 1';
scas(2,:)='Case 2';
scas(3,:)='Case 3';
scas(4,:)='Case 4';
scas(5,:)='Case 5';
scas(6,:)='Case 6';

dech_g=30;
dech_r=10;
doxa=4;
%%DC=0.65;
%%period=0.7;
drel=doxa-1.25;
%Del_su=1.5;
%drsu=doxa*Del_su;
drsu=Del_su;

if isfield(Ps,'iacc_Ndis')==1
 iacc_Ndis= Ps.iacc_Ndis;
 Ndisc_ac = Ps.Ndisc_ac;
else
 iacc_Ndis=0;
end
if isfield(P,'Ndisc')==1
 Ndisc = P.Ndisc;
end

if isfield(Ps,'dar1')==1
 dar1=Ps.dar1;
else
 dar1=0;
end



if isfield(Ps,'BMc')==1
 BMc=Ps.BMc;
else
 BMc=0.5;
end
 BMc=100;


if isfield(Ps,'dar2')==1
 dar2=Ps.dar2;
else
 dar2=0;
end
if isfield(Ps,'dpost_e')==1
 dpost_e=Ps.dpost_e;
else
 dpost_e=0;
end
if isfield(Ps,'dpost_ee')==1
 dpost_ee=Ps.dpost_ee;
else
 dpost_ee=0;
end
if isfield(Ps,'dpost_i')==1
 dpost_i=Ps.dpost_i;
else
 dpost_i=0;
end
if isfield(Ps,'dcontact')==1
 dcontact=Ps.dcontact;
else
 dcontact=0;
end
if isfield(Ps,'Rc')==1
 Rc=Ps.Rc;
else
 Rc=2500;
end

global set_perd0
if isfield(Ps,'set_perd0')==1
 set_perd0=Ps.set_perd0;
else
 set_perd0=0;
end


if isfield(Ps,'N')==1
 N=Ps.N;
 ro=Ps.Nx;
else
 N=0;
 ro=0;
end
T=0;
zeta=0;

sha_ox=Ps.sha;
dr_ox2=Ps.r_ox2;
p_ox2=Ps.p_ox2;
pD1=Ps.D1;
pD2=Ps.D2;
pD3=Ps.D3;
rind=Ps.rind;
Rel=Par.ri(1);

par_ind{1,1}= H;            % H lente
par_ind{2,1}= r_lens(1);       % raggio lente
par_ind{3,1}= Ndisc;        % N discr
par_ind{4,1}= Par.ri(1);            % Larg. lente
par_ind{5,1}= np_mir;       % Paia specchio
par_ind{6,1}= Lbuf(1);          % spacer
par_ind{7,1}= P.att;       % forma
par_ind{8,1}= dox;          % attivo
par_ind{9,1}=Par.re(1);            % altro asse
par_ind{10,1}=p_ox2;           % Delta 2
par_ind{11,1}=P.att;           % Delta 2
par_ind{12,1}=dpost_e;           % Specchio sotto raggio esterno
par_ind{13,1}=dpost_i;           % Specchio sotto raggio interno
par_ind{14,1}=dcontact;          % Contatto
par_ind{15,1}=Npa1(1);           % Specchio sotto raggio esterno oro
par_ind{16,1}=0;
par_ind{17,1}=Par.NpaDie(1);
Datt=dox+P.att;
%' datt', keyboard
LS0=15;
LS0=22;
if isfield(Ps,'delmir')==1
 delmir=Ps.delmir;
else
 delmir=0;
end

LS0=dox*2+delmir;
%'cont', keyboard



par_ind{18,1}=dar1;
par_ind{19,1}=dar2;

Part=Ps.Part;
par_ind{21,1}=Part.Rp(1);
par_ind{22,1}=Part.cp(1);
par_ind{23,1}=Part.Np(1);

if isfield(Par,'oxide')==1
 Refr=Par.oxide;
 par_ind{16,1}=Refr(1);
else
 Refr=1.775;
end 

if isfield(Par,'lossQW')==1
 glass=Par.lossQW;
 Par.glass=glass;
 par_ind{26,1}=glass;
else
 glass=0;
end 

%'Refr'
%keyboard

%Dp=4|P21 Cp=2.8|P22  Ndis=20|P23 

par_ind{26,1}=Par.glass(1);
par_ind{27,1}=Par.air(1);

if isfield(Ps,'ifcav')==0
 ifcav=0;
else
 ifcav=      Ps.ifcav;
 if ifcav==1
  Detch=Ps.Detch;
  D_next=Ps.D_next;
  Dre1=Ps.Dre1;
 end
end


if isfield(Par,'arc')==1
 par_ind{30,1}=Par.arc(1);
 par_ind{29,1}=Par.primarc(1);
else
 par_ind{29,1}=0;
 par_ind{30,1}=0;
end
%  Par.thonret=P0-dret;
 par_ind{31,1}=10000;


if isfield(Par,'thonret')==1
 thonret=Par.thonret;
 par_ind{41,1}=thonret(1);
end
if isfield(Par,'thbelret')==1
 thbelret=Par.thbelret;
 par_ind{40,1}=thbelret(1);
end

%'Par', keyboard

%' parametri'
%keyboard

par_map{1,1}='H specchi ';
par_map{2,1}='raggio specchio ';
par_map{3,1}='Discret. specchio ';
par_map{4,1}='Rilievo ';
par_map{5,1}='Paia specchio ';
par_map{6,1}='Spacer ';
par_map{7,1}='Dimensione attiva ';
par_map{8,1}='Dimensione apertura';
par_map{9,1}='Apertura metallo';
par_map{10,1}='Delta_2';
par_map{11,1}='Arc';
par_map{12,1}='Dim est spec sotto';
par_map{13,1}='Dim int spec sotto';
par_map{14,1}='Dim Contatto';
par_map{15,1}='Paia out half VCSEL';
par_map{16,1}='Ref Oxide';
par_map{17,1}='Np Diel GaAs';
par_map{18,1}='Arc1';
par_map{19,1}='Arc2';
par_map{21,1}='Raggio particle';
par_map{22,1}='Centro particle';
par_map{23,1}='Ndis particle';
par_map{25,1}='Curvilineo sotto';
par_map{26,1}='Perdite QW';
par_map{27,1}='Dente Ox HCG';
par_map{29,1}='Spes prima GaAs';
par_map{30,1}='GaAs spessore';
par_map{31,1}='TOPAS spessore';


if isfield(Par,'grape')==1
  LA=Par.grape;
  DC=Par.gradc;
  Dret=Par.dret;
  gr_acc=Ps.gr_acc;
 par_ind{31,1}=Par.dret(1);
 par_ind{32,1}=LA(1);
 par_ind{33,1}=DC(1);
 par_map{31,1}='RET spessore';
 par_map{32,1}='RET period';
 par_map{33,1}='RET DC';
 end





%global iphase_front  isav_nfcut
global iphase_front  isav_nfcut Ps

iphase_front=Ps.iphase_front;
isav_nfcut=Ps.isav_nfcut;

for kp=1:length(par_ind)
kps{1}=kp;
disp([kps par_ind{kp} par_map{kp}])
end


if ifp==-10
 pausak
end



lvet=[length(dsub) length(dox) length(period)];
global Lmax irec_fie iord_long num_long
Lmax=151000;

if isfield(Ps,'irecfie')==1
 irec_fie=Ps.irecfie;
else
 irec_fie=0;
end
%'irec;,', keyboard
if isfield(Ps,'iord_long')==1
 iord_long=Ps.iord_long;
 num_long=1:iord_long;
else
 iord_long=0;
 num_long=1;
end




data1=Lbuf;
data2=doxv;
data3=r_lens;
data4=Ndisc;
data5=Rc;
data6=BMc;
data7=Par.h;
data8=Par.ri;
data9=Par.re;
data10=Par.Conf;
data11=Par.arc;

if isfield(Ps,'gr_acc')==1
 iacre=Ps.gr_acc;
else
 iacre=0;
end 
iacred=1;

  NPAR=21;
  data12=Dret;
  data13=Par.grape;
  data14=Par.gradc;
  data15=Par.air;
  data16=glass;
  data17=Refr;


if isfield(Par,'misal')
 data18=Par.misal;
else
 data18=0;
end
data19=iord_long;
data20=Par.NpaDie;
data21=Par.NpaBack;



NPAR_ADD=4;
  

 De0=data11(1);

%nloop1=length(data1);
for k=1:NPAR
 kst=num2str(k);
 eval(['nloop',kst,'=length(data',kst,');']);
end


nve=[];
for k=1:NPAR
 kst=num2str(k);
 eval(['nloop',kst,'=length(data',kst,');']);
 eval(['nve=[nve nloop',kst,'];']);
end
[dupar,iso]=sort(nve);

%keyboard

if Rc_acc~=2

np1=iso(end);
np2=iso(end-1);
np3=iso(end-2);

npa1=dupar(end);
npa2=dupar(end-1);
npa3=dupar(end-2);

npacc=[];
else

np3=iso(end);
np2=iso(end-1);
np1=iso(end-2);

npacc=[np2 np3];

npa1=dupar(end-2);
npa2=1;
npa3=1;

end

eval(['nl1=nloop',num2str(np1),';']);
eval(['nl2=nloop',num2str(np2),';']);
siFr=size(FRint);

 iassFR=0;
 pF=prod(siFr);
if pF>1
 iassFR=1;

 if pF~=nl1*nl2
  'errore Fint ', keyboard
  'errore Fint ', keyboard
  'errore Fint ', keyboard
 end
end
%'prima cont ee', keyboard



if ee~=2
%if ee~=3  %ricalcola i parametri
 D1=zeros(nl1,nl2)*NaN;
 for k=1:NPAR
  kst=num2str(k);
  eval(['D',kst,'=D1;']);
  eval(['loind',kst,'=1:nloop',kst,';']);
 end
else
 for k=1:NPAR
  kst=num2str(k);
  eval(['D',kst,'=D',kst,'(1:nl1,1:nl2);']);
  eval(['loind',kst,'=1:nloop',kst,';']);
 end
end

%' cont ee', keyboard




t00=clock;
%' par_ind', keyboard

par_in=par_ind;

Dan=0;
iag=0;
if iag==1
%' ui', keyboard
  if ianti_gui~=0
    Dattg=Datt;
    Dattn=Datt+Dan;
    ro_in=sort([Dattg Dattn]);
    if Dattg<Dattn
     N=[1 0];
     Nref=[0 1];
    elseif Dattg==Dattn
     ro_in=[Dattg Dattg+1e-3];
     N=[1 0];
     Nref=[1 0];
    elseif Dattg>Dattn
     N=[0 1];
     Nref=[1 0];
    end
%    'Datt', keyboard
  else
   N=0;
   ro_in=ro;
  end
end


% vario alim
alim=alim1;
    nk=round(alim/dk);
    nK_dis=nk;

%'nK', keyboard

%np1=iso(end);
%np2=iso(end-1);
%np3=iso(end-2);

%npa1=dupar(end);
%npa2=dupar(end-1);
%npa3=dupar(end-2);
NPA=0;
if npa2>1
 NPA=NPA+1;
end 
if npa3>1
 NPA=NPA+1;
end 
npaver=iso(end-NPA:end);

puind=ones(1,NPAR);
	 for k=1:NPAR
	  kst=num2str(k);
	  eval(['Dati(',kst,')=data',kst,'(1);']);
	 end 

%  addata1=Lbufa;
%  addata12=Part.Rpa;
%  addata15=Par.aira;
%  addata16=Par.glassa;
%  addata12=Part.Rpa;

% controllo somme: 
% posizione: parametro vecchio
% 0: no somma; # parametro da sommare
%            1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20  ...
%Isomma=[

if isfield(Pf,'nKK')==1
 nKK=Pf.nKK;
 npKappa=length(nKK);
 alve=Pf.alve;
  if length(alve)~=npKappa
   'errore nKK in gen_wolu760', keyboard
  end
else 
 npKappa=1;
end

%'keyboard nPkappa', keyboard


for k4=1:npKappa
 if npKappa>1
  alim=alve(k4);
  nk=nKK(k4);
  dk=alim/nk;
  nK_dis=nk;
 end
global iORTA
iORTA=1;

for k3=1:npa3
 puind(np3)=k3; 
 for k2=1:npa2
  puind(np2)=k2; 
  for k1=1:npa1
%  for k1=3
   puind(np1)=k1; 
	 if length(npacc)>0
	   puind(np2)=k1; 
	   puind(np3)=k1; 
	 end
	 for k=npaver
	  kst=num2str(k);
	  kpar=num2str(puind(k));
	  eval(['Dati(',kst,')=data',kst,'(',kpar,');']);
	 end    

	 for k=1:NPAR
	  kst=num2str(k);
	  eval(['df',kst,'=Dati(',kst,');']);
	 end 
  if isfield(Ps,'iprecv')==1
   iprecv=Ps.iprecv;
   numodiacc=abs(iprecv(puind(18)));
  end
  NpaDie=df20;
  par_in{17,1}=NpaDie;
  NpaBack=df21;
  par_in{15,1}=NpaBack;  
  iord_long=df19;
  global cce
  cce=df18;
  if cce==0
   cce=1e-3;
  end
  Ref_sfe=df17;
  par_in{16,1}=Ref_sfe;
  glas=df16; 
  par_in{26,1}=glas;
  air=df15; 
  par_in{27,1}=air;
  if iacred==1
   DC=df14;
  else 
   Np=df14;
   par_in{23,1}=Np;
  end 
  
  if iacred==1
  
   LA=df13;
   par_in{33,1}=LA*DC;
   par_in{32,1}=LA;
%   lam0=.883

   if ikr>=1
    k_ret=lam0/(rr*LA)-fakri;
    kretvero=lam0/(rr*LA);
   end
	alim=alim0;	
	    if ikr>1
	     k_ret=lam0/(rr*LA)-fakri;
	     kretvero=lam0/(rr*LA);
	     if alim0<=0
	      alim=farm*lam0/(rr*LA)+fakrs+fashif;
	     else
	      k_ret=alim0;
	      kretvero=alim0;
	      alim=alim0;
	     end
	     
	     if ikr==3
	      alim=2*lam0/(rr*LA)+fakrs+fashif;
	     end
	    else
	     alim=alim0;
	    end
	   if exist('alim1')
	    if alim1>k_ret
	     alim1=k_ret;
	    end
	    nk=fix((alim-alim1)/dk);
	    nK_dis=nk;
	    alims=alim;
	    nk1=fix(alim1/dk1);
	    nK_dis(1)=nk1;
	    nK_dis(2)=nk;
	    alim(1)=0;
	    alim(2)=alim1;
	    alim(3)=alims;
	  end
	    nk=fix(alim/dk);     

%' qui alim', keyboard


 
  %'akr',keyboard
 else
  cp=df13;
  par_in{22,1}=cp;
 end %iacre

%'ikr', keyboard

 Darc=df11;
 De0=Darc;
 par_in{30,1}=De0;
 
  if iacre==1
   dret=df12;
   par_in{31,1}=dret;
   thonret=thonret0-dret;
   if thonret<0
    'parametro sbagliato, reticolo troppo spesso'
    keyboard
   end
  else
   par_in{21,1}=df12;
  end  %iacre


  Rconf=df10;
  Rext=df9;
  Rint=df8;
  Base=df8;
  Lbuf=df1;
  Rc=df5;
  r_lens=df3;
  if  iBase==1 & Rc==0
   Base=Rint;
   l_minore=Base-r_lens;
   angolo=Lbuf/l_minore;
   Rc=r_lens*sqrt(1+1/angolo^2)
  else
  L2=(Lbuf-Rc)^2+Base^2;
  S2=L2-Rc^2;
  roa=asin(Rc/sqrt(L2));
  aa=asin((Lbuf-Rc)/sqrt(L2));
  abas=roa+aa;
  Tb=sqrt(S2)*sin(pi/2-abas);
  r_lens=Base-Tb;
  df3=r_lens;
  %'sono qui', keyboard
  end
  par_in{4,1}=Rint;
  Hre=df7;
  BMc=df6;
  pppp=LS0-2*BMc;
  Ndisc=df4;
  if length(loind3)>1
   if iacc_Ndis==1
    Ndisc=Ndisc_ac(loop3);
   end
  end

  Dox=df2;
  Dat=Dox+del_att;
  Lbuf0=df1;

  if iacc_Ndis==1
   Ndisc=Ndisc_ac(loop5);
  end
%  if Rc_acc==1
%   Rc=Rc*Lm/Lbuf;
%   df5=Rc;
%  end


   Lbuf=Lbuf0;
   if Ndisc==1
    Nmiru=np_mir;
%    Nmiru=abs(Nmiru);
    H=0;
   else
    if Rc>0
     H=(Rc-sqrt(Rc^2-r_lens^2));
    else
     H=0;
    end
    if imag(H)~=0
     'Errore parametri lente in GEN_voluB', keyboard
    end
    if iBase==0
     Lbuf=Lbuf0-H*1e3;
    else
     Lbuf=Lbuf0-H;
    end
    Nmiru=np_mir;
   end
   par_in{5,1}= Nmiru;       % Paia specchio
%'qui Rc', keyboard
   par_in{1,1}= H;            % H lente
   par_in{6,1}=Lbuf;
   par_in{8,1}=Dox;
   par_in{7,1}=Dat;
   par_in{2,1}=r_lens;
   par_in{3,1}=Ndisc;
   par_in{4,1}=Rel;
   par_in{9,1}=Rext;
   par_in{10,1}=Rconf;
   par_in{20,1}=Base;
   par_in{11,1}=Darc;
  dret=df12;
  par_in{31,1}=dret;   

%keyboard

  global lp1 lp2 lp3 lp1c
  lp1=k1;
  if npa1==1
   lp1=k4;
  end
  lp2=k2;
  lp3=k3;
  lp1c=lp1;
%'qui 760', keyboard
%  if length(find(np1==[6 15 16 17 13]))>0
   lp1c=1;
%  end
if ifp==-10  
'  [lp1 lp2 lp3 lp1c]'
  [lp1 lp2 lp3 lp1c]
%  pausak
end
    


  if ee==2
   if isnan(D1ee(lp1,lp2))==0
    ipro=-1
   else
    ipro=ipro0;
   end
  else
   ipro=ipro0;
  end
  differ=1;
%' differ', keyboard
  if ipro==1 & differ>=0
   if k1==1
   [par_in par_map]
   end
  [lp1 lp2 ]
%   eval(['ddd1=data',num2str(np1),'(',num2str(lp1),');']);
%   eval(['ddd2=data',num2str(np2),'(',num2str(lp2),');']);
Dati(npaver)
%   pausak
  elseif ipro==0
   [par_in par_map]
  end

if iassFR==1  
  eval(['FR1=FRint0(loop',num2str(np1),',loop',num2str(np2),')'])
  eval(['FR2=FRint(loop',num2str(np1),',loop',num2str(np2),')'])
Dlam_mod(1)=FR1;   %in nm
Dlam_mod(2)=FR2;   %in nm
Dlam_mod(3)=NPfreq;
Dlam_mod(4)=shfom;   %in nm
end
  if  ipro==1
%   ' numero tot parametri = ', npa1*npa2
%   ' tempo stimato di calcolo (h) = ', npa1*npa2*Dlam_mod(3)/40
   ' controllo parametri', pausak
  end
%   par_in=[];
   t0=clock;
  if ipro==0 & differ>=0


  if ifp==-10
   ' gen_wol VEL', keyboard
  end

%   keyboard
%   ' gen_mes VEL', keyboard
   if isfield(Pf,'iORTA')==1
    if Pf.iORTA==1
     iORTA=1;
    end
   end
%   keyboard
   [Eqw,Eout,xro,fian,lambda,delf,gain,ord,nrAzim,Cu,PaDy,ADom,ADcm,Plot,Ppol,Dla_out,imod_err]=...
   vel(N,T,ro,zeta,xroI,fimaxi,lfi_inp,par_in,r_pil,fil_str,...
   ipolar,Ev_or_Od,nmasce,numodiacc,Dlam_mod,nK_dis,alim,...
   iplan,iraff,idyn,i2D,iLP,iany,ianys,ifp,Dla,fileeps,nomes,ivfre0,1,[],veac);

   if length(Dla_out)>0
    Dlam_mod=Dla_out;
   end
   if FRadapt~=0
    if iassFR==1, 
     'errore gestione Frequenza', keyboard
    end
    if ifp==-10
    ' aggiusto qui', keyboard
    ' aggiusto qui', keyboard
    end
    Dcomp=delf(1)*1000*lambda;
    Dlam_mod(1:2)=Dcomp+[-1 1]*FRadapt;
    Dlam_mod(3)=Pf.adaN;
    
   end


   if imod_err==0
    for k=1:NPAR
     kst=num2str(k);
     eval(['D',kst,'(lp1,lp2)=df',kst,';']);
    end

  if length(nKK)>1
   Dp=1:length(nKK);
  else
    eval(['Dp1=D',num2str(np1),';']);
    eval(['Dp2=D',num2str(np2),';']);
    eval(['Dp3=D',num2str(np3),';']);
     Dp=Dp1;
    if length(Dp2)>1
    if diff(Dp2(1:2))~=0
     Dp=Dp2;
    end    
    end
    if length(Dp3)>1    
    if diff(Dp3(1:2))~=0
     Dp=Dp3;
    end    
    end
    if Rc_acc==2
     Dp=Dp(:,1);
    end
 end

%    if np1==1
%     Dp=D1;
%    elseif np1==2
%     Dp=D2;
%    elseif np1==3
%     Dp=D3;
%    elseif np1==4
%     Dp=D4;
%    elseif np1==5
%     Dp=D5;
%    elseif np1==6
%     Dp=D6;
%    end

   global  lastimm ilargeint del_Frea gastimm
   
%   ' dopo Vel', keyboard
   if i1D==0
    DlamV{lp1,lp2,lp3}=Dlam_mod;
%    delFre(lp1,lp2,lp3)=del_Frea;
    lamm{lp1,lp2,lp3}=lambda./(1+delf);
    demm{lp1,lp2,lp3}=delf;
    pamm{lp1,lp2,lp3}=ones(size(delf))*lp1;
    gamm{lp1,lp2,lp3}=gain;
    tymm{lp1,lp2,lp3}=PaDy.l';
    M2mm{lp1,lp2,lp3}=PaDy.M2;
    ormm{lp1,lp2,lp3}=ord;
    fatqw=Ppol.fatqw;
    fqwmm{lp1,lp2,lp3}=fatqw;
%    'fatwe'
%    keyboard

    Pdup=Plot;
    Pdup.YP=Pdup.YP{1};
    Pdup.XP=Pdup.XP{1};
    Pdup.X=Pdup.X{1};
    Pdup.Y=Pdup.Y{1};
    Pdup.Cug=Pdup.Cug{1};
    Pdup.ldap=Ppol.ldap;

   
   if isalvasempre==0    
    clear Psav
    Psav.ldap=Pdup.ldap;
    Psav.Ap=Pdup.Ap;
    Psav.ApQ=Pdup.ApQ;
    Psav.parmod=Pdup.parmod;
    Psav.E2xo=Pdup.E2xo;
    Psav.E2yo=Pdup.E2yo;
    Psav.Ef=Pdup.Ef;
    Psav.Cug=Pdup.Cug;
    Psav.ze=Pdup.ze;
    Psav.gg0=Pdup.gg0;
    Psav.ldap=Pdup.ldap;
    if lp1==1 & lp2==1
        Psav.XP=Pdup.XP;
        Psav.YP=Pdup.YP;
        Psav.X=Pdup.X;
        Psav.Y=Pdup.Y;        
    end
    Pdup=Psav;
   end 
    PPlot{lp1,lp2,lp3}=Pdup;    
%    ' PPlot ', keyboard
    clear ee ipro P Pi Pss
%    eval(['save ',fi_sa])

   end
    La1D{lp1,lp2,lp3}=la1Ds;   
    Ga1D{lp1,lp2,lp3}=ga1Ds;   
    La1Dr{lp1,lp2,lp3}=la1Dr;   
    Ga1Dr{lp1,lp2,lp3}=ga1Dr;    
%'dopo', keyboard
varsa=' np1 demm gamm tymm M2mm Ppol lamm ormm Dp  Dp1 Dp2 iLP lambda D1 D2 D3 D4 D5 D6 D7 D8 D9 D10 D11 D12 D13 D14 D15 D16 D17 D18 D19 D20 D21';
%'cosa salvo', keyboard
    eval(['save ',fi_sa])
   ee=ees;
   end

  end
   PA_HCG=[11 12 13 14 20];
   if length(find(PA_HCG-np1==0))==0
    iORTA=0;
   else 
     iORTA=1;
   end
   if isfield(Pf,'iORTA')==1
    if Pf.iORTA==1
     iORTA=1;
    end
   end   
%  'dopo vel ', pausak


   tim=etime(clock,t0);
   ore=fix(tim/3600);
   minuti=fix((tim-ore*3600)/60);
   secondi=(tim-ore*3600-minuti*60);
   format short
   if ipro0==0
   disp(' Partial elapsed time (hours, minutes, seconds) = ');
   disp([ ore minuti secondi]);
   end
end   %fine k1
   if length(find(PA_HCG-np2==0))==0
    iORTA=0;
   else
    iORTA=1;
   end
   if isfield(Pf,'iORTA')==1   
    if Pf.iORTA==1
     iORTA=1;
    end
   end   
     iORTA
%  'dopo vel 2 ', pausak
%'qui prima di salvare', keyboard
end  %fine k2
   if length(find(PA_HCG-np3==0))==0
    iORTA=0;
   else
    iORTA=1;
   end
   if isfield(Pf,'iORTA')==1
    if Pf.iORTA==1
     iORTA=1;
    end
   end   
     iORTA
  %'dopo vel 3', pausak
end  %fine k3
 iORTA=1;
end  %fine k4

   tim=etime(clock,t00);
   ore=fix(tim/3600);
   minuti=fix((tim-ore*3600)/60);
   secondi=(tim-ore*3600-minuti*60);
   format short
   disp(' Total elapsed time (hours, minutes, seconds) = ');
   disp([ ore minuti secondi]);

[' Nome file dati: ',fi_sa]

global Dsav

if length(Dsav)>1
  eval(['cd ',Dsav]);
  eval(['!del *.mat']);
  eval(['cd ..']);
  eval(['!rd ',Dsav]);
end

