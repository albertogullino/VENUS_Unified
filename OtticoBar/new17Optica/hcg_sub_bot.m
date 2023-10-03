%1/8*(115.6/+130.59)+115.6=4547
close all
clear all
ipro=0   %prova loops
inew=0;

load sabard

%'eve', keyboard

irecfie=0
Ps.irecfie=irecfie;  %ricalcola campo

if ifiez==1 
Ps.isav_Az=1;       %0 no campo longitudinale, 1 si
Ps.Lizi=10/1000;    % step longitudinale in micron   %
else
 Ps.isav_Az=0;       %0 no campo longitudinale, 1 si
end
%'ifiez', keyboard


if isfield(Part,'Rp')==0
 Part.Rp=.02;
end

Part.Np=50;
%Part.Np=20;

Ps.iNOsav=0;  %0, salva le matrici, 1, no
Ps.irecfie=1;  %ricalcola campo

%Ps.iord_long=[2];


  Par.misal=Mis;
  Par.air=air;
  Par.glass=glass;
  load struttura
  if ~exist('modacc')
   modacc=0
  end
%  'modacc', keyboard
   Ps.iprecv=-ones(size(Mis))*modacc;
%   Ps.iprecv=-[1 1 2 3 4 5];
%   Ps.iprecv=-[0];



Ps.Npa1=Npa1;
Ps.Part=Part;

Ps.iret_fisso=1;
Ps.ired_ret=1;
Ps.ired_ret=2;
Ps.ired_ret=0;

Ps.iold_1=1;   %1 vecchio metodo chain


%ifp=1
   set_perd0=0;
%   str='bar4cnomet';

load struttura
%   str='tip';
   isho=1;
   isho=0;

%   str='bar40cmodmet';
%   isho=1;  % campo long numerico
   
%   str='barp';
%   str='bar0c';
   
%   str='bar4ccc';
%   str='bar4cccm0';
   
%   str='bar0cc1';
%Lca    25   60.93  Ref  3.476  3  circle 15 Ref 1.5 0
Ps.iord_long=[1];
if isho==1
 Ps.iord_long=[0];
end    

   
  dsub=[.77];
dret=Par.dret;
period=Par.grape;
DC=Par.gradc;
  
  
  iconf=1;
  iconfm=0;
  Deltae=2;
  du=0:.2462:.2462*10;
  du1=.1156+du;
  dut=sort([du du1]);
%  Par.grape=period;

  Par.grasiz=7;




  
  Par.h=dut;
  Par.h=[0];
%'Par', keyboard
  if isfield(Par,'re')==0
   Par.re=0;
  end
  Par.raxdivy=1-linspace(.005,.05,15);
  Par.raxdivy=1;
  dcontact=20*iconfm;
  dpost_e=[24]/2*iconfm;
  dpost_ee=[24+Deltae]/2*iconfm;
  dpost_i=7*iconf;

%'idar', keyboard
  dar1=43*idar;
  dar2=190*idar;

%  Rc=0;

  Ps.gr_acc=0;  % scala reticolo da strato sotto   



  Ndisc=-10
  if isfield(Part,'tip')==1
    Ndisc=Part.tip;
  end
%  Ndisc=-2
%  Ndisc=-20
  
  Nmiru=-1;
  Nmiru=0;

%  ra=30;
%  Ndisc=2;
%  Ndisc=260;
%  Ndisc=1;      %specchio piano
%  ra=40;
%  Ndisc=60;
%  if Ndisc==1
%   Nmiru=abs(Nmiru);
%   H=0;
%  else
%   H=Rc-sqrt(Rc^2-ra^2);
%   Lbuf=Lbuf-H
%  end
% pausak
   H=0;

iglass=0;

iLP=0;
%iLP=1;
iany=0;


dpm=0;

iphase_front=0;
isav_nfcut=0;
rind=1;
D1=0;
D2=0;
D3=0;
sha_ox='circle';
r_ox2=[];
p_ox2=[];


dec=0;
dpm=0;
dmet=0;



  drel=[0];


%' quia', keyboard
 if ~exist('ra')  
  if exist('addox')
   ra=dox+addox;
  else 
   ra=dox+10;
  end
 end 
%' quia dopo', keyboard
  Ps.delmir=2;
  Mtot=(dox+Ps.delmir/2);
  BMcurvo=[3];
  BMcur=Mtot-BMcurvo;





  sha_ox='rhombus';
  sha_ox='circle';
  am_sha=0;
  dp_sha=0;

  r_ox2=am_sha;
  p_ox2=dp_sha;

  amod=-1;
  amod=-2;
  amod=-3;
% amod=3;
% amod=-2;
if isfield(Pf,'ifund')==0
 Pf.ifund=1;
end

if Pf.ifund==0
 amod=-3;
 Ps.veac=-1;  %=-1 accopia TE - TM modo 1; 0 entrambi e disaccop, 1(2) solo 1 o 2
elseif Pf.ifund==1
 amod=-1;
 Ps.veac=-1;  %=-1 accopia TE - TM modo 1; 0 entrambi e disaccop, 1(2) solo 1 o 2
elseif Pf.ifund==2
 amod=-2;
 Ps.veac=-1;  %=-1 accopia TE - TM modo 1; 0 entrambi e disaccop, 1(2) solo 1 o 2
end
%Ps.veac=2;  %=-1 accopia TE - TM modo 1; 0 entrambi e disaccop, 1(2) solo 1 o 2
%Ps.veac=-1;  %=-1 accopia TE - TM modo 1; 0 entrambi e disaccop, 1(2) solo 1 o 2
if iLP==1
 Ps.veac=1;  % 0 Zmedio, 1 Zte, 2 Ztm
end


%'contr ', keyboard

  if isfield(Pf,'sh')==0
    Pf.sh=.2;
  end
  if iLP==1
  Pf.sh=2.5;
  end

 
  Pf.stepF=.02;
  Pf.stepF=.05;   % =0 opera normalmente, altrimenti 2 risonanze (in fsr_sub) e discretizza fra due sol long
  Pf.stepF=.0;   % =0 opera normalmente, altrimenti 2 risonanze (in fsr_sub) e discretizza fra due sol long
  Pf.stepF=8;   % =0 opera normalmente, altrimenti 2 risonanze (in fsr_sub) e discretizza fra due sol long
%  Pf.fi=2;
%  Pf.fu=6;
%  Pf.N=5;

  
%  Pf.fi=-7;
%  Pf.fu=3;
%  Pf.N=10;
%  if iconfm+iconf>0 & iLP==1
  if iconfm+iconf>0 & iLP==2
   Pf.sh=2.;
   Pf.fi=-4.;
   Pf.fu=-2;
  end


 na0=[1 2];
 
 %'qui', keyboard
 
 
 
load struttura 
ize=1;
if ize==1
  Ps.ikr=0; 
  dk0(1)=alim(2)/Nk;  
end

if Ps.igraef_new==0

end

  Ps.Nk=Nk;



   Nrel=0;
   fi_add=rad;

 Ps.set_perd0=set_perd0;
 Ps.dcontact=dcontact;
 Ps.dpost_i=dpost_i;
 Ps.dpost_e=dpost_e;
 Ps.dpost_ee=dpost_ee;
 Ps.BMc=BMcur;
 Ps.dar1=dar1;
 Ps.dar2=dar2;
 Ps.Rc=Rc;
 Ps.dat=Datt;
 Ps.sha=sha_ox;
 Ps.r_ox2=r_ox2;
 Ps.p_ox2=p_ox2;
 Ps.D1=D1*1000;
 Ps.D2=D2*1000;
 Ps.D3=D3*1000;
 Ps.rind=rind;
 Ps.Nrel=Nrel;
 Ps.Nmiru=Nmiru;
 Ps.iphase_front=iphase_front;
 Ps.isav_nfcut=isav_nfcut;
 Ps.Par=Par;

a=Datt(1);
ron=linspace(0,1.5,100);
xx=ron*a;
nv=linspace(.1,1.7,10);
nv=1.08;
ip=0;
for n=nv
ip=ip+1;
f1=3+n+n^2;
f2=2+n^5;
Id=exp(-(ron).^f1+n*ron.^f2);
Nd=2*pi*Id.*xx*[diff(xx) 0]';
I(ip,:)=Id/Nd;
In1(ip,:)=Id/max(Id);
end
gdist=In1;

%figure, plot(ron,In1), pausak

 Ps.N=gdist;
 Ps.Nx=xx;

 Ps.N=0;
 Ps.Nx=0;

 Ps.ifcav=0;
 
 
 tolk(1)=.1; % half of the interval of values for k-vector in the grating region
 tolk(2)=-.03; % offset with respect to lambda0/(period*n_r)
 
 
 if Ps.ikr==2
  dk0=alim/Nk;
  dkLF=dk0(2)*2;
  dk0(1)=dkLF;
 end
 
%load iprec
tolk(1)=.05; % half of the interval of values for k-vector in the grating region
tolk(2)=.0; % offset with respect to lambda0/(period*n_r)

%'tolk' , keyboard

 Ps.tolk=tolk;

%   alim=.05;
%   Nk=[30];
% 
%   dk0=alim./Nk;
%  Ps.Nk=Nk;
if dret>0 
 putype=[5]; 
 iprec=-5;  % normale
 %iprec=-2;  % normale
% iprec=-4;  % normale
 iprec=-2;  % normale
% iprec=0;  % normale
 ipolar=2;
else
 iprec=-1;  % normale
 iprec=0;  % normale
 putype=[0]; 
 ipolar=-1;
end
isoga=1;
Ps.isoga=isoga;
ipolar=1;       % quella OK
ipolar=-1;
%ipolar=2;


pasnu=2;
if isfield(Ps,'pasnu')==1
 pasnu=Ps.pasnu;
end
 
if isfield(Ps,'igr_app')==1
 igr_app=Ps.igr_app;
else
 igr_app=0;
end
 
if igr_app>0
 iprec=0;  % normale
% iprec=-1;  % normale
 putype=0;
 Ps.ired_ret=0;
end 
 if igr_app==-10
  igr_app=0;
  Ps.ired_ret=1;
 end

%ipolar=-1;

%Ps.isav_Az=1;       %0 no campo longitudinale, 1 si
%Ps.Lizi=40/1000;    % step longitudinale in micron   %

%' doxme Bard', pausak
 hcgbot(str,dret,drel,dox,dsub,period,DC,fi_add,putype,ipro,iprec,ifp,alim,iLP,...
 iany,dk0,dec,dpm,H,ra,Ndisc,Lbuf,dmet,Ps,na0,amod,Pf,pasnu,fun_wolu)
 

%' doxme', pausak

	