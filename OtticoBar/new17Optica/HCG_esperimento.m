iscan=1  % 1 esplorazione grossolana parametri
%iscan=0  % 0 calcola tutti gli ORTA
Ps.fast=1;   % velocizza matr. trasm. con lineariz.
Pf.iORTA=1;    %=0 calcola Orta solo quando serve
Pf.NModi=11;    %=0 calcola Orta solo quando serve
Pf.ICRIT=1;
Pf.ICRIT=0;

ifiez=0;   % campo in z ovunque per fz_asim
Mis=0;
Mis=1000;
ifp=-4;
%ifp=-10;
%ifp=-11;
ifpstop=0;
ifpstop=1;
%ifp=0;
ipro=1;
%ipro=2;
ipro=0;


fiload='filoP';
iForm=0;      %=1 forma punta
inotip=1;    % 1, polimero;   2, aria
inotip=2;    % 1, polimero;   2, aria
inotip=0;     %c'e' punta

modacc=3;   %modi accoppiati
%modacc=1;   %modi accoppiati
%modacc=5;   %modi accoppiati
modacc=1;   %modi accoppiati
%modacc=2;   %modi accoppiati
modacc=0;   %modi accoppiati

pasnu=2;
Ps.pasnu=pasnu;
Ps.tuttefreq=0;       %calcola Orta sempre
Ps.tuttefreq=1;       %calcola Orta sempre
Ps.iraff=1;    % =1 ricalcola soluzione vicino alla stima
Ps.iraff=0;    % =1 ricalcola soluzione vicino alla stima


if iscan==1
Ps.tuttefreq=0;       %calcola Orta sempre
Ps.iraff=0;    % =1 ricalcola soluzione vicino alla stima
end

Pf.ipolar=1;   %OK
%Pf.ipolar=-1;   
%Pf.ipolar=2;   

iret=1;

Pf.ifund=1;    % 1 fund, 2 first order
%Pf.ifund=2;    % 1 fund, 2 first order
%Pf.ifund=0;    % tutti

if pasnu==1
 Pf.ifund=2;    % devo partire da mm=0
end

Pf.emme='emme_any0';
%Pf.emme='emme_anyB_bot';
global fun_wolu


Ps.iKexact=0;    % =1 ricalcola soluzione vicino alla stima
Ps.Rc_acc=2;  
Ps.Rc_acc=0;  
Ps.iret=1;  %0 per metodo gamma, 1 metodo Tras

Ps.iret_fisso=1;  % se non cambia reticolo, calcolo K una sola volta
Ps.igraef_new=1;    % =1 usa nuovi indici
Ps.igraef_new=0;    % =1 usa nuovi indici
Ps.igraef_new=2;    % =1 usa nuovi indici
Ps.ired_ret=0;
igr_app=1; % =0: nothing, =1: effective refractive index with infinite lateral extension, =2: effective refractive index with lateral extension defined by Delret
igr_app=0; % =0: nothing, =1: effective refractive index with infinite lateral extension, =2: effective refractive index with lateral extension defined by Delret
Ps.igr_app=igr_app;


iULM=0;
iULM=1;
fun_wolu='GEN_berlin';
%fun_wolu='GEN_woluB';

%Pf.ipolar=-1;
if iret==0
%modacc=0;   %modi accoppiati
end

iSALVA=1;
nosa='iskV';



iBEL=1995
iBEL=2012
tv=0;
save tv tv


Ps.igacrit=1;
 
RAD{1}='esp';  % dret-DC
%RAD{1}='MueMapLar2ver';  % dret-DC


%alim=[0 .4];
if inotip==1
 alim=[0 .15];
 Nk=[80]; %OK
 RAD{1}=[RAD{1},'0'];   
 alim=[0 .32];
 Nk=[100]; %OK 
elseif inotip==2
 alim=[0 .15];
 Nk=[50]; %OK
 RAD{1}=[RAD{1},'2'];  
 ifiez=0;
elseif inotip==0 
 alim=[0 .25];
 alim=[0 .28];
 alim=[0 .1];
 %alim=[0 .15];
 if iret==0
  %alim=alim*.8;
 end 
 %alim=[0 .07];
 Nk=[15]; %OKR 
 %Nk=[30]; %OKR 
 %Pf.nKK=[45 60 75 100];
 %Pf.alve=[.15 .2 .25 .3];
 Pf.nKK=[75 100 125 150];
 Pf.alve=[.3 .3 .3 .3];
 Pf.nKK=[75 100 125 150];
 Pf.alve=[.26 .26 .26 .26]; 
 Pf.nKK=Nk;
 Pf.alve=[.3];  
end
%Nk=[50]; %OKR



%' qui', keyboard

Ps.dissfun='diss_nst1';
Ps.dissfun='diss_hcg';
Ps.dissfun='diss_new';
%Ps.igacrit=2;   % solo plot Geq
rdisp=3;
Ps.rdisp=rdisp;





%       1   2   3   4   5   6  7  
iBEv=[iBEL 10 10 10 10 11 11];
iTov=[  0  0   0   0  0 .1 .1];

Pf.nmasce=-2;
if pasnu==1
 Pf.nmasce=6;
end

    Pf.Ntol=5;
    
    Pf.tol=.1;

% isk=[1 1 1 1 0 1 1];
 isk=[0 1 1 1 1 1 1]; 
 isum=[0 0 0 0 0 0 0 ]+10;
 eval(['save ',nosa,' isk Mis ifp ipro iSALVA fiload iBEL isum Pf RAD iBEv iTov nosa fun_wolu'])





N_sim=1;
if isk(N_sim)==0
N_simf=N_sim+isum(N_sim);

i1D=1;
i1D=0;

rad=RAD{N_sim};
iBEL=iBEv(N_sim);
Pf.tol=iTov(N_sim);

Rc=[0];  % da calcolare dentro con punta in base a base maggiore e minore

Rc=[0.4 0.5 1];
Lbuf=[6.5 26 33];
Base=-[3 4 4.5]/2;

Rc=[0.4];
Lbuf=[6.5];
Base=-[3]/2;


Rc=[.3:.2:2.5];
Lbuf=15*ones(size(Rc));
Base=-3*ones(size(Rc));
Par.ri=Base;

dLb=[0:.05:1];
Lbuf=.7+dLb
%Lbuf=2.6
Lbuf=1.2
%Lbuf=1.
Base=0;
Rc=0;
Par.ri=Base;

if Ps.Rc_acc==2
 c1=length(Rc)==length(Lbuf);
 c2=length(Rc)==length(Base);
 if c1*c2~=1
  'errore punta ', keyboard
 end
end

Metal=2.5;
Conf=0;
Pf.iplan=0;
Pf.ianti_gui=0;
Par.arc=[110:20:510]; 
Par.arc=[180:20:220]; 
Par.arc=[60]; 
DenteOxHCG=[20];
Par.arc=[80]; 
DenteOxHCG=[0];

if iULM==1
      Pf.sh=NaN;      
       Pf.fi=0;
       Pf.fu=1;       
       Pf.N=5;      
      if iret==1
%       str='hcg_VDAYflat';         
       str='hcg_Ber0';         
       str='hcg_Ber0LZero';         
       str='hcg_mue0';         
%       str='hcg_mueCrit';        
%       str='hcg_mue';        
       %str='hcg_Ber0LZeroC';         
       %str='hcg_mueCrit1';        
%       str='hcg_BerN';         

 
       %if Pf.ipolar==-1
       %Pf.fi=0;
       %Pf.fu=3;        
       %Pf.N=20;
       %end

      else
       str='bottouCo';         
       str='MIMSco_hcgNO';    
       str='DBR_VDAYflat';              
       %Pf.fi=0;
       %Pf.fu=2;
       %Pf.N=10;
%       Pf.sh=3;            
      end
%      str='realULMtipL';         
    if iForm==1
     str='realULMtipPol';         
    end
     if inotip==1
      str='realULMpol';   
     elseif inotip==2
      str='realULMair';   
     end  
 Metal=3.3;
 Metal=3.7;
 Metal=0;
 %Metal=6.5;
 Conf=1.8;
 %Conf=0;
 dox=[2.8]/2; 
 dox=[3 6 9]; 
 dox=3; 
 Conf=[10 12 14];
 Conf=[0];
 Datt=0.5; 
 Datt=0.; 
 Pf.iplan=0;
else
      Pf.fi=0.;
      Pf.fu=0.8;
      Pf.N=3;
      Datt=0;
     if inotip==0
      str='tip760MetVero';      
     end  
      Mesa.Isha_me=0;   %pareti verticali
      Mesa.Drme=3;
      Mesa.Rame=7.5;
      Mesa.StDmin=.4;
      Mesa.IMesa=1;
      Mesa.n_ext=1.6;
      Ps.Mesa=Mesa;      
end
Par.re=Metal;
Par.Conf=Conf;
save struttura str Nk alim modacc
%' cont', keyboard

rplot=max(dox)*2.5;		% r_max per plot campo
Ps.rplot=rplot;
r=1.5;		% rastremazione punta
%r=[.3:.1:1.5];		% rastremazione punta
if inotip~=0
r=[1.];		% rastremazione punta
end

if Ps.Rc_acc==2
 r=[0];		% rastremazione punta
end
ra=r;


air=DenteOxHCG;

glass=0;

Par.NpaDie=[0];        % paia dielettriche prima

Part.Rp=[.01];     

Part.tip=-(-20+j*40);       %numero step radiali sfera punta: intervalli se positivo, nm se negativo
Part.tip=(100-j*10);       %numero step z sfera punta: intervalli di h se positivo, nm se negativo
                                          % parte immaginaria negativa: num intervalli radiali: FISSA il n. di strati
Part.cp=-10;       %passo discretizzazione cono punta: num intervalli se negativo

Np_u=[1 5];
Np_u=[1:2:9];
Np_u=[1:3:7];
Np_u=[3];
Np_b=31;
Par.NpaDie=Np_u;
Par.NpaBack=Np_b;

  DC=[.63 .65 .67];
  DC=[5.5:.25:7.5]/10;
  DC=[50]/100;
  period=[.7];  
  %period=linspace(.4,.8,15);  
  %period=period(5);
  dret=linspace(275,315,13);    
  %dret=dret(11);    
  dret=285;
  Par.dret=dret;
  Par.gradc=DC;
  Par.grape=period;
  noxide=[1.6 1.8 2];
  noxide=[1.75];
  Par.oxide=noxide;
  lossQW=[2000];
  Par.lossQW=lossQW;

%'vedo', keyboard
hcg_bot

end
% varia airgap fitto

