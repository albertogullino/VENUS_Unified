mode.flgBTJ=0;  % flag to TJ presence
mode.flgBTJ_lithographic=0; % 0: no etching; 1: etching in VELM only; 2: etching in VELM and in DD

irel=0; % if 1, relief; if 0, standard VCSEL

mode.NumJac=0;  % flag for numerical jacobian computation
mode.Vnum=0;

% Shunt resistance section - by D'Alessandro
mode.VmaxStable=4.5; % turns off the recombination after VmaxStable
mode.Shunt=0;  % to enable Shunt generation rate
mode.K=-1e24;  % coefficient U=K(phi-phi_0)
mode.Ymat=0;   % adds a parallel conductance when mode.IdriveON=1;
mode.MoveON=0; % impose res=0 after a max number of iterations
mode.ImoveON=0.5; % current threshold to disable mode.MoveON
mode.MAXiterMoveOn=8; % max iteration when mode.MoveON is enabled
mode.TQWFake=0; 
mode.Elementi=1;    % DD with temperature on the elements (1) or nodes (0) 
mode.ThermalDB=0; 
mode.ThermalFake=0;  %interpolated Temperature

Effetti=0;
STIMA_Temp=1;
%'passo dot', pausak
% Restart feature
irest=IREST; % restart flag; 0 restarts from 0, 2 restarts from last point
if IPAR~=11
    T0=273+Temperature;         % environment temperature, K
end

% Needed for parametric loops
T0a=T0-273;
if T0a==110
    Vadd=[1.55:0.02:2.2];
    Imassimo=8;  % massima corrente analizzata
elseif T0a==80
    Vadd=[1.55:0.03:2.4];
    Imassimo=12;  % massima corrente analizzata
elseif T0a==50
    Vadd=[1.55:0.04:2.7];
    Imassimo=14;  % massima corrente analizzata
elseif  T0a==20
    Imassimo=16;  % massima corrente analizzata
    Vadd=[1.55:0.05:2.8] ;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VCSEL, and environmental parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T300=293;         % environment temperature, K

% mode.quasi1D=1;
mode.quasi1D=0;
mode.stairs=0;

%  mode.heteromode=2; % if 2 puts the "heteropoints"; if 1 doesn't
mode.heteromode=1; % if 2 puts the "heteropoints"; if 1 doesn't

% Scaling of equations/matrices/unknowns
mode.CarrierNorm=1; % 1D: 1e8; 3D: 1
mode.CarrierNorm2D=mode.CarrierNorm;

mode.C1=1;   % matrix scaling (in place of mode.CarrierNorm)
mode.EqPermutazione=0;

mode.Cpot=1e3; % 1e3 original VENUS

if mode.quasi1D==0
    % pin
    MEsa=[22.5:1:27.5]+.35; % VENUS3 mesa vector (JSTQE)
    % MEsa=[27.5:1:32.5]; % mesa vector for oxide diameters [6:11] um, VENUS1
    OX=MEsa-10.75*2;
    
    % DelOx=2
    DelOx=1.5;
%     DelOx=1 % smaller value --> earlier threshold of superior modes
    
    % Radial grading of the oxide
    OX1=OX+0.4;
%     OX1=OX+0.2
    OX2=OX1+DelOx;
    OX3=OX2+DelOx;
    
%     OX1=OX;
%     OX2=OX1;
%     OX3=OX2;
    
    OxTot=[OX; OX1; OX2; OX3;];
    OxT=OxTot(:,Isize)/2;
    
    Width_Contact=6;
    Relief=0;      % oxide radius, um
    
    COnt=MEsa-Width_Contact*2; % VENUS3 mesa vector
    Contact=COnt(Isize)/2;      % Contact radius, um
    Contact_e=Contact+Width_Contact;    
else
    % quasi 1D simulation
    OX=9.52*ones(1,Isize);    % DIAMETERS (D1ANA - MDPIAS)
    %             OX=10.50*ones(1,Isize);    % DIAMETERS (D1ANA - MDPIAS) - TJ fitting od 3D case
    
    Contact=OX;
    Contact_e=OX;
    
    Relief=0;      % oxide radius, um
end


mode.OptScaling=0;    % area reduction w.r.t. to the eletrical area
if mode.quasi1D==1 && mode.flgBTJ==1
    mode.fPdifScaling=1.15; % 1 for Oxide or 3D simulations; 1.2 for TJ (1D)
else
    mode.fPdifScaling=1; % 1 for Oxide or 3D simulations; 1.2 for TJ (1D)
end

mode.fPdifScaling=1; % 1 for Oxide or 3D simulations; 1.2 for TJ (1D)

List_strNames

strSave=STR{iStruttura}  % for saving purpose

for ks=1:length(STR)
    STR{ks}=[nomeSR,STR{ks}];
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gain and quantum models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% epsNLg=4e-17;
% epsNLg=3e-17;
epsNLg=1.5e-17;

Qc=0.6;
FLos=1;
iTappo=8;       % setting quantum capture model (see assem.m) ; ==8, no tappo in Fermi
Fat_N2=1;
Fat_P2=1;
% tauH=20;         % hole capture time, ps
tauH=5;         % hole capture time, ps
FattoreTauE=2;  % electron/hole capture time ratio
T_tauscat=-50;   % Fat_cap=exp((T-T0)/T_tauscat0).^Tcap_EXP; Inf: exponential
% Tcap_EXP=1;  % Fat_cap=(1+(T-T0)/T_tauscat0).^Tcap_EXP; Inf: exponential
Tcap_EXP=0;  % Fat_cap=(1+(T-T0)/T_tauscat0).^Tcap_EXP; Inf: exponential


Tstart=250;     % Temp inizio coseno Rialzato
tL_high=1e-7;
mode.tL_high=tL_high;
mode.leakage=60; % intervallo T di Raised Cos per leakage
mode.tL_low=tL_high/200;
mode.leakage=0; % intervallo T di Raised Cos per leakage

% nuovi Flag
mode.GapTemp=1;    % variazione del Gap con temperatura
mode.KT=1;		% variazione di KT con temperatura

if mode.quasi1D==0
    mode.BGN=0;     % Enables band gap narrowing (BGN)
else
%     mode.BGN=1;     % Enables band gap narrowing (BGN)
    mode.BGN=0;     % Enables band gap narrowing (BGN)
end

%Tcap_EXP=0;     % Fat_cap=(1+(T-T0)/T_tauscat0).^Tcap_EXP; Inf: exponential
tauRat=1000;      % taucapture/tauescap ratio, only for iTappo=2 (Debernardi)
fat_gain=1;     % factor to be multiplied times LUT parameters (gain, Rsp)
% fat_gain=1e-3;     % factor to be multiplied times LUT parameters (gain, Rsp), EXCLUDE OPTICAL simulation
CN_Auger=.5;        % moltiplica 1e-30 in mw_phmat, original
% CN_Auger=.1        % moltiplica 1e-30 in mw_phmat
% CN_Auger=0        % moltiplica 1e-30 in mw_phmat
FatNP_Auger=1;      % Cppn=FatNP_Auger*Cnnp, vedi mw_phmat
CTemp_Auger=1.;     % beta_aug in (8) of 2019Debernardi_JSTQE, original
% CTemp_Auger=.2;     % beta_aug in (8) of 2019Debernardi_JSTQE
FatAuger23D=1;
% FatAuger23D=5
C_Temp_DD=1;


mode.AugerExpulsion=0; % se 1, attiva il modello di Auger expulsion
mode.AugerExpulsion3D=0;
Fat_regeneration=1;

Auger_broad=10;   %(nm)

iLUT=2;         % LUT index (old LUT)
% iLUT=5;         % LUT index (new LUT)

LUT{1}='LUT4D_Jun_Markus_Lorentzian';
LUT{2}='LUT4D_Jun_Markus_nMark_40';
LUT{3}='LUT4D_Julian_nMark';
LUT{4}='LUT4D_Julian_nMark_xmol=0.06';
LUT{5}='LUT4D_Stephan_nMark_xmolQW=0_xmolBarr=0.3QWthick=7.7';


for ks=1:length(LUT)
    LUT{ks}=[nomeLUT,LUT{ks}];
end

Deltalam=3;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thermal and optical coefficients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C_Temp=1;  % Coeff. Temperatura totale
C_TempGain=1;  % Coeff. Temperatura Gain
dndT=2.3e-4;  % dn/dT  2.3e-4 da dati sperimentali
dndT1D=4e-4;  % for 1D fitting
fat_RAD=0.50;   % self-absorption heating from spont. recombination
fat_RAD=5e-2   % self-absorption heating from spont. recombination

mode.FatQtot=115;    % scaling factor for Q to fit 3D and 1D deltaT
mode.FatV=46;   % kT scaling between substrate and VCSEL regions
mode.FatQcontact=1.00; % 1D simulation: scaling at the contact (0.65 for Oxide; 1.00 for TJ)

% fCondTer=1;   % transverse thermal conducibility
% fCondTerZ=.8;   % Longitudinal thermal conducibility
fCondTer=0.9;   % transverse thermal conducibility
fCondTerZ=0.8*fCondTer;   % Longitudinal thermal conducibility


AlTarocco=1;    % multiplication factor for optical absorption heating
TAROCCO=1;      % correction factor for Pst/I slope

if mode.quasi1D==1
    Exp_Temp0=-1.20;
    fatt_dndT=1;  % dn/dT  2.3e-4 da dati sperimentali
else
    Exp_Temp0=-1.30;	% VENUS
    Exp_Temp0=-1.1
%     fatt_dndT=0.95;  % dn/dT  2.3e-4 da dati sperimentali
    fatt_dndT=1.03  % dn/dT  2.3e-4 da dati sperimentali
end
TARde=1;
mode.ABSe=5;
mode.ABSh=11;
mode.ABSe0=3;   % these are multiplied by Fat_Perd0 in ASSEGNO_mode: f_alpha in 2019Debernardi_JSTQE, (7)
mode.ABSh0=7;

% Fat_Perd0=2.6;  % con Log 1
Fat_Perd0=2.9  % con Log 1

% increase ABS_Texp with Tvelm=DeltaT+T0, in f_CallVELM:
% ABS_Texp=mode.ABS_Texp+mode.PerCoefExT*Tvelm;
PerCoefExT=0
%PerCoefExT=2e-3 

% changes Fat_Perd0 in ASSEGNO_mode: Fat_Perd_mod=Fat_Perd0+Fat_PerCoefTemp*(mode.T0-T300);
% Fat_PerCoefTemp=(.9-Fat_Perd0)/90;
% Fat_PerCoefTemp=0   
Fat_PerCoefTemp=0.002
if IPAR==42
    Fat_PerCoefTemp=0;
end


% ABS_Texp=2.4
ABS_Texp=0      % in VELM: ABS.eleccentro=ABS.eleccentro.*(1+Tvelm/T300).^ABS_Texp;
% ABS_Texp=1.2      % in VELM: ABS.eleccentro=ABS.eleccentro.*(1+Tvelm/T300).^ABS_Texp;

ABS_Apor=0;   %dipendenza lineare !!!!!!!
ABS_Ader=0;

mode.ifFC=1;   % =1 usa Free Carriers; =0 usa Droganti; =2: media dei due

mode.xmolPiatto=0;   % =1 toglie tutta la variazione di x_mol

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electron/hole mobilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ExpH=2.5;       % exponent for thermal dependence of hole mobility
% FattoreExpE=0.8;  % electron/hole ratio for thermal dependence of mobility

% ExpH=1;       % exponent for thermal dependence of hole mobility
ExpH=0.4;       % exponent for thermal dependence of hole mobility
FattoreExpE=1;  % electron/hole ratio for thermal dependence of mobility

if mode.quasi1D==0
    IHILS=0;   % 0 fisso, 1 variabile, original VENUS
else
    IHILS=1;   % 0 fisso, 1 variabile
end
% N_X=1.5e17;      % Hilsum model parameter
N_X=2.5e17      % Hilsum model parameter
NxCoe=.011;
if IPAR==4
    NxCoe=0;
end

Fat_Dop=1.;
CTemp_Ion=0;
% mode.Zmat=8.5; % linear network embedding nonlinear device
Zmat=0; % linear network embedding nonlinear device

% mode.DopingMobilityModel='none'; % 'none' or 'Hilsum'
% mode.DopingMobilityModel='Hilsum'; % 'none' or 'Hilsum'
% mode.DopingMobilityModel='Roland'; % 'none' or 'Hilsum' or 'Roland' or 'RolandIon'
mode.DopingMobilityModel='RolandIon'; % 'none' or 'Hilsum' or 'Roland' or 'RolandIon'
%mode.DopingMobilityModel='RolandFC'; % 'none' or 'Hilsum' or 'Roland' or 'RolandIon'
mode.ireno=0;

mode.FatMob=1;
% cot=[3.5e-2 1.2];   % factor for mobility dependence on T: f(T)=cot(1)*T+cot(2)
load COT

 FAT_Diff_E=0.4;   % factor to be multiplied times QW electron mobility
% FAT_Diff_E=0.5   % factor to be multiplied times QW electron mobility
FAT_Diff_H=1;   % factor to be multiplied times QW hole mobility
% FAT_Diff_H=0.2   % factor to be multiplied times QW hole mobility
mode.idiffusioneQW=3;   % 0: no QW diffusion; 1: QW diffusion; 2: NO; 3: drift-diffusion in QW
% mode.idiffusioneQW=2;   % 0: no QW diffusion; 1: QW diffusion; 2: NO; 3: drift-diffusion in QW
mode.iambipolarQW=0;    % use ambipolar mobility in QW
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recombination parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if mode.flgBTJ==1
    mode.GR={'SRH','rad','Auger','BTBT'}; % generation/recombination, {} for none
else
    mode.GR={'SRH','rad','Auger'}; % generation/recombination, {} for none
end
% mode.GR={}; % generation/recombination, {} for none
% mode.taun = 1e-9; % electron SRH time, s
% mode.taup = 1e-9; % hole SRH time, s
mode.taun = 10e-9; % electron SRH time, s
mode.taup = 10e-9; % hole SRH time, s
mode.taunQW = 100e-9; % QW electron SRH time, s
mode.taupQW = 100e-9; % QW hole SRH time, s
% Brad = 1.8e-10; % Brad 3D
Brad = 1.0e-10; % Brad 3D
%
mode.tauExp = 1; % exponent for thermal dependence of SRH times
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Parametri VELM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ipolar=-1; % con 3 calcola le due polarizzazioni separatamente: necessario per forti anisotropie
ipolar=1; %-1 TE, 1 TM. Con 2 calcola le due polarizzazioni separatamente: necessario per forti anisotropie
% ipolar=1; %-1 TE, 1 TM. Con 2 calcola le due polarizzazioni separatamente: necessario per forti anisotropie
VelmOptions.ipolar=ipolar; % 1: modo coseno; se vuoi il modo "seno", metti ipolar=-1; (1, -1, 2 (both))
% VelmOptions.ivett=1;  % VELM vectorial mode
VelmOptions.ivett=0;    % VELM scalar mode
VelmOptions.iany=0;
VelmOptions.imod_acc=0;  % 0 per LP

fat_ag=1;     % fattore antiguiding (per ridurlo o aumentarlo); per toglierlo, agire su ianti_gui [0 o 1]
VelmOptions.ianti_gui=1; % 1 in VENUS quasi 1D (as in D1ANA); 0 in VENUS 3D
% VelmOptions.ianti_gui=0;
VelmOptions.gain_gui=1;

mode.verbVELM=0;
if mode.quasi1D==1
    mode.verbVELM=-1;   % <0 to see VELM results only the first time; >0: always
end
mode.verbVELM=-2;   % <0 to see VELM results only the first time; >0: always

itutmir=0; % if 1, the "entire" optical structure is studied with thermal (strong discretization)


if mode.flgBTJ==1
    NUMERO_MODI_v=[1 1 1 2 3 4 4]; % numero modi in VENUS
else
    NUMERO_MODI_v=[1 1 1 3 4 4 4];
end
NUM_Azim_v=[1 1 1 2 3 3 3]; % modi azimutali

% NUMERO_MODI_v=[1 1 1 1 1 1 1];
% NUM_Azim_v=[1 1 1 1 1 1 1];

if isfield(mode,'quasi1D') & mode.quasi1D==1
    NUMERO_MODI_v=ones(size(NUMERO_MODI_v));
    NUM_Azim_v=ones(size(NUM_Azim_v));
end

% In questo caso considero NUM_Azim = 3 modi azimutali nu. In PALLE: 1  2 e 4
% Per ogni nu, prendo i primi due modi (nmasce=-2), IN FREQUENZA sul grafico (Fint,alvet)
% Li impilo tutti, ordinandoli poi in FREQUENZA  (si tratta del parametro isoga=0)
% Ne prendo i primi NUMERO_MODI, che sono i modi dati in pasto a VENUS.

NUM_Azim=NUM_Azim_v(Isize); % numero massimo variazioni radiali, modi azimutali nu (if =3 -> palle: 1,2,4)
NUMERO_MODI=NUMERO_MODI_v(Isize);   % numero modi in VENUS (nmodes)
Pf.nmasce=-2;  % modi radiali per ogni modo azimutale, in FREQUENZA sul plot(Fint,alvet)
 
VelmOptions.isoga=0;    % order modes in VELM, based on wavelength (0), gain (1)

if NUMERO_MODI>1
    % 0: ordina in FREQUENZA; 1: ordina in GAIN
    VelmOptions.isoga=0;    % SACRO ed intoccabile per VENUS
end
NP_k=[30 30 25 20 20 20];
VelmOptions.NP_k=NP_k(Isize);
% VelmOptions.NP_k=20;
% VelmOptions.NP_k=60; % 2 etched TJ case
VelmOptions.num_azim=NUM_Azim; % example, 3 is 0, 1, 2
VelmOptions.Dlam=[-.5 4.5 5 0 .4]; % ORIGINAL
if NUMERO_MODI<=2
    VelmOptions.Dlam=[-.5 3.5 5 0 .4]; % ORIGINAL
%     VelmOptions.Dlam=[1 5 5 0 .4]; % ORIGINAL, 2 modes, dox=2 um
%     VelmOptions.Dlam=[-10 5 5 0 .4]; % ORIGINAL, Isise=4, 2 TJ (etched)
end
% KMax=[.21 .18 .15 .12];
% KMax=[.18 .15 .14 .12 .11 .10];
KMax=[.25 .20 .16 .12 .11 .10];
VelmOptions.krel_max=KMax(Isize);               %kmax; 0.1 va bene per aperture normali (3-4 um)
% mode.mintempVELM=200; % degrees, minimum temperature such that VELM is called
mode.mintempVELM=.5; % degrees, minimum temperature such that VELM is called
mode.IsoThermal=0;
mode.DT0=100; % degrees, minimum temperature for estimating Dlam_mod
mode.frsp=1/90; %.5 x sopra e sotto,  .3 per cono ricezione, 100 fattore tras specchio
mode.frsp=1/150; %.5 x sopra e sotto,  .3 per cono ricezione, 100 fattore tras specchio
%

% Altri settings per VELM diversi da default
%VelmOptions.dissfun='diss_fillRET1';
%VelmOptions.igraef_new=2;
%VelmOptions.iBWnew=2;
%VelmOptions.iany=3;
Pf.ipolar=ipolar;
%Pf.emme='emme_navyNEWlastAny';
VelmOptions.Pf=Pf;

% Sub Losses, Temperature
mode.ABSTlosses=1; % 1, "new" (T in FCA); 0, "old", (T in denne)

% gestione DD

mode.maxiter=15; % maximum number of iterations (DD+PB)
mode.Verbose=0; % 0: stop if not convergent; 1: verbose; 2: always move on
ContaBisezMax=10;
% mode.tolconv=1e-6; % expected relative tolerance
mode.tolconv=1e-5; % expected relative tolerance
mode.dlossflg=1; % dielectric losses (MW solver)
mode.report=1; % verbose mode
mode.nflg=1; % include electron continuity equation
mode.pflg=1; % include hole continuity equation
mode.tflg=0; % include trap equations
mode.Tflg=1; % include thermal effects
mode.oflg=1; % include quantum effects and stimulated recombination
mode.RAD_spalmato=1;
if mode.oflg==0
    mode.RAD_spalmato=0;
end
mode.BULK=1; % include quantum effects and stimulated recombination

mode.iTfig=0; % if 1 the thermal simulator plots intermediate results, -1 solo la prima volta
mode.maxScheckRepeat=0; % maximum times of Scheck>1 condition before acting
mode.ScheckMultiplicationFactor=200; % multiplication factor to reset Pst
mode.minthermalvoltage=0.01; % min. voltage to activate thermic simulator
mode.Pmin_Pfit=.1; % prediction of PDiss (PDissPred); set 50 to avoid it; 100 for old StimaTempWU=1.05
mode.minPorVELM=.5e12; % 1/cm2, minimum 2D carrier such that VELM is called

mode.tolconv_neutr=1e-13; % expected relative tolerance
mode.tolinvferdr=1e-13; % tolerance on relative invferdr function value
mode.tolGji=1e-12; % tolerance on cylindrical geometry factor (cm)
mode.BankRose=0; % 1: Bank-Rose is ON; 0: Bank-Rose is OFF
mode.maxiter_t=6; % maximum number of inner iterations (Bank-Rose)
mode.mobility='none'; % mobility model
mode.ntrap=0; % number of trap levels
mode.stats='Fermi'; % "Fermi" or "Boltzmann" are available
mode.ionization='Incomplete'; % "Full" or "Incomplete" are available
mode.symmetry='Cylindrical-Y'; % rotation around y axis
mode.taucarrier=0;      %mette Vallone se 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Contact section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mode.ContactPosition='left'; % put "line" contact at "left" or "right".
mode.ContactPosition='right'; % put "line" contact at "left" or "right".

mode.Idrive=1;	% current driving flag
mode.Ilog=0;

% minimum power or applied bias at which current driving is turned ON
mode.Pmin=0.5;
mode.Vmin=1.7;

mode.Istep=0.2e-3/mode.CarrierNorm; % A, current step for current driving below Imin
mode.IstepLarge=0.5e-3/mode.CarrierNorm; % A, current step for current driving above Imin

if T0a==20
    mode.Imin=1.5e-3/mode.CarrierNorm;  % 1.5 in case of "etched" TJ
elseif T0a==50
    mode.Imin=3.5e-3/mode.CarrierNorm;  % 1.5 in case of "etched" TJ
elseif T0a==80
    mode.Imin=3e-3/mode.CarrierNorm;  % 1.5 in case of "etched" TJ
elseif T0a==110
    mode.Imin=2.5e-3/mode.CarrierNorm;  % 1.5 in case of "etched" TJ
    mode.IstepLarge=mode.IstepLarge/2; % A, current step for current driving above Imin
end

if mode.Idrive==1
    Imassimo=18;
end

if strcmp(mode.ContactPosition,'left')==1
    DD0=-DD0;
    Imassimo=-Imassimo;
    mode.Imin=-mode.Imin;
    mode.Istep=-mode.Istep;
    mode.IstepLarge=-mode.IstepLarge;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multiprecision section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mode.mp=0;
mode.Vmp=5;
if mode.mp==1
    addpath('mp')
    addpath('C:\Users\albig\Documents\Multiprecision Computing Toolbox 2019\')
    % Advanpix MultiPrecision Toolbox flags
    %     mp.Digits(150); % number of significat digits
    mp.Digits(100); % number of significat digits
    %     mp.Digits(70); % number of significat digits
    %     mp.Digits(40); % number of significat digits
    %     mp.Digits(34); % number of significat digits
    %     mp.Digits(16); % number of significat digits
    mp.FollowMatlabNumericFormat(true); % if "true" uses standard MATLAB format
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Custom DD start - by Torrelli
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mode.DDLater = 0;   % flag to enable this feature
% Voltage from which you want to start the DD. It must be one of the
% voltages contained in mode.v0_dd. 4 V for 3 TJ, 3 V for 2 TJ.
mode.WantedV = 3;
% Fraction of voltage drop that falls over the AR. If the pin structure is
% considered this is set to 1.
mode.FracAR = 0.93;
% Fraction of voltage drop that falls over the TJ. If the pin structure is
% considered this is set to 0.
mode.FracTJ = 1-mode.FracAR;
% Multiplicative parameter which adjust the actual width of th AR, since
% the spatial extension of the region where the potential drops is lower
% with respect to the actual AR width.
mode.AdjustWAR = 0.57;
% Multiplicative parameter which adjust the actual width of th TJ, since
% the spatial extension of the region where the potential drops is lower
% with respect to the actual AR width. If the pin structure is considered,
% this parameter is not used
mode.AdjustWTJ = 0.32;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adiabatic section (not working)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mode.Vadiab=1e3; % bias at which the adiabatic procedure begins
mode.AdiabVec=logspace(-2,0,6); % logspace applied to the adiabatic quantity (mode.FatAdiab)

if ~exist('ILOOP')
    ILOOP=0;
end
if ILOOP==0
    save Vadd Vadd
    ASSEGNO_mode
end
if exist('kpvet')
    %if IPAR==0 || IPAR==30 || IPAR==14
    irest=IREST; % restart flag; 0 restarts from 0, 2 restarts from last point
    %end
end

mode.ioldTemp=0;