% lista parametri

Svar{1}='Relief';
Svar{2}='Deltalam';
Svar{3}='ExpH';
Svar{4}='N_X';
Svar{5}='tauH';
Svar{6}='T_tauscat';
Svar{7}='Tcap_EXP';
Svar{8}='tauRat';
Svar{9}='fat_gain';
Svar{10}='fCondTer';
%Svar{22}='Fat_CondZ';
Svar{22}='fCondTerZ';

Svar{11}='T0';
Svar{12}='iLUT';
Svar{13}='iStruttura';
Svar{36}='iStruttura';
Svar{37}='iStruttura';
Svar{14}='iTappo';
Svar{15}='FattoreTauE';
Svar{16}='fatt_dndT';
Svar{17}='fat_RAD';
Svar{18}='AlTarocco';
Svar{19}='TAROCCO';
Svar{20}='RatHE';
Svar{21}='Oxide';

Svar{23}='FAT_Diff_E';
Svar{24}='FAT_Diff_H';
Svar{25}='CN_Auger';
Svar{26}='FatNP_Auger';
Svar{27}='CTemp_Auger';
Svar{28}='CTemp_Ion';
Svar{29}='Zmat';

Svar{30}='Effetti';

Svar{31}='Fat_regeneration';
Svar{32}='Auger_broad';
Svar{33}='C_Temp';
Svar{34}='Brad';
Svar{35}='C_TempGain';
Svar{38}='FatAuger23D';
Svar{39}='C_Temp_DD';
Svar{40}='nLC';
Svar{41}='Isize';
Svar{42}='Fat_Perd0';
Svar{43}='Fat_Dop';
Svar{44}='ABS_Texp';
Svar{45}='Width_Contact';
Svar{46}='Isize';
Svar{47}='DelOx';
Svar{48}='Qc';
Svar{49}='Fat_N2';
Svar{50}='Fat_P2';
Svar{51}='Fat_PerCoefTemp';
Svar{52}='FLos';
Svar{53}='ABS_Apor';
Svar{54}='ABS_Ader';
Svar{55}='Exp_Temp0';
Svar{56}='TARde';
Svar{57}='NxCoe';
Svar{58}='fat_ag';
Svar{59}='epsNLg';
Svar{60}='ifFC';
Svar{61}='PerCoefExT';
Svar{62}='STIMA_Temp';

% Gullino
Svar{63}='Fat_VTJ';
% Svar{100} = 'Non_loop';



IPAR=30;
PMAT{IPAR}=[0:7];
PMAT{IPAR}=[0 8:10];
PMAT{IPAR}=[8:10];
NP(IPAR)=length(PMAT{IPAR});

% 1 Relief
IPAR=1;
PMAT{IPAR}=[1.5 2 3 5];
PMAT{IPAR}=[4 2];
PMAT{IPAR}=[2.25];
% PMAT{IPAR}=[3];
% PMAT{IPAR}=[1.5];
NP(IPAR)=length(PMAT{IPAR});

% 2 Delta Lam
IPAR=2;
PMAT{IPAR}=[-5 5];
%PMAT{IPAR}=10;
%PMAT{IPAR}=[2 5 8];
NP(IPAR)=length(PMAT{IPAR});

% 3 Exp Mobilità
IPAR=3;
% FattoreExpE=0.8; % l'esponenziale degli elettroni ï¿½ quello delle lacune per FattoreExpE
PMAT{IPAR}=[0 .75 1.5];
PMAT{IPAR}=[0 0.3 1];
NP(IPAR)=length(PMAT{IPAR});

% 4 Hilsum 
IPAR=4;
PMAT{IPAR}=[5 10 20]*1e18;
PMAT{IPAR}=[1 10]*1e17;
NP(IPAR)=length(PMAT{IPAR});

% 5 taus  cattura
IPAR=5;
%PMAT{IPAR}=[1 5 10];
PMAT{IPAR}=[2 25];
PMAT{IPAR}=[1 10];

NP(IPAR)=length(PMAT{IPAR});

% 6 Temperatura modello tau
IPAR=6;
PMAT{IPAR}=[-10000 -50 -30];
PMAT{IPAR}=[-200 -50 -20];
NP(IPAR)=length(PMAT{IPAR});

% 7 Exp modello tau  (1+DeltaT/taus)^Exp
IPAR=7;
PMAT{IPAR}=[1 1.5 2];
PMAT{IPAR}=[1];
NP(IPAR)=length(PMAT{IPAR});
% 8 Rapporto t_esc / t_cap

IPAR=8;
PMAT{IPAR}=[20 25 30];
PMAT{IPAR}=[50 35 20];
NP(IPAR)=length(PMAT{IPAR});

% 9 fat_gain
IPAR=9;
PMAT{IPAR}=[.3 .6 1];
NP(IPAR)=length(PMAT{IPAR});

% 10 fCondTer
IPAR=10;
PMAT{IPAR}=[1 .6];
PMAT{IPAR}=[.9 .8];
NP(IPAR)=length(PMAT{IPAR});

% 22 fCondTerZ
IPAR=22;
PMAT{IPAR}=[1  .6];
NP(IPAR)=length(PMAT{IPAR});


% 11 Temperatura Peltier
IPAR=11;
% PMAT{IPAR}=273+[20 50 70];
% PMAT{IPAR}=273+[80 50];
PMAT{IPAR}=273+[110 80 50 20];
%PMAT{IPAR}=273+[20 50 80 110];
%PMAT{IPAR}=273+[20 110];
PMAT{IPAR}=273+20;
% PMAT{IPAR}=273+[110 50 20];
if exist('Tpelt')
    PMAT{IPAR}=273+Tpelt;
end 

NP(IPAR)=length(PMAT{IPAR});

% 12 LUT (vedi settings vari)
IPAR=12;
PMAT{IPAR}=[2 1];
NP(IPAR)=length(PMAT{IPAR});
% 13 Struttura (vedi settings vari)
IPAR=13;
PMAT{IPAR}=[10 9 7];
%PMAT{IPAR}=[10];
%PMAT{IPAR}=[9];
%PMAT{IPAR}=[2 6 7];
%PMAT{IPAR}=[2];
NP(IPAR)=length(PMAT{IPAR});

% 14 iTappo  0 (no Tappo), 1 (Tappo), 2 (cap-esc), 3 (0-esc), 4 (1-esc), 6 Mariangela
IPAR=14;
PMAT{IPAR}=[0 1 6 8];
%PMAT{IPAR}=[0];
NP(IPAR)=length(PMAT{IPAR});

% 15 FattoreTauE
IPAR=15;
PMAT{IPAR}=[.5 1 5];
NP(IPAR)=length(PMAT{IPAR});

% 16 fatt_dndT0=1.4; % prova
IPAR=16;
PMAT{IPAR}=[1 1.4];
NP(IPAR)=length(PMAT{IPAR});

% 17 fat_RAD0=0.50;
IPAR=17;
PMAT{IPAR}=[0.5 .7 .9];
NP(IPAR)=length(PMAT{IPAR});

% 18 AlTarocco0=.7;
IPAR=18;
PMAT{IPAR}=[.7 1];
NP(IPAR)=length(PMAT{IPAR});

% 19 TAROCCO0=1;
IPAR=19;
PMAT{IPAR}=[1 0.9];
NP(IPAR)=length(PMAT{IPAR});

% 20 RatHE0=2;
IPAR=20;
PMAT{IPAR}=[1 2];
NP(IPAR)=length(PMAT{IPAR});

% 21 Oxide0=3.25;
IPAR=21;
PMAT{IPAR}=[2.5 3.5 4.5];
%PMAT{IPAR}=[2.5];
NP(IPAR)=length(PMAT{IPAR});


% 23 FAT_Diff_E0=1;
IPAR=23;
PMAT{IPAR}=[.1 .2 .5 1];
PMAT{IPAR}=[.1 .15 .2];
PMAT{IPAR}=[.05 .2 1];
PMAT{IPAR}=[.1 1];
NP(IPAR)=length(PMAT{IPAR});

% 24 FAT_Diff_H0=1;
IPAR=24;
PMAT{IPAR}=[1 2 3];
NP(IPAR)=length(PMAT{IPAR});

% 25 CN_Auger=1;
IPAR=25;
PMAT{IPAR}=[2 3 4];
PMAT{IPAR}=[3.5 5 2.5];
PMAT{IPAR}=[1 10 100];
PMAT{IPAR}=[.3 .5 .8];
PMAT{IPAR}=[.2 1 5];
PMAT{IPAR}=[.25 .5  1];
PMAT{IPAR}=[.2 2];
%PMAT{IPAR}=[.2 .4 .8];
%PMAT{IPAR}=[.1 1 10];
NP(IPAR)=length(PMAT{IPAR});

% 26 FatNP_Auger : moltiplicato per CN dà CP
IPAR=26;
PMAT{IPAR}=[2 3 4];
NP(IPAR)=length(PMAT{IPAR});

% 27 CTemp_Auger=1;
IPAR=27;

PMAT{IPAR}=[0.5  2];
NP(IPAR)=length(PMAT{IPAR});

% 28 CTemp_Ion;
IPAR=28;
PMAT{IPAR}=[1];
NP(IPAR)=length(PMAT{IPAR});

% 29 Z_substrato;
IPAR=29;
PMAT{IPAR}=[-3 0 20 30];      %negativo prende abs(iStruttura)
PMAT{IPAR}=[0 20 30];      %negativo prende abs(iStruttura)
NP(IPAR)=length(PMAT{IPAR});


% 31 Fat_regeneration
IPAR=31;
PMAT{IPAR}=[-1 0 1];
NP(IPAR)=length(PMAT{IPAR});

% 32 Auger_broad (nm)
IPAR=32;
PMAT{IPAR}=[3 5 10];
NP(IPAR)=length(PMAT{IPAR});

% 33 C_Temp
IPAR=33;
PMAT{IPAR}=[.8 .5 0];
PMAT{IPAR}=[1 0.6 0.1];
%PMAT{IPAR}=[0 0.5 1];
%PMAT{IPAR}=[.05];
NP(IPAR)=length(PMAT{IPAR});

% 34 Brad
IPAR=34;
PMAT{IPAR}=[.5 1 1.5]*1e-10;
NP(IPAR)=length(PMAT{IPAR});

% 35 C_TempGain
IPAR=35;
PMAT{IPAR}=[.8 .4 0];
NP(IPAR)=length(PMAT{IPAR});

% 36 Struttura 1 (vedi settings vari)
IPAR=36;
PMAT{IPAR}=[8 9];
NP(IPAR)=length(PMAT{IPAR});

% 37 Struttura 2 (vedi settings vari)
IPAR=37;
PMAT{IPAR}=[8 9];
NP(IPAR)=length(PMAT{IPAR});

% 38 FatAuger23D
IPAR=38;
PMAT{IPAR}=[.1 .3 .6 ];
NP(IPAR)=length(PMAT{IPAR});

% 39 C_Temp_DD
IPAR=39;
PMAT{IPAR}=[.8 .4 0];
NP(IPAR)=length(PMAT{IPAR});

% 40 nLC
IPAR=40;
% PMAT{IPAR}=[1.60 1.65 1.70];
PMAT{IPAR}=[1.60 1.625 1.65 1.675 1.70 1.725];
NP(IPAR)=length(PMAT{IPAR});

% 41 Oxide0 size;O
IPAR=41;
% PMAT{IPAR}=[6 5 4];
PMAT{IPAR}=[1 2 3 4 5 6];
PMAT{IPAR}=[4 5];
% PMAT{IPAR}=[4 6];
%PMAT{IPAR}=[1 3  5];
%PMAT{IPAR}=[1];
NP(IPAR)=length(PMAT{IPAR});




% 42 perdite  Fat_perd
IPAR=42;
PMAT{IPAR}=[0:3];
%PMAT{IPAR}=[1 4];
PMAT{IPAR}=[1.2 1.3 1.4];
PMAT{IPAR}=[1.5 1.75];
PMAT{IPAR}=[2 1.8 1.6];
PMAT{IPAR}=[1 2.5];
%PMAT{IPAR}=[10 7 4];
NP(IPAR)=length(PMAT{IPAR});


% 43 Doping
IPAR=43;
PMAT{IPAR}=[1 0.2 .5 2];
NP(IPAR)=length(PMAT{IPAR});


%44 mode.ABS_Texp=2;
IPAR=44;
PMAT{IPAR}=[0 2 ];
PMAT{IPAR}=[0.5 ];
%PMAT{IPAR}=[4 5];
%PMAT{IPAR}=[1 2 3];
PMAT{IPAR}=[0 1.5 3];
%PMAT{IPAR}=[1.8 1.7 1.6];
%PMAT{IPAR}=[1.7 1.5 1.3];
NP(IPAR)=length(PMAT{IPAR});

% 45 Width_Contact
IPAR=45;
PMAT{IPAR}=[4 6];
%PMAT{IPAR}=[6];
NP(IPAR)=length(PMAT{IPAR});

% 46 Oxide0 size;O
IPAR=46;
PMAT{IPAR}=[1:6];
%PMAT{IPAR}=[1 3  5];
%PMAT{IPAR}=[1];
NP(IPAR)=length(PMAT{IPAR});

% 47 Oxide0 shape
IPAR=47;
PMAT{IPAR}=[.5 1 2];
PMAT{IPAR}=[1 1.5 2];
NP(IPAR)=length(PMAT{IPAR});

% 48 Qc energy bands
IPAR=48;
PMAT{IPAR}=[0.2 0.6];
NP(IPAR)=length(PMAT{IPAR});

% 49 Fat_N2
IPAR=49;
PMAT{IPAR}=[.5  1 10];
NP(IPAR)=length(PMAT{IPAR});

% 50 Fat_P2
IPAR=50;
PMAT{IPAR}=[1 .8 ];
NP(IPAR)=length(PMAT{IPAR});

% 51 Coef FatPerdT
IPAR=51;
PMAT{IPAR}=-[0 1]/100;
NP(IPAR)=length(PMAT{IPAR});

% 52 Ottico
IPAR=52;
PMAT{IPAR}=[1 1000];
NP(IPAR)=length(PMAT{IPAR});

% 53 ABS_Apor
IPAR=53;
PMAT{IPAR}=[0 .1 1 2];
PMAT{IPAR}=[0 .1 .5 1];
%PMAT{IPAR}=[0];
NP(IPAR)=length(PMAT{IPAR});

% 54 ABS_Ader
IPAR=54;
PMAT{IPAR}=[0  1];
NP(IPAR)=length(PMAT{IPAR});

% 55 Exp_Temp
IPAR=55;
PMAT{IPAR}=[-1 -1.5  -2];
PMAT{IPAR}=[-1.5 -1.3  -1.1];
PMAT{IPAR}=[-1.7 -1.3 -0.9];
PMAT{IPAR}=[ -1.7:.2:-.5];
PMAT{IPAR}=[ -2 -0.6];
NP(IPAR)=length(PMAT{IPAR});

% 56 TAR ag
IPAR=56;
PMAT{IPAR}=[-1 -1.5  -2];
PMAT{IPAR}=[-1.5 -1.3  -1.1];
PMAT{IPAR}=[1 10 50];
NP(IPAR)=length(PMAT{IPAR});

% 57 Coef N_X
IPAR=57;

PMAT{IPAR}=[9 6 3 0]/1000;
NP(IPAR)=length(PMAT{IPAR});

% 58 fat_ag
IPAR=58;
PMAT{IPAR}=[0.01 1 ];
NP(IPAR)=length(PMAT{IPAR});

% 59 epsNLg
IPAR=59;
PMAT{IPAR}=[0 3 6]*1e-17;
NP(IPAR)=length(PMAT{IPAR});

% 60 ifFC
IPAR=60;
PMAT{IPAR}=[2 1];
NP(IPAR)=length(PMAT{IPAR});

% 61 CoefExp
IPAR=61;
PMAT{IPAR}=-[2 1 0]/1000;
NP(IPAR)=length(PMAT{IPAR});

% 62 Stima_Temp
IPAR=62;
PMAT{IPAR}=[1 0];
NP(IPAR)=length(PMAT{IPAR});

FattoreExpE=1; % l'esponenziale degli elettroni = quello delle lacune per FattoreExpE

% 63 Fat_VTJ
IPAR=63;
PMAT{IPAR}=linspace(0,1,9);
NP(IPAR)=length(PMAT{IPAR});



for kt=1:length(Svar)
 if kt==4
  Fpar=1e-18;
 else
  Fpar=1;
 end
 if kt==59
  Fpar=1e17;
 end 
 param=Fpar*PMAT{kt};
 if kt==11
  param=param-273;
 end
 eval(['tit{',num2str(kt),'}=[''',Svar{kt},'  ',num2str(param),'''];'])
end 


