function [TeLR,ThRL]=MatProiezPerRLskew(ky,n1R,n2R,d1,d2,VTE1L,VTE2L,VTE1LH,VTE2LH,ITE1L,ITE2L,ITE1LH,ITE2LH,ITM1L,ITM2L,ITM1LH,ITM2LH,...
    VTE1R,VTE2R,VTE1RH,VTE2RH,VTM1R,VTM2R,VTM1RH,VTM2RH,ITM1R,ITM2R,ITM1RH,ITM2RH,...
    kzL_TE,ZL_TE,kzR_TM,ZR_TM)


%       R. Orta  Dicembre 2010
% --------------------------------------------------------------

% --------------------------------------------------------------
% Funzione per il calcolo della matrice di proiezione tra i modi
% di una regione stratificata periodica Left (L) e quelli di una regione stratificata
% periodica Right (R), con la stessa dimensione, lo stesso numero di modi e lo stesso numero di polinomi.
% --------------------------------------------------------------
% Variabili in ingresso:
% n1R, d1R, n2R , d2R indice e spessore dei due strati della regione R
% VTE1L, VTE2L ecc matrici contenenti i coefficienti degli sviluppi in 
%               serie di pol. di Legendre dei modi  della regione stratificata
%               periodica L 
%             
% VTE1R, VTE2R  ecc. matrici contenenti i coefficienti degli sviluppi in 
%               serie di pol. di Legendre dei modi  della regione stratificata
%               periodica R 
%
% Variabile d'uscita:
%
%                   |  TeLTERTE        TeLTERTM |
%% MATRICE TeLR=    |                           |
%                   |  TeLTMRTE        TeLTMRTM |

% indice di riga n indica il modo della regione L
% indice di colonna m indica il modo della regione R
% prima i modi TE poi i TM
%  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%
%                   |  ThRTELTE        ThRTELTM |
%% MATRICE ThRL=    |                           |
%                   |  ThRTMLTE        ThRTMLTM |
%
% indice di riga m indica il modo della regione R 
% indice di colonna n indica il modo della rregione L
% prima i modi TE poi i TM
%  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

% si usa lo stesso numero di polinomi nella regione 1Left e nella 1Right
% si usa lo stesso numero di polinomi nella regione 2Left e nella 2Right
% si usa lo stesso numero di modi nella regione Left e nella  Right

M1=size(VTE1L,1)-1;
M2=size(VTE2L,1)-1;
Nmodi =  size(VTE1L,2);
N1vet = 0:M1;
N2vet = 0:M2;

weight1=(2./(2*N1vet+1))'*ones(1,Nmodi)*d1/2;
weight2=(2./(2*N2vet+1))'*ones(1,Nmodi)*d2/2;


% Calcolo TeLR
% Calcolo proiezioni V e I

VTELH_VTER=VTE1LH.'*(weight1.*VTE1R)+VTE2LH.'*(weight2.*VTE2R);

VTELH_VTMR=VTE1LH.'*(weight1.*VTM1R)+VTE2LH.'*(weight2.*VTM2R);

ITELH_ITMR_n2=1/n1R^2*ITE1LH.'*(weight1.*ITM1R)+1/n2R^2*ITE2LH.'*(weight2.*ITM2R);

ITMLH_ITMR_n2=1/n1R^2*ITM1LH.'*(weight1.*ITM1R)+1/n2R^2*ITM2LH.'*(weight2.*ITM2R);

TeLTERTE=VTELH_VTER;

col=(ky*ZL_TE./kzL_TE);
coeffkyZsukzL=col*ones(1,Nmodi);
rig=(ky./(kzR_TM.*ZR_TM)).';
coeffkykzZR=ones(Nmodi,1)*rig;
TeLTERTM= -coeffkyZsukzL.*ITELH_ITMR_n2-coeffkykzZR.*VTELH_VTMR;

TeLTMRTE=zeros(Nmodi);

TeLTMRTM=ITMLH_ITMR_n2;

TeLR=[TeLTERTE,TeLTERTM;TeLTMRTE,TeLTMRTM];


% Calcolo ThRL
% Calcolo proiezioni V e I

VTERH_VTEL=VTE1RH.'*(weight1.*VTE1L)+VTE2RH.'*(weight2.*VTE2L);

VTMRH_VTEL=VTM1RH.'*(weight1.*VTE1L)+VTM2RH.'*(weight2.*VTE2L);

ITMRH_ITEL_n2=1/n1R^2*ITM1RH.'*(weight1.*ITE1L)+1/n2R^2*ITM2RH.'*(weight2.*ITE2L);

ITMRH_ITML_n2=1/n1R^2*ITM1RH.'*(weight1.*ITM1L)+1/n2R^2*ITM2RH.'*(weight2.*ITM2L);

ThRTELTE=VTERH_VTEL;

ThRTELTM=zeros(Nmodi);

ThRTMLTE=coeffkykzZR.'.*VTMRH_VTEL+coeffkyZsukzL.'.*ITMRH_ITEL_n2;

ThRTMLTM=ITMRH_ITML_n2;
ThRL=[ThRTELTE,ThRTELTM;ThRTMLTE,ThRTMLTM];

