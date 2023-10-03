function [TeGF, ThFG] = MatProiezPerUnifAnalitSkewRight(phi,csi,ky,n1,d1,n2,d2,VTE1, VTE2, VTE1H, VTE2H,...
    ITE1, ITE2, ITE1H, ITE2H ,kzG_TE,ZG_TE,VTM1, VTM2, VTM1H, VTM2H,ITM1, ITM2, ITM1H, ITM2H,kzG_TM,ZG_TM,segem)

if exist('segem')==1
seg=segem;
else
seg=1;
end
% segno relativo TE TM; =1, cambio segno relativo fra TE e TM rispetto convenzione fun. gen. scalare.

% Matrici di proiezione modi TEx, TMx della PSW su modi di Floquet TEz, TMz
% --------------------------------------------------------------
%       R. Orta Dicembre 2010
% --------------------------------------------------------------
% Funzione per il calcolo della matrice di proiezione tra i modi
% del semispazio di destra e quelli di una regione stratificata
% periodica (Modi di Floquet analitici TEz e TMz).
% Calcola il prodotto interno simmetrico sul modo di Floquet exp(-j*csi*x) ossia
% int(exp(-j*csi*x)*Pn(x)dx)
%
% --------------------------------------------------------------
% Variabili in ingresso:
% csi   costanti trasversali kx nella regione uniforme
% n1, d1, n2 , d2 indice e spessore dei due strati
% VTE1, VTE2 ecc matrici contenenti i modi della regione stratificata
%               periodica ottenuti come sviluppo in serie di polinomi
%               di Legendre
%
% Variabile d'uscita:
%
%
%                   |  TeGTEFTE        TeGTEFTM |
%% MATRICE TeGF=    |                           |
%                   |  TeGTMFTE        TeGTMFTM |
%
% indice di riga m indica il modo della regione periodica
% indice di colonna p indica il modo di Floquet
% prima i modi TE poi i TM
%  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%                   |  ThFTEGTE        ThFTEGTM |
%% MATRICE ThFG=    |                           |
%                   |  ThFTMGTE        ThFTMGTM |

% indice di riga p indica il modo di Floquet
% indice di colonna m indica il modo della regione periodica
% prima i modi TE poi i TM
%  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

%
d=d1+d2;
M1=size(VTE1,1)-1; % grado max polinomi in dielettrico 1
M2=size(VTE2,1)-1; % grado max polinomi in dielettrico 2
Nmodi =  size(VTE1,2);
% % Calcolo TeGF
ktF=psqrt(csi.^2+ky^2);

col=cos(phi)*ones(Nmodi,1);% per risolvere la singolarita' per theta=0
ind0=find(ktF~=0);
col(ind0)=(csi(ind0)./ktF(ind0)).';
coeffkxkt=col*ones(1,Nmodi);
% 
col=sin(phi)*ones(Nmodi,1);% per risolvere la singolarita' per theta=0
ind0=find(ktF~=0);
col(ind0)=(ky./ktF(ind0)).';
coeffkykt=col*ones(1,Nmodi);
% 

% Calcolo TeGF
% calcolo trasformate
VTEHtransf = Transform(csi,n1,d1,VTE1H,n2,d2,VTE2H,'plain');
ITMHtransf = Transform(csi,n1,d1,ITM1H,n2,d2,ITM2H,'plain');
ITEHtransf = Transform(csi,n1,d1,ITE1H,n2,d2,ITE2H,'plain');

rig=(ky*ZG_TE./kzG_TE).';
coeffkyZsukz=ones(Nmodi,1)*rig;

TeGTEFTE=(j*coeffkykt.*coeffkyZsukz.*ITEHtransf+j*coeffkxkt.*VTEHtransf).';
TeGTEFTM=(-j*coeffkxkt.*coeffkyZsukz.*ITEHtransf+j*coeffkykt.*VTEHtransf).';
TeGTMFTE=(-j*coeffkykt.*ITMHtransf).';
TeGTMFTM=(j*coeffkxkt.*ITMHtransf).';

TeGF=[seg*TeGTEFTE,TeGTEFTM; seg*TeGTMFTE,TeGTMFTM];

% Calcolo ThFG
% calcolo trasformate
ITMtransf = Transform(-csi,n1,d1,ITM1,n2,d2,ITM2,'plain');
ITEtransf = Transform(-csi,n1,d1,ITE1,n2,d2,ITE2,'plain');
VTEtransf = Transform(-csi,n1,d1,VTE1,n2,d2,VTE2,'plain');

rig=(ky*ZG_TE./kzG_TE).';
coeffkyZsukz=ones(Nmodi,1)*rig;

ThFTEGTE=-j*coeffkxkt.*VTEtransf+j*coeffkykt.*coeffkyZsukz.*ITEtransf;
ThFTEGTM=j*coeffkykt.*ITMtransf;
ThFTMGTE=-j*coeffkykt.*VTEtransf-j*coeffkxkt.*coeffkyZsukz.*ITEtransf;
ThFTMGTM=-j*coeffkxkt.*ITMtransf;


ThFG=[seg*ThFTEGTE, seg*ThFTEGTM;ThFTMGTE,ThFTMGTM];

