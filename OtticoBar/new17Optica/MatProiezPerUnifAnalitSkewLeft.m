function [TeFG,ThGF] = MatProiezPerUnifAnalitSkewLeft1(phi,csi,ky,n1,d1,n2,d2,VTE1, VTE2, VTE1H, VTE2H,...
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
% del semispazio di sinistra e quelli di una regione stratificata
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
%                   |  TeFTEGTE        TeFTEGTM |
%% MATRICE TeFG=    |                           |
%                   |  TeFTMGTE        TeFTMGTM |

% indice di riga m indica il modo di Floquet
% indice di colonna n indica il modo della regione periodica
% prima i modi TE poi i TM
%  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%
%                   |  ThGTEFTE        ThGTEFTM |
%% MATRICE ThGF=    |                           |
%                   |  ThGTMFTE        ThGTMFTM |
%
% indice di riga m indica il modo della regione periodica
% indice di colonna n indica il modo di Floquet
% prima i modi TE poi i TM
%  &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

% --------------------------------------------------------------
%
d=d1+d2;
M1=size(VTE1,1)-1; % grado max polinomi in dielettrico 1
M2=size(VTE2,1)-1; % grado max polinomi in dielettrico 2
Nmodi =  size(VTE1,2);
% Calcolo TeFG
% Calcolo trasformate
VTEtransf = Transform(-csi,n1,d1,VTE1,n2,d2,VTE2,'plain');
VTMtransf = Transform(-csi,n1,d1,VTM1,n2,d2,VTM2,'plain');
ITM_n2transf = Transform(-csi,n1,d1,ITM1,n2,d2,ITM2,'weighted');
ktF=psqrt(csi.^2+ky^2);

col=cos(phi)*ones(Nmodi,1); % per risolvere la singolarita' per theta=0
ind0=find(ktF~=0);
col(ind0)=(csi(ind0)./ktF(ind0)).';
coeffkxkt=col*ones(1,Nmodi);
TeFTEGTE=-j*coeffkxkt.*VTEtransf;

col=sin(phi)*ones(Nmodi,1); % per risolvere la singolarita' per theta=0
ind0=find(ktF~=0);
col(ind0)=(ky./ktF(ind0)).';
%col=(ky./ktF).';
coeffkykt=col*ones(1,Nmodi);
rig=(ky./(kzG_TM.*ZG_TM)).';
coeffkykzZ=ones(Nmodi,1)*rig;
TeFTEGTM=j*coeffkykt.*ITM_n2transf+j*coeffkxkt.*coeffkykzZ.*VTMtransf;

TeFTMGTE=-j*coeffkykt.*VTEtransf;

TeFTMGTM=-j*coeffkxkt.*ITM_n2transf+j*coeffkykt.*coeffkykzZ.*VTMtransf;

TeFG=[seg*TeFTEGTE, seg*TeFTEGTM; TeFTMGTE,TeFTMGTM];

% Calcolo ThGF
% Calcolo trasformate
VTEHtransf = Transform(csi,n1,d1,VTE1H,n2,d2,VTE2H,'plain');
VTMHtransf = Transform(csi,n1,d1,VTM1H,n2,d2,VTM2H,'plain');
ITMH_n2transf = Transform(csi,n1,d1,ITM1H,n2,d2,ITM2H,'weighted');

ThGTEFTE=(j*coeffkxkt.*VTEHtransf).';
ThGTEFTM=(j*coeffkykt.*VTEHtransf).';
ThGTMFTE=(j*coeffkxkt.*coeffkykzZ.*VTMHtransf-j*coeffkykt.*ITMH_n2transf).';
ThGTMFTM=(j*coeffkykt.*coeffkykzZ.*VTMHtransf+j*coeffkxkt.*ITMH_n2transf).';

ThGF=[seg*ThGTEFTE, ThGTEFTM; seg*ThGTMFTE, ThGTMFTM];
