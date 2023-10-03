function P = MatProiezPerUnifAnalit(csi,n1,d1,P1,n2,d2,P2, polariz,formulaz)

%
% --------------------------------------------------------------
% Comando:      P = MatProiezPerUnifAnalit(csi,n1,d1,P1,n2,d2,P2,
% polariz,formulaz);
% --------------------------------------------------------------
% Funzione per il calcolo della matrice di proiezione tra i modi
% di una regione uniforme e quelli di una regione stratificata
% periodica (Modi di Floquet analitici).
% Calcola il prodotto interno simmetrico sul modo di Floquet exp(-j*csi*x) ossia
% int(exp(-j*csi*x)*Pn(x)dx)
%
% --------------------------------------------------------------
% Variabili in ingresso:
% csi   costanti di propazione nella regione uniforme
% n1, d1, n2 , d2 indice e spessore dei due strati
% P1, P2 matrici contenenti i modi della regione stratificata
%               periodica ottenuti come sviluppo in serie di polinomi
%               di Legendre
% polariz='TE" o 'TM'
%  formulaz='HFIE' o 'EFIE'
%
% Variabile d'uscita:
% P(m,n)             matrice di proiezione  =<eFloq,efun> se TM HFIE
%                                           =<hFloq,hfun> se TM EFIE
%                                           =<eFloq,efun> se TE HFIE
%                                           =<hFloq,hfun> se TE EFIE
%                                           (uguale alla precedente)
% indice di riga m indica il modo di Floquet
% indice di colonna n indica il modo della regione periodica
% --------------------------------------------------------------
%
d=d1+d2;
M1=size(P1,1)-1; % grado max polinomi in dielettrico 1
M2=size(P2,1)-1; % grado max polinomi in dielettrico 2
Nmodi =  size(P1,2);

if strcmp(polariz,'TM')==1&strcmp(formulaz,'HFIE')==1
%    In questo caso la matrice di proiezione e' <eFloq,efun>
        fat1 = 1/n1^2;
        fat2 = 1/n2^2;
else
% qui si passa se polariz=TM e formulaz=EFIE, per cui calcola <hFloq,hfun>
% oppure se polariz= TE, e ogi formulaz. Infatti <eFloq,efun>=<hFloq,hfun>
        fat1 = 1;
        fat2 = 1;
end

alpha1 = -csi*d1/2;
alpha2 = -csi*d2/2;

%k = i*sqrt(2*pi/d);

N1vet = 0:M1;
N2vet = 0:M2;
Jn1 = SphericalBesselJ(N1vet,alpha1.');
Jn2 = SphericalBesselJ(N2vet,alpha2.');

Poly1 = (2 ./(-i).^N1vet);
Poly2 = (2 ./(-i).^N2vet);
Floq1 = exp(-i*csi*d1/2);
Floq2 = exp(-i*csi*(d1+d)/2);
coeffPoly1=ones(Nmodi,1)*Poly1;
coeffPoly2=ones(Nmodi,1)*Poly2;
coeffFloq1=Floq1.'*ones(1,Nmodi);
coeffFloq2=Floq2.'*ones(1,Nmodi);
T1=d1/2*fat1/sqrt(d)*coeffFloq1.*((coeffPoly1.*Jn1)*P1);
T2=d2/2*fat2/sqrt(d)*coeffFloq2.*((coeffPoly2.*Jn2)*P2);
P = T1 + T2;        % Matrice di proiezione complessiva

