function P = MatPr(csi, P1, P2, M1, M2, polariz)

%
% --------------------------------------------------------------
% Comando:      P = MatPr(csi, P1, P2, M1, M2);
% --------------------------------------------------------------
% Funzione per il calcolo della matrice di proiezione tra i modi
% di una regione uniforme e quelli di una regione stratificata
% periodica.
% --------------------------------------------------------------
% Variabili in ingresso:
% csi   costanti di propazione nella regione uniforme
% P1, P2 matrici contenenti i modi della regione stratificata
%               periodica ottenuti come sviluppo in serie di polinomi
%               di Legendre
% M1, M2 ordine massimo dei polinomi di Legendre usati per lo
%               sviluppo nei due strati della regione periodica
% Le altre variabili d'ingresso entrano attraverso il file
% datret2g.m.
%
% Variabile d'uscita:
% P             matrice di proiezione
% --------------------------------------------------------------
%
global HE teta err numbet c mu0 divis spp spd sps lp ls ld d d1 d2 n1 n2 n3 n4 eps1 eps2
%eval(file);
%%datret2g;

if ~HE & polariz
        fat1 = 1/eps1;
        fat2 = 1/eps2;
else
        fat1 = 1;
        fat2 = 1;
end

y1 = -conj(csi)*d1/2;
y2 = -conj(csi)*d2/2;

k = i*sqrt(2*pi/d);

N1 = 0:M1;
N2 = 0:M2;
Mdiag1 = ((-i).^N1).';
Mdiag2 = ((-i).^N2).';

%
% CREAZIONE DELLA MATRICE DI PROIEZIONE
%

T1 = zeros(numbet, numbet);
T2 = zeros(numbet, numbet);
P = zeros(numbet,numbet);
IM1 = 1./(N1 + 0.5);
IM2 = 1./(N2 + 0.5);

meny1 = find(y1 < 0);
piuy1 = find(y1 > 0);
Jn1 = zeros(length(y1),M1+1);
Jn1(piuy1,:) = besselj(N1+1/2,y1(piuy1).');
Jn1(meny1,:) = i*ones(size(meny1,2),1)*(-1).^N1.*besselj(N1+1/2,abs(y1(meny1)).');

meny2 = find(y2 < 0);
piuy2 = find(y2 > 0);
Jn2 = zeros(length(y2),M2+1);
Jn2(piuy2,:) = besselj(N2+1/2,y2(piuy2).');
Jn2(meny2,:) = i*ones(size(meny2,2),1)*(-1).^N2.*besselj(N2+1/2,abs(y2(meny2)).');

nzy1 = find(y1 ~= 0);
nzy2 = find(y2 ~= 0);
zery1 = find(y1 == 0);
zery2 = find(y2 == 0);
Nd1 = exp(i*conj(csi(nzy1))*d1/2)./sqrt(y1(nzy1));
Nd2 = exp(i*conj(csi(nzy2))*(d1+d)/2)./sqrt(y2(nzy2));
T1(:,nzy1) = -d1/2*fat1*(k * diag(Nd1) * (Jn1(nzy1,:) * (diag(Mdiag1) * P1))).';
T2(:,nzy2) = -d2/2*fat2*(k * diag(Nd2) * (Jn2(nzy2,:) * (diag(Mdiag2) * P2))).';
if size(zery1) == 1
        T1(:, zery1) = -d1*fat1*i/sqrt(d)*P1(1,:).';
        T2(:, zery2) = -d2*fat2*i/sqrt(d)*P2(1,:).';
end

P = (T1 + T2).';        % Matrice di proiezione complessiva

if ~polariz
        P = -P;
end
