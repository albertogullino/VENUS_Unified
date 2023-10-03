%function [Beff, M1, M2, P1, P2] = stimbet(polariz, lambda, file)
function [Beff, M1, M2, P1, P2] = stimbet(polariz, lambda)

%
% -------------------------------------------------------
% Riga di comando:
%       [Beff, M1, M2, P1, P2] = stimbet(polariz, lambda)
% -------------------------------------------------------
% Funzione per il calcolo delle costanti di propagazione
% in una regione stratificata periodica e delle relative
% autofunzioni. Il tutto e' ricavata da un'equazione agli
% autovalori.
% Si calcola inoltre il numero di polinomi di Legendre in
% ogni strato in modo tale da avere un certo errore
% massimo sulle costanti di propagazione tramite una
% formula empirica ricavata da Vito Lancellotti (vedi
% articolo).
% Le autofunzioni in uscita sono normalizzate in potenza.
% -------------------------------------------------------
% Variabili d'ingresso:
% polariz       =0 --> TE
%                       =1 --> TM
% lambda        lunghezza d'onda a cui vengono fatti i calcoli
% Le altre variabili utili vengono lette dal file
% datret2g.m.
%
% Variabili d'uscita:
% Beff          vettore delle costanti di propagazione al
%                       quadrato
% M1, M2        ordine massimo dei polinomi di Legendre per
%                       lo sviluppo nei due strati della regione period.
% P1, P2        matrici degli autovettori relativi alle auto-
%                       funzioni modali nei due strati della regione
%                       stratificata.
% -------------------------------------------------------
%

global HE teta err numbet c mu0 divis spp spd sps lp ls ld d d1 d2 n1 n2 n3 n4 eps1 eps2
%eval(file);
%datret2g;

LAMBDA = lambda;
KBd = 2*pi*d/LAMBDA*n3*sin(teta*pi/180);
% K di Bloch moltiplicata per d

k02 = (2*pi/LAMBDA)^2;

fd1 = 1-(1-eps1)*polariz;
fd2 = 1-(1-eps2)*polariz;

%
% Scelta del numero di polinomi di Legendre da usare nei due strati
% della regione stratificata per calcolare i beta quadro con un errore
% inferiore a err (formula di Lancellotti - vedi articolo relativo)
%
kx = KBd/d + pi/d*(numbet-1);
z = kx * d1/2;
M1 = round(4.79*z^(7/24)*sqrt(((log(err)+0.69)/14.52).^2+1)+0.96*z-1.14);

z = kx * d2/2;
M2 = round(4.79*z^(7/24)*sqrt(((log(err)+0.69)/14.52).^2+1)+0.96*z-1.14);

% M1+1 e M2+1 sono il numero di polinomi di
% Legendre presi rispettivamente nelle regioni
% larghe d1 e d2 presenti nella struttura
% periodica stratificata

alfa = 2/d1;
beta = 2/d2;
%
% CREAZIONE MATRICE pincTMx PER LA DEFINIZIONE DEI p
% DI ORDINE SUPERIORE
%
L = ones(4);
L(1,3) = -(-1)^(M2 -1);
L(1,4) = -(-1)^M2;
L(2,1) = alfa/fd1 * (M1-1)*M1/2;
L(2,2) = alfa/fd1 * (M1+1)*M1/2;
L(2,3) = -beta/fd2 * (-1)^M2*(M2-1)*M2/2;
L(2,4) = -beta/fd2 * (-1)^(M2+1)*(M2+1)*M2/2;
L(3,1) = -exp(-i*KBd)*(-1)^(M1-1);
L(3,2) = -exp(-i*KBd)*(-1)^M1;
L(4,1) = -exp(-i*KBd)*alfa/fd1 * (-1)^M1 * (M1-1)*M1/2;
L(4,2) = -exp(-i*KBd)*alfa/fd1 * (-1)^(M1+1) * (M1+1)*M1/2;
L(4,3) = beta/fd2 * (M2-1)*M2/2;
L(4,4) = beta/fd2 * (M2+1)*M2/2;

A = -ones(4,M1+M2-2);

m1 = 0:M1-2;
m2 = 0:M2-2;
A(1,m2+M1) = (-1).^m2;
A(2,m1+1) = -alfa/fd1*m1.*(m1+1)/2;
A(2,M1+m2) = beta/fd2*(-1).^(m2+1).*m2.*(m2+1)/2;
A(3,m1+1) = exp(-i*KBd)*(-1).^m1;
A(4,m1+1) = exp(-i*KBd) * alfa/fd1 * (-1).^(m1+1).*(m1+1).*m1/2;
A(4,M1+m2) = -beta/fd2 * m2.*(m2+1)/2;

pinc = inv(L) * A;

%
% CALCOLO MATRICE B GENERALE DALLA QUALE SI RICAVANO I BETA
% COME SUOI AUTOVALORI
%

B = zeros(M1+M2-2, M1+M2-2);

for n = 0:M1-2
        m = n:2:M1-2;
        B(n+1,m+1) = alfa^2 * (n+1/2) * (m+n+1).*(m-n);
        if max(m) == M1-2
                B(n+1,:) = B(n+1,:) + alfa^2 * (n+1/2) * (M1+n+1)*(M1-n)*pinc(2,:);
        else
                B(n+1,:) = B(n+1,:) + alfa^2 * (n+1/2) * (M1+n)*(M1-1-n)*pinc(1,:);
        end
end

for n = 0:M2-2
        m = n:2:M2-2;
        B(n+M1,m+M1) = beta^2 * (n+1/2) * (m+n+1).*(m-n);
        if max(m) == M2-2
                B(n+M1,:) = B(n+M1,:) + beta^2 * (n+1/2) * (M2+n+1)*(M2-n)*pinc(4,:);
        else
                B(n+M1,:) = B(n+M1,:) + beta^2 * (n+1/2) * (M2+n)*(M2-1-n)*pinc(3,:);
        end
end

De = diag(k02*[eps1*ones(1,M1-1), eps2*ones(1,M2-1)]);
B = De + B;

[VET, Beff] = eig(B);
Beff0 = diag(Beff);
Beff = diag(Beff);

%' cont Bet 1', keyboard

%idu=find(imag(Beff)<=0);
%Beff=Beff0(idu);
%VET=VET(:,idu);


[BE, IN] = sort(real(Beff));
IN = flipud(IN);
IN = IN(1:numbet);

Beff = Beff(IN);

% CONTROLLO SUI BETA2: DEVONO AVERE TUTTI PARTE
% IMMAGINARIA MINORE DI 0. ALTRIMENTI POSSONO
% PROVOCARE DEI PROBLEMI SULLE AMMETTENZE DELLA
% REGIONE PERIODICA: POTREBBERO RISULTARE CON
% PARTE REALE MINORE DI 0 E QUINDI SI OTTERREB-
% BE UNA STRUTTURA ATTIVA ==> ELEMENTI DELLA
% MATRICE SCATTERING MAGGIORI DI 1

%inbet = find(angle(Beff) > 0 & angle(Beff) < 1e-5);
%Beff(inbet) = real(Beff(inbet));
%inbet = find(angle(Beff) > (pi-1e-5) & angle(Beff) < pi);
%Beff(inbet) = real(Beff(inbet));
inbet = find(angle(Beff) > 1e-5 & angle(Beff) < (pi-1e-5));
if length(inbet) > 0
        disp('WARNING!!! ');
        disp('BETA2 CON PARTE IMMAGINARIA POSITIVA');
        disp('LA STRUTTURA E'' VISTA COME ATTIVA');
        disp(' ');
%        pausak
end

P1 = zeros(M1+1,numbet);
P2 = zeros(M2+1,numbet);

pin = pinc * VET(:,IN);
P1 = VET(1:M1-1, IN);
P1(M1:M1+1,:) = pin(1:2,:);
P2 = VET(M1:M1+M2-2,IN);
P2(M2:M2+1,:) = pin(3:4,:);

Norm = zeros(numbet,1);
N1 = 0:M1;
N2 = 0:M2;
IM1 = 1./(N1+1/2);
IM2 = 1./(N2+1/2);
Norm = 1./sqrt(1/fd1*d1/2*IM1*(abs(P1).^2)+1/fd2*d2/2*IM2*(abs(P2).^2));

P1 = P1 * diag(Norm);
P2 = P2 * diag(Norm);


%'stimbet1', keyboard