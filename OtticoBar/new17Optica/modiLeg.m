function [beta2, M1, M2, P1, P2, P1H, P2H] = modiLeg(n1,n2,d1,d2,n3,theta,polariz,lambda,Nmodi,errore)

%
% -------------------------------------------------------
% Riga di comando:
%	[Beff, M1, M2, P1, P2] = stimbet(polariz, lambda)
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
% polariz	='TE' --> TE
%			='TM' --> TM
% lambda	lunghezza d'onda a cui vengono fatti i calcoli
% Le altre variabili utili vengono lette dal file
% datret2g.m.
%
% Variabili d'uscita:
% Beff		vettore delle costanti di propagazione al
%			quadrato
% M1, M2	ordine massimo dei polinomi di Legendre per
%			lo sviluppo nei due strati della regione period.
% P1, P2 	matrici degli autovettori relativi alle auto-
%			funzioni modali nei due strati della regione
%			stratificata.
% P1H,P2H   idem per il problema aggiunto
%  ind. riga --> ordine polinomio di Legendre
%  ind. colonna --> ordine del modo, quindi P1 e P2 hanno lo stesso numero
%  di colonne

% -------------------------------------------------------
%

%datret2g;
d=d1+d2; %periodo
KBd = 2*pi*d/lambda*n3*sin(theta);   	
% K di Bloch moltiplicata per d
Phasefact=exp(-j*KBd);

k02 = (2*pi/lambda)^2;

if strcmp(polariz,'TE')
    fd1=1;
    fd2=1;
elseif strcmp(polariz,'TM')
    fd1=n1^2;
    fd2=n2^2;
end


%
% Scelta del numero di polinomi di Legendre da usare nei due strati
% della regione stratificata per calcolare i beta quadro con un errore
% inferiore a err (formula di Lancellotti - vedi articolo relativo)
%
kx = KBd/d + pi/d*(Nmodi-1);
z = real(kx) * d1/2;
M1 = round(4.79*z^(7/24)*sqrt(((log(errore)+0.69)/14.52).^2+1)+0.96*z-1.14);

z = real(kx) * d2/2;
M2 = round(4.79*z^(7/24)*sqrt(((log(errore)+0.69)/14.52).^2+1)+0.96*z-1.14);

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
L(3,1) = -Phasefact*(-1)^(M1-1);
L(3,2) = -Phasefact*(-1)^M1;
L(4,1) = -Phasefact*alfa/fd1 * (-1)^M1 * (M1-1)*M1/2;
L(4,2) = -Phasefact*alfa/fd1 * (-1)^(M1+1) * (M1+1)*M1/2;
L(4,3) = beta/fd2 * (M2-1)*M2/2;
L(4,4) = beta/fd2 * (M2+1)*M2/2;

A = -ones(4,M1+M2-2);

m1 = 0:M1-2;
m2 = 0:M2-2;
A(1,m2+M1) = (-1).^m2;
A(2,m1+1) = -alfa/fd1*m1.*(m1+1)/2;
A(2,M1+m2) = beta/fd2*(-1).^(m2+1).*m2.*(m2+1)/2;
A(3,m1+1) = Phasefact*(-1).^m1;
A(4,m1+1) = Phasefact * alfa/fd1 * (-1).^(m1+1).*(m1+1).*m1/2;
A(4,M1+m2) = -beta/fd2 * m2.*(m2+1)/2;

pinc = inv(L) * A;

%
% CALCOLO MATRICE B GENERALE DALLA QUALE SI RICAVANO I BETA2
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

De = diag(k02*[n1^2*ones(1,M1-1), n2^2*ones(1,M2-1)]);
B = De + B;

[VET, beta2mat] = eig(B);
beta2 = diag(beta2mat);
clear beta2mat
[BE, IN] = sort(real(beta2));
IN = flipud(IN);
IN = IN(1:Nmodi);
beta2 = beta2(IN);

% CONTROLLO SUI BETA2: DEVONO AVERE TUTTI PARTE
% IMMAGINARIA MINORE DI 0. ALTRIMENTI POSSONO
% PROVOCARE DEI PROBLEMI SULLE AMMETTENZE DELLA
% REGIONE PERIODICA: POTREBBERO RISULTARE CON
% PARTE REALE MINORE DI 0 E QUINDI SI OTTERREB-
% BE UNA STRUTTURA ATTIVA ==> ELEMENTI DELLA
% MATRICE SCATTERING MAGGIORI DI 1
indbet = find(angle(beta2) > 0 & angle(beta2) < 1e-5);
beta2(indbet) = real(beta2(indbet));
indbet = find(angle(beta2) > (pi-1e-5) & angle(beta2) < pi);
beta2(indbet) = real(beta2(indbet));
indbet = find(angle(beta2) > 1e-5 & angle(beta2) < (pi-1e-5));
%if length(indbet) > 0 
%	disp('WARNING!!! ');
%	disp('BETA2 CON PARTE IMMAGINARIA POSITIVA');
%	disp('LA STRUTTURA E'' VISTA COME ATTIVA');
%	disp(' ');
%end

P1 = zeros(M1+1,Nmodi);
P2 = zeros(M2+1,Nmodi);

pin = pinc * VET(:,IN);
P1 = VET(1:M1-1, IN);
P1(M1:M1+1,:) = pin(1:2,:);
P2 = VET(M1:M1+M2-2,IN);
P2(M2:M2+1,:) = pin(3:4,:);


% Calcolo dei modi del problema aggiunto
if (imag(n1)==0)&(imag(n2)==0)&(imag(n3)==0)
    % caso senza perdite
    beta2H=beta2;
    P1H=conj(P1);
    P2H=conj(P2);
else
    
    % ci sono perdite in qualche mezzo
    Phasefact=exp(j*KBd);
    
    % CREAZIONE MATRICE pincTMx PER LA DEFINIZIONE DEI p 
    % DI ORDINE SUPERIORE
    %
    % L = ones(4);
    % L(1,3) = -(-1)^(M2 -1);
    % L(1,4) = -(-1)^M2;
    % L(2,1) = alfa/fd1 * (M1-1)*M1/2;
    % L(2,2) = alfa/fd1 * (M1+1)*M1/2;
    % L(2,3) = -beta/fd2 * (-1)^M2*(M2-1)*M2/2;
    % L(2,4) = -beta/fd2 * (-1)^(M2+1)*(M2+1)*M2/2;
    L(3,1) = -Phasefact*(-1)^(M1-1);
    L(3,2) = -Phasefact*(-1)^M1;
    L(4,1) = -Phasefact*alfa/fd1 * (-1)^M1 * (M1-1)*M1/2;
    L(4,2) = -Phasefact*alfa/fd1 * (-1)^(M1+1) * (M1+1)*M1/2;
    % L(4,3) = beta/fd2 * (M2-1)*M2/2;
    % L(4,4) = beta/fd2 * (M2+1)*M2/2;
    
    % A = -ones(4,M1+M2-2);
    % 
    % m1 = 0:M1-2;
    % m2 = 0:M2-2;
    % A(1,m2+M1) = (-1).^m2;
    % A(2,m1+1) = -alfa/fd1*m1.*(m1+1)/2;
    % A(2,M1+m2) = beta/fd2*(-1).^(m2+1).*m2.*(m2+1)/2;
    A(3,m1+1) = Phasefact*(-1).^m1;
    A(4,m1+1) = Phasefact * alfa/fd1 * (-1).^(m1+1).*(m1+1).*m1/2;
    %A(4,M1+m2) = -beta/fd2 * m2.*(m2+1)/2;
    
    pinc = inv(L) * A;
    
    %
    % CALCOLO MATRICE B GENERALE DALLA QUALE SI RICAVANO I BETA2
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
    
    De = diag(k02*[n1^2*ones(1,M1-1), n2^2*ones(1,M2-1)]);
    B = De + B;
    
    [VET, beta2mat] = eig(B);
    beta2H = diag(beta2mat);
    clear beta2mat
    [BE, IN] = sort(real(beta2H));
    IN = flipud(IN);
    IN = IN(1:Nmodi);
    beta2H = beta2H(IN); % beta2 del problema aggiunto
    
    % CONTROLLO SUI BETA2: DEVONO AVERE TUTTI PARTE
    % IMMAGINARIA MINORE DI 0. ALTRIMENTI POSSONO
    % PROVOCARE DEI PROBLEMI SULLE AMMETTENZE DELLA
    % REGIONE PERIODICA: POTREBBERO RISULTARE CON
    % PARTE REALE MINORE DI 0 E QUINDI SI OTTERREB-
    % BE UNA STRUTTURA ATTIVA ==> ELEMENTI DELLA
    % MATRICE SCATTERING MAGGIORI DI 1
    indbet = find(angle(beta2H) > 0 & angle(beta2H) < 1e-5);
    beta2H(indbet) = real(beta2H(indbet));
    indbet = find(angle(beta2H) > (pi-1e-5) & angle(beta2H) < pi);
    beta2H(indbet) = real(beta2H(indbet));
    indbet = find(angle(beta2H) > 1e-5 & angle(beta2H) < (pi-1e-5));
 %   if length(indbet) > 0 
 %       disp('WARNING!!! ');
%        disp('BETA2 CON PARTE IMMAGINARIA POSITIVA');
%        disp('LA STRUTTURA E'' VISTA COME ATTIVA');
%        disp(' ');
%    end
    
    P1H = zeros(M1+1,Nmodi);
    P2H = zeros(M2+1,Nmodi);
    
    pin = pinc * VET(:,IN);
    P1H = VET(1:M1-1, IN);
    P1H(M1:M1+1,:) = pin(1:2,:);
    P2H = VET(M1:M1+M2-2,IN);
    P2H(M2:M2+1,:) = pin(3:4,:);
    % P1H e P2H sono la rappresentazione dei modi del problema aggiunto nella base di Legendre
    
end %if imag(n1)&imag(n2)&imag(n3)==0

%Normalizzazione
N1vet = 0:M1;
N2vet = 0:M2;
IM1 = 1./(N1vet+1/2);
IM2 = 1./(N2vet+1/2);
Norm = 1./sqrt(1/fd1*d1/2*IM1*(P1H.*P1)+1/fd2*d2/2*IM2*(P2H.*P2));
Normmat1=ones(M1+1,1)*Norm;
Normmat2=ones(M2+1,1)*Norm;
P1 = P1 .* Normmat1;
P2 = P2 .* Normmat2;
P1H = P1H .* Normmat1;
P2H = P2H .* Normmat2;
% verifica ortonormalita' funzioni modali 
% N1vet=0:Mdie1;
% N2vet=0:Mdie2;
weight1=(2./(2*N1vet+1))'*ones(1,Nmodi);
weight2=(2./(2*N2vet+1))'*ones(1,Nmodi);

%%% togliere i commenti per analizzare le proprieta' delle funzioni modali

% deltaKron significa cio' che dovrebbe essere idealmente una delta di Kronecker
%  deltaKron=1/fd1*d1/2*P1H.'*(weight1.*P1)+...
%      1/fd2*d2/2*P2H.'*(weight2.*P2);
%  offdiag=deltaKron-diag(diag(deltaKron));
%  norm=max(abs(diag(deltaKron)-ones(Nmodi,1)))
%  orto=max(max(abs(offdiag)))
%  errorebeta2=max(abs(beta2H-beta2))

