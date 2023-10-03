function [betaz2, kTG2, M1, M2, P1, P2, P1H, P2H, Q1, Q2, Q1H, Q2H] =...
    modiLeg(n1,n2,d1,d2,n3,kx,ky,polariz,lambda,Nmodi,errore)

% R.Orta Dicembre 2010
% -------------------------------------------------------
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
% Gli indici di rifrazione possone essere complessi. Sono ricavati i modi
% della struttura originale e di quella aggiunta.
% Viene verificata la correttezza della biortonormalita' 
% delle autofunzioni del problema originale e di quello aggiunto.
% -------------------------------------------------------
% Variabili d'ingresso:
% polariz	='TE' --> TE
%			='TM' --> TM
% lambda	lunghezza d'onda a cui vengono fatti i calcoli
%
% Variabili d'uscita:
% kzG_TE2		vettore delle costanti di propagazione al
%			quadrato
% M1, M2	ordine massimo dei polinomi di Legendre per
%			lo sviluppo nei due strati della regione period.
% P1, P2 	matrici degli autovettori relativi alle auto-
%			funzioni modali nei due strati della regione
%			stratificata. Nel caso TE rappresentano V(x), nel TM
%			rappresentano I(x)
% P1H,P2H   idem per il problema aggiunto
% Q1,Q2     matrice autovettori: Nel caso TM rappresentano V(x), nel TE
%			rappresentano I(x)
% Q1H,Q2H   idem per il problema aggiunto
%  ind. riga --> ordine polinomio di Legendre
%  ind. colonna --> ordine del modo, quindi VTE1 e VTE22 hanno lo stesso numero
%  di colonne

% -------------------------------------------------------
%

Clight=2.99792458e8; %nm / ns, quindi frequenze in GHz
mu0=4*pi*1e-7;

d=d1+d2; %periodo
k0=2*pi/lambda;
omega=k0*Clight;
eps0=k0/(Clight*omega*mu0);

%KBd = k0*d*n3*sin(theta)*cos(phi);   	
KBd = kx*d;   	
% K di Bloch moltiplicata per d
Phasefact=exp(-j*KBd);
%ky= k0*n3*sin(theta)*sin(phi);
k02 = (k0)^2;

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
if length(indbet) > 1e10 
	disp('WARNING!!! ');
	disp('BETA2 CON PARTE IMMAGINARIA POSITIVA');
	disp('LA STRUTTURA E'' VISTA COME ATTIVA');
	disp(' ');
end
kTG2=beta2;
betaz2=beta2-ky^2; % skew incidence
betaz=psqrt(betaz2);

P1 = zeros(M1+1,Nmodi);
P2 = zeros(M2+1,Nmodi);

pin = pinc * VET(:,IN);
P1 = VET(1:M1-1, IN);
P1(M1:M1+1,:) = pin(1:2,:);
P2 = VET(M1:M1+M2-2,IN);
P2(M2:M2+1,:) = pin(3:4,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    if length(indbet) > 1e10 
%    if length(indbet) > 0 
        disp('WARNING!!! ');
        disp('BETA2 CON PARTE IMMAGINARIA POSITIVA');
        disp('LA STRUTTURA E'' VISTA COME ATTIVA');
        disp(' ');
    end
    
    P1H = zeros(M1+1,Nmodi);
    P2H = zeros(M2+1,Nmodi);
    
    pin = pinc * VET(:,IN);
    P1H = VET(1:M1-1, IN);
    P1H(M1:M1+1,:) = pin(1:2,:);
    P2H = VET(M1:M1+M2-2,IN);
    P2H(M2:M2+1,:) = pin(3:4,:);
    % P1H e P2H sono la rappresentazione dei modi del problema aggiunto nella base di Legendre
    
end %if imag(n1)&imag(n2)&imag(n3)==0
betaz2H=beta2H-ky^2; %skew incidence

%Normalizzazione
N1vet = 0:M1;
N2vet = 0:M2;
IM1 = 1./(N1vet+1/2);
IM2 = 1./(N2vet+1/2);
Norm = 1./sqrt(1/fd1*d1/2*IM1*(P1H.*P1)+1/fd2*d2/2*IM2*(P2H.*P2));
Normmat1=ones(M1+1,1)*Norm;
Normmat2=ones(M2+1,1)*Norm;
P1 = P1 .* Normmat1; % se pol=TE rappresenta V(x), se pol=TM rappresenta I(x) 
P2 = P2 .* Normmat2; % se pol=TE rappresenta V(x), se pol=TM rappresenta I(x) 
P1H = P1H .* Normmat1; % se pol=TE rappresenta V(x), se pol=TM rappresenta I(x) 
P2H = P2H .* Normmat2; % se pol=TE rappresenta V(x), se pol=TM rappresenta I(x) 
% determinazione di I(x) (se pol=TE) o V(x) (se pol=TM)
% Costruzione matrice di derivazione
matD=zeros(M1+1,M1+1);
for ind=1:2:M1
    matD=matD+diag(ones(1,M1+1-ind),ind);
end
col=(2*(0:M1)+1)';
riga=ones(1,M1+1);
coeff=col*riga;

Q1=(matD.*coeff)*P1*alfa;
Q1H=(matD.*coeff)*P1H*alfa;
if strcmp(polariz,'TE')
    Q1=j/(k0*Clight*mu0)*Q1;  %%rappresenta I(x)
    Q1H=j/(k0*Clight*mu0)*Q1H;  %%rappresenta I(x)
elseif strcmp(polariz,'TM')
    Q1=j/(k0*n1^2/(Clight*mu0))*Q1;  %%rappresenta V(x)
    Q1H=j/(k0*n1^2/(Clight*mu0))*Q1H;  %%rappresenta V(x)
end

matD=zeros(M2+1,M2+1);
for ind=1:2:M2
    matD=matD+diag(ones(1,M2+1-ind),ind);
end
col=(2*(0:M2)+1)';
riga=ones(1,M2+1);
coeff=col*riga;

Q2=(matD.*coeff)*P2*beta;
Q2H=(matD.*coeff)*P2H*beta;

if strcmp(polariz,'TE')
    Q2=j/(k0*Clight*mu0)*Q2;  %%rappresenta I(x)
    Q2H=j/(k0*Clight*mu0)*Q2H;  %%rappresenta I(x)
elseif strcmp(polariz,'TM')
    Q2=j/(k0*n2^2/(Clight*mu0))*Q2;  %%rappresenta V(x)
    Q2H=j/(k0*n2^2/(Clight*mu0))*Q2H;  %%rappresenta V(x)
end


% verifica ortonormalita' funzioni modali 
% N1vet=0:Mdie1;
% N2vet=0:Mdie2;
weight1=(2./(2*N1vet+1))'*ones(1,Nmodi);
weight2=(2./(2*N2vet+1))'*ones(1,Nmodi);

% deltaKron significa cio' che dovrebbe essere idealmente una delta di Kronecker
 deltaKron=1/fd1*d1/2*P1H.'*(weight1.*P1)+...
     1/fd2*d2/2*P2H.'*(weight2.*P2);
 offdiag=deltaKron-diag(diag(deltaKron));
 norm=max(abs(diag(deltaKron)-ones(Nmodi,1)));
 ortoP=max(max(abs(offdiag)));
 errorebeta2=max(abs(beta2H-beta2));
 if norm >= 1.e-13
     norm
 end
 if ortoP >= 1.e-2
     ortoP
     keyboard
 end
% if errorebeta2 >= 1.e-13
 if errorebeta2 >= 1.e-5
     errorebeta2
 end
 % verifica ortonormalita' funzioni modali di corrente (caso TE) 

 % weight1=(2./(2*N1vet+1))'*ones(1,Nmodi);
% weight2=(2./(2*N2vet+1))'*ones(1,Nmodi);

% %  Questo pezzetto serve a calcolare il prodotto interno delle correnti nel
% %  solo caso TE
%  innerprod=1/fd1*d1/2*Q1H.'*(weight1.*Q1)+...
%      1/fd2*d2/2*Q2H.'*(weight2.*Q2);
%  %offdiag=deltaKron-diag(diag(deltaKron));
%  
%  kx12=k0^2*n1^2-kTG2;
%  coeff1=ones(M1+1,1)*kx12.';
%  kx22=k0^2*n2^2-kTG2;
%  coeff2=ones(M2+1,1)*kx22.';
%  term1=1/fd1*d1/2*(coeff1.*P1H).'*(weight1.*P1);
%  term2=1/fd2*d2/2*(coeff2.*P2H).'*(weight2.*P2);
%  innerprodAn=-(term1+term2)/(omega*mu0)^2;
% 
%  
 
% %  Questo pezzetto serve a calcolare il prodotto interno delle tensioni per le correnti nel
% %  solo caso TE
%  innerprod1=1/fd1*d1/2*P1H.'*(weight1.*Q1)+...
%      1/fd2*d2/2*P2H.'*(weight2.*Q2);
%  %offdiag=deltaKron-diag(diag(deltaKron));
%  innerprod2=1/fd1*d1/2*Q1H.'*(weight1.*P1)+...
%      1/fd2*d2/2*Q2H.'*(weight2.*P2);
 
% 

% % Calcolo Potenza trasportata da controllare!!!!!!!
%  crosspower=1/fd1*d1/2*P1'*(weight1.*P1)+...
%      1/fd2*d2/2*P2'*(weight2.*P2);
%  selfpower=diag(crosspower);
% 
%  beta = psqrt(beta2);
% if strcmp(polariz,'TM')
%         ymodale = k0./(beta.'*c*mu0);
% else % polariz='TE'
%         ymodale = beta.'./(k0*c*mu0);
% end
% testprog=real(conj(ymodale).*selfpower.');
% 
