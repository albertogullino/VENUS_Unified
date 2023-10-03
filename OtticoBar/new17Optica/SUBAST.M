
global ico

if length(ico)==0
 ico=0;
end 
%'ico ', keyboard
%ico=0;   %controllo Pierluigi

% Bastonero Sergio
% Torino, 11 maggio 1999
% ---------------------------------------------------------------
% ***************************************************************
% ---------------------------------------------------------------
%      BETA OTTENUTI CON UNO SVILUPPO SU POLINOMI DI LEGENDRE
% ---------------------------------------------------------------
% ***************************************************************
% ---------------------------------------------------------------
% Programma completo utile a ricavare le efficienze di trasmissio-
% ne e riflessione relative al caso di un'onda piana incidente su
% un reticolo a profilo rettangolare e con regione stratificata
% composta solo da due strati a indice diverso (reticolo rettango-
% lare classico).
% E` possibile considerare anche il caso con perdite (quindi indi-
% ci di rifrazione complessi).
% Il metodo usato e' quello del MODE-MATCHING in cui le funzioni
% in cui e' stato espanso il campo sull'apertura sono i modi della
% regione PERIODICA.
%
% I dati della struttura devono essere inseriti nel file testo
% datret2g.m, che deve quindi essere presente nella directory
% di lavoro. Inoltre si puo' scegliere la polarizzazione (TE o TM)
% dell'onda piana incidente e il tipo di formulazione per ricavare
% la matrice scattering (HFIE o EFIE).
%
% Le altre funzioni necessarie sono:
% msqrt.m               per il calcolo della radice quadra di un numero
%                               in modo tale che la parte immaginaria sia
%                               negativa
% stimbet.m     per il calcolo delle costanti di propagazione
%                               e dei modi nella regione stratificata periodica
% MatPr.m               per il calcolo della matrice di proiezione tra
%                               i modi della regione uniforme e quella strati-
%                               ficata
% GSMgiunzU_S.m per il calcolo della matrice scattering della
%                               giunzione reg.unif.-reg.stratif.
% GSMgiunzS_U.m per il calcolo della matrice scattering della
%                               giunzione reg.stratif.-reg.unif.
% andmodp.m             per calcolare l'andamento dei modi della
%                               regione stratificata periodica
% calcpotf.m    per il calcolo della potenza complessa agli
%                               estremi della matrice scattering completa
% EHplot.m              per plottare i campi elettrico e magnetico in
%                               modulo e fase su una giunzione
%
% Il programma permette di fare i calcoli al variare della
% lunghezza d'onda dell'onda piana incidente e della profondita'
% del reticolo.
% ---------------------------------------------------------------
% ***************************************************************
% ---------------------------------------------------------------
%




pl = 1;

%&&&&&
%&&&&&
%for LAMBDA = ls:lp:ld
for LAMBDA = lambda
%&&&&&
%&&&&&

%disp(['LAMBDA = ', num2str(LAMBDA)])

k0 = 2*pi/LAMBDA;
KBd = k0*n3*d*sin(teta*pi/180);

%     CALCOLO DEI BETA
% ************************
%'quia', keyboard
global x1 x2
divis=1000;
x1 = 0:d/divis:d1;
x2 = d1+d/divis:d/divis:d;
x = [x1 x2];

if ico==0
 [Beff, M1, M2, P1, P2] = stimbet1(polariz, LAMBDA);
% [Beff, M1, M2, P1, P2] = stimbet2(polariz, LAMBDA, ico);
else 
 [Beff, M1, M2, P1, P2] = stimbet2(polariz, LAMBDA, ico);
end

BEFF(:,pl) = Beff;

% P1 e P2 CONTENGONO GLI AUTOVETTORI NORMALIZZATI;
% Beff CONTIENE GLI AUTOVALORI ORDINATI IN SENSO DECRESCENTE
%
% M1+1 e M2+1 sono il numero di polinomi di Legendre utilizzati
% rispettivamente sugli strati [0 d1] e [d1 d] della region
% stratificata periodica
M1v(pl) = M1;
M2v(pl) = M2;
%
% COSTANTI DI PROPAGAZIONE LUNGO z NELLA REGIONE UNIFORME
%
q = (numbet-1)/2;
k = 1:q;
Csi(2*k) = -k;
Csi(2*k+1) = k;
csi = KBd/d + 2*pi*Csi/d;
clear k

% ***********************************
% CALCOLO DELLA MATRICE DI PROIEZIONE
% ***********************************
%P = MatPr(csi, P1, P2, M1, M2, file, polariz);
P = MatPr1(csi, P1, P2, M1, M2, polariz);

% *****************************************************
%             CALCOLO MATRICE SCATTERING
% *****************************************************
%fi=find(abs(imag(Beff))<1e-8 | imag(Beff)>0);
fi=find(abs(imag(Beff))<1e-8);
Beff0=Beff;
Beff0i=imag(Beff);
Beff0i(fi)=0;
Beff=real(Beff)+i*Beff0i;
bet = msqrt(Beff);
K1 = msqrt((k0*n3)^2 - csi.^2);
K2 = msqrt((k0*n4)^2 - csi.^2);
%'Beff', keyboard

if polariz
        ypn = k0./(bet*c*mu0);
        yi1 = k0*n3^2./(K1*c*mu0);
        yi2 = k0*n4^2./(K2*c*mu0);
else
        ypn = bet./(k0*c*mu0);
        yi1= K1./(k0*c*mu0);
        yi2= K2./(k0*c*mu0);
end

Yp = diag(ypn);
zpn = 1./ypn;

Yi1 = diag(yi1);
zi1 = 1./yi1;
Yi2 = diag(yi2);
zi2 = 1./yi2;

if ~HE
        [S111, S112, S121, S122] = GSMgiunzU_S(P, yi1, ypn, HE);
        [S2110, S2120, S2210, S2220] = GSMgiunzS_U(P, yi2, ypn, HE);
else
        [S111, S112, S121, S122] = GSMgiunzU_S(P, zi1, zpn, HE);
        [S2110, S2120, S2210, S2220] = GSMgiunzS_U(P, zi2, zpn, HE);
end

%
% CALCOLO DELLA MATRICE SCATTERING COMPLETA
% AL VARIARE DELLA PROFONDITA' DEL RETICOLO
% SCELTA DALL'UTENTE
%
%'bet'
%pausak
punt = 1;
for sp = sps:spp:spd

        if (ld-ls)<lp
%                disp(['h = ',num2str(sp)])
        end
        %
        % Sposto gli assi della seconda matrice verso sinistra
        %
        ex = exp(-i*bet*sp);
        ex = diag(ex);
        S211 = ex*S2110*ex;
        S212 = ex*S2120;
        S221 = S2210*ex;
        S222 = S2220;
        
%        'ex', keyboard
        %
        % Matrice scattering completa derivante
        % dalla cascata di S1 e S2
        %
        I1 = eye(numbet);
        A = inv(I1 - S122*S211);
        S11 = S111+S112*S211*A*S121;
        S12 = S112*(S211*A*S122 + I1)*S212;
        S21 = S221*A*S121;
        S22 = S222+S221*A*S122*S212;

        a1 = eye(2*q+1,1);
        b1 = S11*a1;
        b2 = S21*a1;
        b3 = S12*a1;
        b4 = S22*a1;
        
        s11=b1(1);
        s21=b2(1);
        s12=b3(1);
        s22=b4(1);
        
        sinv=inv(s12);
	Tr(1,1)=s21-s22*sinv*s11;
	Tr(1,2)=s22*sinv;
	Tr(2,1)=-sinv*s11;
	Tr(2,2)=sinv;
	
  St=[S11 S12; S21 S22];
        at = eye(2*(2*q+1),1);
        
        be1=K1;
        be1(1)=j;
        fita1=find(imag(be1)==0);
        be2=K2;
        be2(1)=j;
        fita2=find(imag(be2)==0);        
        
        
        bdif=sum(abs(b1(fita1)).^2)+sum(abs(b2(fita2)).^2);
        
        b0=sum(abs(b1(1).^2)+abs(b2(1).^2));
        
        s11=b1(1);
        s21=b2(1);
        s12=b3(1);
        s22=b4(1);
        
        Tras(1,1)=s21-s22*s11/s12;
        Tras(1,2)=s22/s12;
        Tras(2,1)=-s11/s12;
        Tras(2,2)=1/s12;
%        SW=1-sum(abs(b1(2:3).^2)+abs(b2(2:3).^2))/sum(abs(b1(1)).^2+abs(b2(1).^2));
        SW=1-bdif;	
%' sQw', keyboard
        BR(pl,punt) = (b1(1));
        BRvero(pl,punt) = (b1(1));
        BRm1(pl,punt) = (b1(2));
        BR1(pl,punt) = (b1(3));
        BT(pl,punt) = (b2(1));
        BT1(pl,punt) = (b2(3));
        BTm1(pl,punt) = (b2(2));
        punt = punt + 1;
        
        
        
        
% if isto==1
% ' dentro BAST ', keyboard
% end
end %ciclo for per sp
        pl = pl+1;

%if LAMBDA+lp <= ld
%    save archmm pl BEFF LAMBDA BR BRm1 BR1 BT BTm1 BT1 M1v M2v file polariz
%    clear all
%    load archmm
%    eval(file)
%end

%&&&&&
%&&&&&
end     %ciclo FOR per LAMBDA
%&&&&&
%&&&&&
%' end beta;', keyboard
if ico==0
 return
end 
% ***************************************************
% Controllo della continuità dei campi alla giunzione
% ***************************************************

bp = A * S121 * a1;
ap = S211*A*S121*a1;
as = ex * bp;
TOTpEs = (I1 + S2110)*as;
TOTpHs = (I1 - S2110)*as;
TOTHs = b2;
TOTEs = b2;

TOTE = a1 + b1;
TOTH = a1 - b1;
TOTpE = S121*a1 + (S122 + I1)*ap;
TOTpH = S121*a1 + (S122 - I1)*ap;


[VTh, VTe] = andmodp(x1, d1, eps1, x2, d2, eps2, [M1 M2], P1, P2, polariz);

if ~polariz
        modoe = -i/sqrt(d)*(exp(-i * x.' * csi));
        modoh = i/sqrt(d)*(exp(-i * x.' * csi));
else
        modoe = i/sqrt(d)*(exp(-i * x.' * csi));
        modoh = modoe;
end
%
% CAMPI NELLE SEZIONI 0+ E sp-
%
EP = TOTpE.' * diag(1./sqrt(ypn)) * VTe;
HP = TOTpH.' * diag(sqrt(ypn)) * VTh;
EPs = TOTpEs.' * diag(1./sqrt(ypn)) * VTe;
HPs = TOTpHs.' * diag(sqrt(ypn)) * VTh;
%
% CAMPI NELLE SEZIONI 0- E sp+
%
ET = modoe * diag(1./sqrt(yi1)) * TOTE;
HT = modoh * diag(sqrt(yi1)) * TOTH;
ETs = modoe * diag(1./sqrt(yi2)) * TOTEs;
HTs = modoh * diag(sqrt(yi2)) * TOTHs;

xn=x-x(1);
%figure, plot(xn,ET,xn+d,ET), pausak
%figure, plot(xn,ETs,xn+d,ETs), pausak


if ico==1
  figure, semilogy(xn,abs(VTh(1:4,:).^2)) , 
  a=axis;
  a(3)=1e-4;
  a(4)=max(max(abs(VTh(1:4,:).^2)));
  axis(a)
  pausak
figure, plot(xn,ET,xn+d,ET), pausak
'controoll old', keyboard
end

Kr1 = msqrt((k0*n1)^2 - bet.^2);
Kr2 = msqrt((k0*n2)^2 - bet.^2);

[du1,i1]=sort(abs(imag(Kr1)));
[du2,i2]=sort(abs(imag(Kr2)));
return
%
% Controllo della conservazione della potenza attiva
% tramite il calcolo della potenza complessa
[Pot1, Pot2] = calcpotf(S11, S21, yi1, yi2);
disp(['Errore sulla conservazione della potenza reale = ', num2str(real(Pot1-Pot2))])
%
% Plot dei campi sulla giunzione
EHplot(HT, HP, ET, EP, x, d, numbet, HE, polariz);
