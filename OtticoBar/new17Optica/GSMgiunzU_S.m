function [S11, S12, S21, S22] = GSMgiunzU_S(P, mi1, mpn, HE)

%
% ----------------------------------------------------------
% Riga di comando
%		[S11, S12, S21, S22] = GSMgiunzU_S(P, mi1, mpn, HE);
% ----------------------------------------------------------
% Funzione per il calcolo della matrice scattering della
% giunzione reg.unif-reg.stratif.periodica.
% ----------------------------------------------------------
% Variabili d'ingresso:
% P		matrice di proiezione tra i modi della regione unif.
%		e i modi della regione stratif.
% mi1, mpn	vettori delle ammettenze o delle impedenze (a
%		seconda della formulazione desiderata) delle due 
%		regioni
% HE	variabile che distingue una formulazione dall'altra
%			0 --> HFIE
%			1 --> EFIE
%
% Variabili d'uscita:
% S11, S12, S21, S22 elementi della matrice scattering
%		della giunzione
% ----------------------------------------------------------
%

Iu = eye(length(mi1));
Ip = eye(length(mpn));
Mp = diag(mpn);
RMp = diag(sqrt(mpn));
Mi1 = diag(mi1);
RMi1 = diag(sqrt(mi1));

INV1 = inv(P'*Mi1*P + Mp);

S21 = 2*RMp*INV1*P'*RMi1;
S12 = 2*RMi1*P*INV1*RMp;
S11 = 2*RMi1*P*INV1*P'*RMi1;
S22 = 2*RMp*INV1*RMp;

if HE
	S11 = Iu - S11;
	S22 = Ip - S22;
else
	S11 = S11 - Iu;
	S22 = S22 - Ip;
end