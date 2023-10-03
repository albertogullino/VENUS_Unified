function [S11, S12, S21, S22] = GSMgiunzU_Sorta(Teab,Thba, mFloq3, mmodale, formulaz)

%RO Apr 2010 apparentemente usa la versione ab del mode matching, pag 25
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
Nmodi=length(mmodale);

mFloq3mat=mFloq3.'*ones(1,Nmodi);
Q=Thba.'*(mFloq3mat.*Teab);
ind=1:Nmodi;
Q=Q+diag(mmodale);
Qinv=inv(Q);
Rmodmat=sqrt(mmodale.')*ones(1,Nmodi);
RFloqmat=sqrt(mFloq3.')*ones(1,Nmodi);
S21=2*(Rmodmat.*(Qinv*Thba.')).*RFloqmat.';
S12 = 2*(RFloqmat.*(Teab*Qinv)).*Rmodmat.';
S11 = 2*(RFloqmat.*(Teab*Qinv*Thba.')).*RFloqmat.';
S22 = 2*(Rmodmat.*Qinv).*Rmodmat.';

if strcmp(formulaz,'EFIE')
	S11(ind,ind) = -(S11(ind,ind)-eye(Nmodi));
	S22(ind,ind) = -(S22(ind,ind)-eye(Nmodi));
else % formulaz='HFIE'
	S11(ind,ind) = S11(ind,ind) - eye(Nmodi);
	S22(ind,ind) = S22(ind,ind) - eye(Nmodi);
end