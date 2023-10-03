function [S11, S12, S21, S22] = GSMgiunzS_U(Teab,Thba, mFloq4, mmodale, formulaz)

%
% ----------------------------------------------------------
% Riga di comando
%               [S11, S12, S21, S22] = GSMgiunzS_U(P, mi1, mpn, HE);
% ----------------------------------------------------------
% Funzione per il calcolo della matrice scattering della
% giunzione reg.stratif.periodica-reg.unif.
% ----------------------------------------------------------
% Variabili d'ingresso:
% P             matrice di proiezione tra i modi della regione unif.
%               e i modi della regione stratif.
% mi1, mpn      vettori delle ammettenze o delle impedenze (a
%               seconda della formulazione desiderata) delle due
%               regioni
% HE    variabile che distingue una formulazione dall'altra
%                       0 --> HFIE
%                       1 --> EFIE
%
% Variabili d'uscita:
% S11, S12, S21, S22 elementi della matrice scattering
%               della giunzione
% ----------------------------------------------------------
%
Nmodi=length(mmodale);

mmodalemat=mmodale.'*ones(1,Nmodi);
%Q=Teab.'*(mmodalemat.*Thba); %errato??
Q=Thba*(mmodalemat.*Teab.');
ind=1:Nmodi;
Q=Q+diag(mFloq4);
Qinv=inv(Q);
Rmodmat=sqrt(mmodale.')*ones(1,Nmodi);
RFloqmat=sqrt(mFloq4.')*ones(1,Nmodi);
S21=2*(RFloqmat.*(Qinv*Thba)).*Rmodmat.';
%S12 = 2*(Rmodmat.*(Teab.'*Qinv)).*RFloqmat.';
S12 = 2*(Rmodmat.*(Teab.'*Qinv)).*RFloqmat.';
%S11 = 2*(Rmodmat.*(Teab.'*Qinv*Thba)).*Rmodmat.';
S11 = 2*(Rmodmat.*(Teab.'*Qinv*Thba)).*Rmodmat.';
S22 = 2*(RFloqmat.*Qinv).*RFloqmat.';

if strcmp(formulaz,'EFIE')
	S11(ind,ind) = -(S11(ind,ind)-eye(Nmodi));
	S22(ind,ind) = -(S22(ind,ind)-eye(Nmodi));
else % formulaz='HFIE'
	S11(ind,ind) = S11(ind,ind) - eye(Nmodi);
	S22(ind,ind) = S22(ind,ind) - eye(Nmodi);
end

