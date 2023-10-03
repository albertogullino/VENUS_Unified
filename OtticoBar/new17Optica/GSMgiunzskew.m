function [S11, S12, S21, S22] = GSMgiunzskew(TeLR,ThRL, immittanceL, immittanceR)

%       R. Orta Dicembre 2010 
% apparentemente usa la versione ab del mode matching, pag 25 (tesi
% Bastonero)
% ----------------------------------------------------------
% Funzione per il calcolo della matrice scattering della
% giunzione tra due PSW a partire dalle matrici di proiezioni.
% ----------------------------------------------------------
%
NmodiTETM=length(immittanceR);

immittanceLmat=immittanceL.'*ones(1,NmodiTETM);
Q=ThRL*(immittanceLmat.*TeLR);
Q=Q+diag(immittanceR);
Qinv=inv(Q);

SqRmat=sqrt(immittanceR.')*ones(1,NmodiTETM);
SqLmat=sqrt(immittanceL.')*ones(1,NmodiTETM);
S21=2*(SqRmat.*(Qinv*ThRL)).*SqLmat.';
%S12 = 2*(SqLmat.*(TeLR*Qinv)).*SqRmat.';
S11 = 2*(SqLmat.*(TeLR*Qinv*ThRL)).*SqLmat.'-eye(NmodiTETM);
%S22 = 2*(SqRmat.*Qinv).*SqRmat.';

Q=ThRL*(immittanceLmat.*TeLR);
Q=diag(immittanceR)-Q;
S22=SqRmat.*(Qinv*Q)./SqRmat.';

S12=SqLmat.*(TeLR*(eye(NmodiTETM)+Qinv*Q))./SqRmat.';

