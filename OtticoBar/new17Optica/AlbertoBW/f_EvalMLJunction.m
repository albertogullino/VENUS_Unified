function [S11,S21,S12,S22]=f_EvalMLJunction(ZL,ZR)

S11=diag((ZR-ZL)./(ZL+ZR));
S22=-S11;
%-- Fondamentale usare sqrt(ZL).*sqrt(ZR), per i problemi di segno della
%radice; in questo modo il coefficiente di trasmissione ha segno giusto e
%siamo tutti felici.
S21=diag(2*sqrt(ZL).*sqrt(ZR)./(ZL+ZR));
S12=S21;