function r = msqrt(x)

% Funzione per il calcolo della radice 
% quadrata tale per cui la parte immaginaria
% della radice di un numero complesso
% e' negativa
% RIga di comando
%          r = msqrt(x);
% ------------------------------------------
% Matlab infatti calcola le radici quadrate
% scegliendo tra i due possibili risultati
% quello con parte immaginaria positiva
% ------------------------------------------
% Autore: Sergio Bastonero
% Torino, 23 gennaio 1998
% ------------------------------------------

r = sqrt(x);
indice = find(imag(r)>0);
r(indice) = - r(indice);
