function [S11,S12,S21,S22] = f_CascS(S11p,S12p,S21p,S22p,S11s,S12s,S21s,S22s,Kz,l,AttAccdB)

%------------------------------------------------------
% Calcola la matrice scattering della cascata di due
% matrici date in ingresso
%
% Sintassi:
% 	[S11,S12,S21,S22] =	f_CascS(S11p,S12p,S21p,S22p,S11s,S12s,S21s,...
%						S22s,Kz,l,AttAccdB2)

% Input : S**p = parametri scattering della prima guida
%			 S**s = parametri scattering della sec. guida
%			 Kz = costante di propagazione nella linea tra
%						le due matrici [mm^-1]
%			 l = lunghezza della linea tra le due matrici
%			 AttAccdB = attenuazione massima per i modi
%						considerati		
% Output : Parametri della matrice scattering cascata 
%						delle due in ingresso
%------------------------------------------------------

% Ricerca dei modi con attenuazione minore della soglia

D = exp(-j*Kz*l);
I = find(20*log10(abs(D))>=-AttAccdB); 
NAcc = length(I);
D = D(I);
D = diag(D);
if(l==0)
   I = ':';
   D = eye(length(Kz));
end

S12p = S12p(:,I);
S21p = S21p(I,:);
S22p = S22p(I,I);

S11s = S11s(I,I);
S12s = S12s(I,:);
S21s = S21s(:,I);


Kz = Kz(I);

% Traslazione della  prima matrice 

S22p = D*S22p*D;
S21p = D*S21p;
S12p = S12p*D;

% Cascata delle matrici scattering

ML=inv(eye(size(S22p*S11s))-S22p*S11s);

S11=S11p+S12p*S11s*ML*S21p;

S21=S21s*ML*S21p;

S22=S22s+S21s*ML*S22p*S12s;

ML=inv(eye(size(S11s*S22p))-S11s*S22p);

S12=S12p*ML*S12s;

return