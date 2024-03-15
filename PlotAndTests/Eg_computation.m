s_LoadConstants
heV=h/qel;
lambda=850e-9;
c=Clight/1e2;
E=heV*c/lambda

xmol=0.12

Tg=300;
alpha_G=5.41e-4; beta_G=204;
alpha_X=4.6e-4; beta_X=204;
alpha_L=6.05e-4; beta_L=204;


Eg0_G=1.519+1.155.*xmol+0.370.*xmol.^2;
Eg_G=Eg0_G-(alpha_G.*Tg.^2)./(beta_G+Tg);
%
Eg0_X=1.981+0.124.*xmol+0.144.*xmol.^2;
Eg_X=Eg0_X-(alpha_X.*Tg.^2)./(beta_X+Tg);
%
Eg0_L=1.815+0.69.*xmol; % mind the error in Ioffe!
Eg_L=Eg0_L-(alpha_L.*Tg.^2)./(beta_L+Tg);
%
%' qui Eg', keyboard
Eg = min(Eg_G,Eg_X)

if E-Eg>0
    fprintf('Absorption (Eg w/out BGN)\n')
end

% BGN
NA=8e19;
% NA=2e20;

% ND=1e19;
ND=2e19;
% ND=3e19;
% ND=4e19;

dEgA=3e-8.*(NA).^(1/3)
dEgD=3.5e-8.*(ND).^(1/3)

EgA=Eg-dEgA
EgD=Eg-dEgD

if E-EgA>0
    fprintf('Absorption (Eg with BGN - acceptor)\n')
end

if E-EgD>0
    fprintf('Absorption (Eg with BGN - donor)\n')
end