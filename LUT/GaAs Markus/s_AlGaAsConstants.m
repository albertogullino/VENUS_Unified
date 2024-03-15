XMOL=mesh.xmol_well;

ER_STA = 12.90-2.84*XMOL;                  % Al(x)GaAs, Ioffe
ER_INF = 10.89-2.73*XMOL;                  % Al(x)GaAs, Ioffe
%
Nb = sqrt(ER_STA);               % refractive index
%
engy_LO = (36.25+1.83*XMOL+17.12*XMOL^2-5.11*XMOL^3)*1e-3;  % Al(x)GaAs, Ioffe, eV
%
alpha_G=5.41e-4; beta_G=204; % coefficients for bandgap temperature shift (GaAs)
