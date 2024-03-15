ER_STA = 15.1-2.87*mesh.xmol_well+0.67*mesh.xmol_well^2;            % GaAs, Ioffe
ER_INF = 12.3-1.4*mesh.xmol_well;                  % GaAs, Ioffe
%
Nb = sqrt(ER_STA);               % refractive index
%
% mesh.gammak = 6e-3*qel/hbar;     % 1/s
% mesh.gammak = 10e-3*qel/hbar;     % 1/s
mesh.gammak = 1e13;     % 1/s
mesh.gammak = 1e13;     % 1/s
 mesh.tnm=150e-15;  %s  for non markovian effects
 mesh.tnm=100e-15;  %s  for non markovian effects
%mesh.tnm=220e-15;  %s  for non markovian effects
%mesh.tnm=40e-15;  %s  for non markovian effects
%
engy_LO = 34e-3;                 % InGaAs, Ioffe, eV
% Much more complicate expression for InGaAs, include? yes later
alpha_G=(6*(1-mesh.xmol_well)^2-8.6*(1-mesh.xmol_well)+5.2)*10^-4; 
beta_G=337*(1-mesh.xmol_well)^2-455*(1-mesh.xmol_well)+196; % coefficients for bandgap temperature shift

%Eg (0) + (6x2- 8.6x +5.2)·10-4 T2/(337x2- 455x +196)
% GaxIn1-xAs on Ga As Karachevtseva et al.(1994) inverted on ioffe so
% Determination of the Temperature Dependence of the Band Gap Energy of Semiconductors from Transmission Spectra
%JEAN WEI,1,3 JOEL M. MURRAY,1 JACOB BARNES,1 LEONEL P. GONZALEZ,2 and SHEKHAR GUHA2

