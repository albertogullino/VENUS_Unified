ER_STA = 12.90;                  % GaAs, Ioffe
ER_INF = 10.89;                  % GaAs, Ioffe
%
Nb = sqrt(ER_STA);               % refractive index
%
% mesh.gammak = 6e-3*qel/hbar;     % 1/s
% mesh.gammak = 10e-3*qel/hbar;     % 1/s
mesh.gammak = 1e13;     % 1/s
% mesh.tnm=180e-15;  %s  for non markovian effects
%mesh.tnm=220e-15;  %s  for non markovian effects
mesh.tnm=40e-15;  %s  for non markovian effects
%
engy_LO = 35e-3;                 % GaAs, Ioffe, eV
%
alpha_G=5.41e-4; beta_G=204; % coefficients for bandgap temperature shift
