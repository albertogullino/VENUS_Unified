%% Adjust Gap value value of InAs and GaP
xwell=xmol_well;
xbarr=xmol_barr; % given data from reference; % givenm data from reference
%starting point
refEgGaAs=1.412;
refEgInAs=0.354;
refEgAlAs= 3.03 ;
%% well
x=xwell;
Egi=linspace(0,refEgInAs*2,1000);
Egg=linspace(0,refEgGaAs*2,1000);
EgQ= 0.36 +0.629*x + 0.426*x^2; 
EgL= x*Egg+ (1-x)*Egi; 
[val pos_E]=min(abs(EgQ-EgL));
mesh.EgInAs=Egi(pos_E);
mesh.EgGaAs=Egg(pos_E);
% %% barr no adjustement for AlGaAs it is already linear
% x=xbarr;
mesh.EgAlAs=1.424 + 1.247*x ;
%mesh.EgAlAs=2.6590;
% 