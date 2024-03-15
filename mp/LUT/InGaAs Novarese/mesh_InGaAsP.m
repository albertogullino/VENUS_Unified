xmol_well=1-0.083; % 
%xmol_well=0.8; % 980 nm
xmol_barr= 1-0.37; % givenm data from reference

%% Mesh parameters InGaAs
gap_adjust

QW_strain_set
 %mesh.Lz=6*10^-9; % quantum well width, m too small due to strain compensation,
  %mesh.L=mesh.Lz+2*10*10^-9; % Total length of the domain, m AlGaAs
mesh.L=40*10^-9; % Total length of the domain, m AlGaAs 17 sigola barriera( si accoppia con altri well, completo)
mesh.Lz=WW*10^-9;
step=0.1*1e-9;
mesh.nn=round(mesh.L/step+1); % Number of spatial mesh nodesn
mesh.xmol_barrier=xmol_barr; % As molar fraction (barrier)
mesh.xmol_well=xmol_well; % Ga molar fraction (well)
mesh.Eg=gap; % 
mesh.C0=Delta_Ec;
mesh.V0=Delta_Ev; % minus if using chuang model
% 
%  strain=0;
%  strain_b=0;
mesh.Delta=well.so; % spin-orbit coupling, eV in the well
mesh.strain_w=strain;
mesh.strain_b=strain_b;
% effective mass, instead of this use the normal formula from chuang
mesh.meffn_w=well.mn;
mesh.meffn_b=barr.mn;
%%% full formula vangard, from chuang
% %InGaAs
% mesh.meffn_w=0.08-0.116*xmol_well+0.026-0.059*xmol_well+0.064-0.02*xmol_well+(0.06+0.032)*xmol_well^2;
% %GaAsP
% mesh.meffn_b=0.08-0.116+0.026*xmol_barr-0.059*xmol_barr+(0.064-0.02)*xmol_barr^2+0.06+0.032*xmol_barr;
% GaAsP(0.10), 30 nm
% In(0.27)GaAs, 6 nm
% GaAsP(0.10), 10 nm
% In(0.27)GaAs, 6 nm
% GaAsP(0.10), 10 nm
% In(0.27)GaAs, 6 nm
% GaAsP(0.10), 30 nm

mesh.BB=BB;
mesh.WW=WW;
mesh.NQW=NQW;

mesh.barr1QW=30*10^-9;
mesh.well1QW=WW*10^-9;

mesh.barr2QW=BB*10^-9;
mesh.well2QW=WW*10^-9;
mesh.barr2QW=BB*10^-9;
mesh.well3QW=WW*10^-9;
mesh.barr3QW=30*10^-9

mesh.L=mesh.barr1QW+mesh.well1QW+mesh.barr2QW+mesh.well2QW+mesh.barr2QW+mesh.well3QW+mesh.barr3QW;
step=0.1*1e-9;
mesh.nn=round(mesh.L/step)+1; % Number of spatial mesh nodesn









