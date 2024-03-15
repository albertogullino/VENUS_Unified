%% Mesh parameters ALGaAs
mesh.L=25e-9; % Total length of the domain, m AlGaAs
mesh.Lz=7.7e-9; % quantum well width, m
% mesh.xmol_barrier=0.286; % Al molar fraction (barrier)
mesh.xmol_barrier=0.21; % Al molar fraction (barrier)
mesh.xmol_well=0.0; % Al molar fraction (well)
mesh.Eg=1.412; % fitted from Gerlach PL measurements
mesh.DeltaEg=1.247; % DeltaEg from pag 666 chuang
mesh.Qc=+0.62; % conduction band-offset percentage of DeltaEg
mesh.Delta=0.34; % spin-orbit coupling, eV
mesh.meffn=0.067; % electron conduction mass
mesh.meffn_w=0.067;
mesh.meffn_b=0.067+(.124-.067)*mesh.xmol_barrier;
%  mesh.L=100*10^-9; % Total length of the domain, m AlGaAs
%  mesh.Lz=mesh.L*7.7/25;
 mesh.nn=round(251*mesh.L/25); % Number of spatial mesh nodesn equal spatial discretization