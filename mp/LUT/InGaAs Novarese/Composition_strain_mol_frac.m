function [Thick_barr, strain_well,strain_barrier]= Composition_strain_mol_frac(xmol_InGaAs,Thick_well,xmol_barr);
%=================================================================
%% costants
GapGaAs=1.421; %eV
aGaAs=5.65325; %A3
V_well=Thick_well;   %nm^
%% Molar fraction InGaAs vs GaAs substrate lattice costant In(1-x)Ga(x)
%http://www.ioffe.ru/SVA/NSM/Semicond/GaInAsP/basic.html
% a~= 5.8688-0.4176x+0.1896y+0.0125xy  for GaxIn1-xAsyP1-y so by putting
% y=1-> InGaAs x=1->GaAsP lattice costants
% x=linspace(0,1,10000);
x=xmol_InGaAs;

% Lamd=1.2408*Eg.^(-1); % um
% [v1 pos_min_InGaAs ]=min(abs(Lamd-Emis_lambd));
% xmol_InGaAs=x(pos_min_InGaAs);

%% Strain and GaAsP
% with the known molar fraction we can calculate the 
%lattice costant of  InGaAs
a_InGaAs=5.8688-0.4176*xmol_InGaAs+0.1896+0.0125*xmol_InGaAs;
strain_well=(-a_InGaAs+aGaAs)/a_InGaAs;



a_AlGaAs=5.6533+0.0078*xmol_barr;% check if it becomes 
%indirect! approximately at x=50%
%now we see at which point we have the same strain for the GaAsP (opposite
%sign, reason for the minus sign in the equation)

strain_barrier=(-a_AlGaAs+aGaAs)/a_AlGaAs; % takes into account the difference in volume of the wells with repsect to the barrier well*t1+barr*t2=0
%strain_barrier=0; %980nm GaAs

%ekins-daukes2002_strain
%==========================================================================
%strain compensation
%==========================================================================
Thick_barr=-Thick_well*strain_well/strain_barrier;
end

% figure
% plot(x,Lamd,'linewidth',2)
% xlabel('Molar fraction Ga');
% ylabel('\lambda  [ \mu m]');
% hold on
% plot(x,Emis_lambd*ones(1,1000),'--','linewidth',2);
% grid on
% box on
