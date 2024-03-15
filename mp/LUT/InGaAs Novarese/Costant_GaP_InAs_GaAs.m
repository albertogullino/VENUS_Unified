function [GaAs,InAs,AlAs]=Costant_GaP_InAs_GaAs(mesh);

%costants from zwang
%GaAs note that the AlAs used here (Ec-Ev) gives 1.42 and not the one from
%Gerlach  used in the main code, al
bin_number=3;
GaAs.Ev_av=-6.92;  %eV
GaAs.Ev=0.111; %eV
GaAs.Ec=1.53; % eV
GaAs.so=0.34; %eV
GaAs.C11=11.879; % 10^11 dyn/cm^2
GaAs.C12=5.376; % 10^11 dyn/cm^2
GaAs.ac=-7.17; %eV
GaAs.av=1.16; %eV
GaAs.b= -1.7; %eV
GaAs.mn=0.067; %effective mass
GaAs.Eg=mesh.EgGaAs;
% GaAs.Eg=1.412;
% GaAs.Eg=1.311647647647648;
% %GaAs.Eg= 1.337089089089089;%980 nm
% GaAs.Eg=1.297513513513513; % 0.31 xmol
%InAs
InAs.Ev_av=-6.67;
InAs.Ev=0.441;
InAs.Ec=0.81;
InAs.so=0.39;
InAs.C11=8.329;
InAs.C12=4.526;
InAs.ac=-5.08;
InAs.av=1;
InAs.b=-1.8;
InAs.mn= 0.0224;
InAs.Eg=mesh.EgInAs;
% InAs.Eg=0.354;
% InAs.Eg=0.328840840840841;
% %InAs.Eg=0.335219219219219; %980 nm
% InAs.Eg=0.325297297297297; % 0.31 xmol
% AlAs
AlAs.Ev_av=-7.49;
AlAs.Ev=-0.388;
AlAs.Ec=2.5225;
AlAs.so=0.28;
AlAs.C11=5.34;
AlAs.C12=12.05;
AlAs.ac=-5.64;
AlAs.av=2.47;
AlAs.b=-1.5;
AlAs.mn= 0.124;
AlAs.Eg=mesh.EgAlAs;
% AlAs.Eg=2.776;
% AlAs.Eg=2.706530530530530;
%AlAs.Eg=2.706530530530530; % 0.31 xmol
end