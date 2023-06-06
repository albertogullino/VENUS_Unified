addpath(genpath('Termico'))
addpath(genpath('dati'))
addpath(genpath('Gain'))
addpath(genpath('out'))
addpath(genpath('generageomBar'))

addpath(genpath('OtticoBar'))
rmpath('OtticoBar/new23OpticaGR')
rmpath('OtticoBar/old22Optica')
% rmpath('OtticoBar/new22Optica')

addpath(genpath('ExperimentalData'))
addpath(genpath('PlotAndTests'))
addpath(genpath('LUT'))

%% Path and input/output files definition
% Change the path for the geometry generation and for Optical simulation
% (keep inewpat = 1, by default)
Dir=cd;
datiOut='out';
datiIn='dati';
datiLUT='LUT';
nomeSW=[Dir,'\',datiOut,'\'];
nomeSR=[Dir,'\',datiIn,'\'];
nomeLUT=[Dir,'\',datiLUT,'\'];