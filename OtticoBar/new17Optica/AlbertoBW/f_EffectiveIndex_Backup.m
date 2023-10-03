
function [neff,nBW,flag] = f_EffectiveIndex(Period,lambda,DutyCycle,nin,nout,n1,n2,thickness,NModes)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Geometry of the structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-- Global structure parameters
Geometry.d=Period; % nanometers
Geometry.n_in=nin; %-- Refractive index left half-space
Geometry.th_in=0; %-- Thickness (z axis) of the left half-space
Geometry.DC=DutyCycle; %-- Duty cycle of the bar
Geometry.n1=[n1];
Geometry.n2=[n2];
Geometry.Displacement=[0]; %-- Displacement of the grating
Geometry.RollOff=[0]; %-- Roll-off of the raised-cosine profile
Geometry.th=[thickness]; %-- Thickness of the tooth (z axis)
Geometry.n_out=nout; %-- Refractive index right (outer) half-space
Geometry.th_out=0; %-- Thickness (z axis) of the right half-space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options of the program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Options.RebuildProfile=0; %-- If 1, the layer "n" profile is re-built
Options.PlotGeometry=0; %-- If 1, the geometry is plot
Options.PlotFieldJunctions=0; %-- If 1, junction fields are plotted
Options.PlotTotalField=0; %-- If 1, 3-D field plots are produced
Options.PlotForcedHarmonicContent=0; %-- If 1, 3-D field plots are produced
Options.PlotModeHarmonicContent=0; %-- If 1, 3-D field plots are produced
Options.VerifyPower=0; %-- If 1, the power conservation is verified
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters of the program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Parameters.theta=0;
Parameters.phi=0;
Parameters.vLam=lambda; % nanometers
Parameters.NFFT=2^11;
%-- ALWAYS ODD numbers (generation of Floq. modes)
Parameters.NModes=NModes; %-- inner hemispace
Parameters.NHarmonicsTE=2*Parameters.NModes+1; %-- Harmonics used to build TE PSWW modes
Parameters.NHarmonicsTM=2*Parameters.NModes+1; %-- Harmonics used to build TM PSWW modes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Programma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta=Parameters.theta*pi/180;
phi=Parameters.phi*pi/180;

Nlam=length(Parameters.vLam);
vLambda=logspace(0,3,101);

RTE0=zeros(length(vLambda),Nlam);
RTM0=zeros(length(vLambda),Nlam);

indMultiModeTE=[];
indMultiModeTM=[];

DC=Geometry.DC;
n1=Geometry.n1;
n2=Geometry.n2;

% Born-Wolf 
nTEBornWolf=sqrt(n1.^2.*DC+n2.^2.*(1-DC));
nTMBornWolf=sqrt(1./(DC./(n1.^2)+(1-DC)./(n2.^2)));

k0=2*pi/lambda;
n_in=Geometry.n_in;

ky=k0*n_in*sin(theta)*sin(phi);
kx0=k0*n_in*sin(theta)*cos(phi); %-- kx for Floq. harmonic 0 = KBd/d

[S11,S21,S12,S22,JunctionInfo,LayerInfo,HalfSpaceInfo]=f_EvalSMatrixRCWA(Parameters.vLam,Parameters,Geometry,Options,kx0,ky);

S11Tot(:,:)=S11;
S12Tot(:,:)=S12;
S21Tot(:,:)=S21;
S22Tot(:,:)=S22;

kzTE=JunctionInfo(1).kzL(1:Parameters.NModes);
kzTM=JunctionInfo(1).kzL(Parameters.NModes+1:end);
indPropTE=find(abs(imag(kzTE))<=1e-8);
indPropTM=find(abs(imag(kzTM))<=1e-8);

RTE0=S11(1,1);
RTM0=S11(Parameters.NModes+1,Parameters.NModes+1);
flag=1;
if(length(indPropTE)>1)
    flag=0;
end

n2TE=fzero(@(X)f_WeightFunctionEffectiveMedium(X,RTE0,Geometry.n_in,Geometry.n_out,Geometry.th,Parameters.vLam,Parameters.theta*pi/180,'TE'),nTEBornWolf(1));
n2TM=fzero(@(X)f_WeightFunctionEffectiveMedium(X,RTM0,Geometry.n_in,Geometry.n_out,Geometry.th,Parameters.vLam,Parameters.theta*pi/180,'TM'),nTMBornWolf(1));
% n2TE=fzero(@(X)f_WeightFunctionEffectiveMedium(X,RTE0,Geometry.n_in,Geometry.n_out,Geometry.th,Parameters.vLam,Parameters.theta*pi/180,'TE'),[n2,n1]);
% n2TM=fzero(@(X)f_WeightFunctionEffectiveMedium(X,RTM0,Geometry.n_in,Geometry.n_out,Geometry.th,Parameters.vLam,Parameters.theta*pi/180,'TM'),[n2,n1]);

neff=[n2TE,n2TM];
nBW=[nTEBornWolf,nTMBornWolf];

return