%==========================================================================
% Initialization
%==========================================================================
clear
clear global
close all
clc

global NQW
Qw_sel=1;
%==========================================================================
% Quantum well parameters
%==========================================================================

NQW=3;

BBoxSide=31;  % spessore barriera
BB=4;  % spessore barriera
WW=4;  % spessore barriera

%==========================================================================
% Quantum well parameters
%==========================================================================



if Qw_sel==0
    mesh_ALGaAs;
end

if Qw_sel==1
  mesh_InGaAsP_GEN;
end

%==========================================================================
% Optical response: settings
%==========================================================================
mode.fileNameLUT='LUT5QW-Bar4nm_Jun19'; % filename for LUT
%mode.fileNameLUT=['PRO',num2str(BB),'_Jun19']; % filename for LUT
mode.LUT=1; % if 0 enables additional saves/plots
mode.DeltaDensity_Perc=1+0.05; % increment of carrier density for Jacobian derivatives
mode.iplot=0; % if 1, several intermediate plots are produced
mode.fgain=0; % if equal to 1 computes gain by summing the gain from Upper and Lower matrix of 4x4 obtained from 8x8 Hamiltonian.
mode.flagLimitMaxDensity=0; % 1: limits max carrier density; 0: doesn't
mode.iline='nMark'; % 'nMark' ; 'lorentzian'; 'Landsberg'; Landsberg_extended' and Landsberg_modified  models 
mode.Expgammak=6; % exponent for carrier-carrier scattering corrections
mode.iren=1; % enable gap renormalization computation
mode.ieh_equal=1; % 1: e- and h- densities are assumed equal; 0: they aren't
mode.ifit=1; % 1: use a larger k grid by fitting it with spline; 0: don't
mode.vic=[]; % if empty, all transitions are computed
mode.viv=[]; % if empty, all transitions are computed
% % to compute the contributions to gain and spontaneous emission with
mesh.ncinit=5; % number of conduction subbands to be computed
mesh.nvinit=18; % number of valence subbands to be computed
mesh.num_kvectors=22; % number of k points
mesh.max_k=0.09*2; % maximum k point, angstrom
mesh.parabolic=0; % if 1 the parabolic approximation is applied

% % few transitions, use the following lines
% % mode.vic=1;
% % mode.viv=1:2;
%==========================================================================
% Optical response: parameters
%==========================================================================
% Tvet=[300 350];
% Densityv=[ 0.01 5 10]*1e12;
 Tvet=[290:10:520];
 lambdavet=linspace(800,900,101)*1e-9;
 Densityv=linspace(0,20,31); % for PL plot
 
 Densityv(1)=4;
 Densityv=Densityv*1e12;

Tvet=290;
lambdavet=850*1e-9;

%lambdavet=[800:1:900]*1e-9;% linspace(830,870,41)*1e-9;
%==========================================================================
%
%
%
%
%
%==========================================================================
% Computing or loading subbands
%==========================================================================
mesh.fileName=['subbands_',num2str(NQW),'QW-w=',num2str(mesh.Lz*1e9),'-',num2str(BB),'Bar-',num2str(mesh.num_kvectors),'k_',num2str(mesh.max_k),'kMax.mat'];
%
if(exist(mesh.fileName,'file'))
    disp('Loading existing subbands file')
    load(mesh.fileName)
else
    [Ban,mesh]=f_ComputeQWSubbands(mesh,mode,Qw_sel);
end
%==========================================================================
% Loading constants
%==========================================================================
s_LoadConstants
if Qw_sel==0
s_GaAsConstants
else
 s_InGaAsConstants
end

%==========================================================================
% Program: initialization of constants
%==========================================================================
omegavet=2*pi*c_light./lambdavet;
Evet=omegavet*hbar/qel;
%
Ban=f_Refine_kgrid(Ban,mesh,mode); % refinement of k grid
%
lP=length(Densityv);
lT=length(Tvet);
lL=length(lambdavet);
%
EFcv=zeros(lP,lT);
EFvv=zeros(lP,lT);
eD1=zeros(lP,lT);
hD1=zeros(lP,lT);
%
Cost_Rsp=pi/(pi^2*(c_light/Nb).^3)*1e-6;
%
%==========================================================================
% Preliminary temperature loop: initialization of variables
%==========================================================================
for indT=1:length(Tvet)
    T=Tvet(indT);
    
    kBT=kB*T;
    kBTev=kBT/qel;
    
    DeltaE_Temp=(alpha_G.*T.^2)./(beta_G+T)-(alpha_G.*300.^2)./(beta_G+300);
    ECV_tot0=qel*(Ban.ECVf-DeltaE_Temp)/hbar;
    
    Pargain(indT).kBT=kBT;
    Pargain(indT).kBTev=kBTev;
    Pargain(indT).DeltaE_Temp=DeltaE_Temp;

    Pargain(indT).M2d=Ban.M2df;
    Pargain(indT).M2esd=Ban.M2esdf;
    Pargain(indT).ECV_tot0=ECV_tot0;
    Pargain(indT).Cost_Rsp=Cost_Rsp;
    
end
%
itrans=[];
if(isempty(mode.vic) & isempty(mode.vic))
    mode.vic=1:mesh.ncb;
    mode.viv=1:mesh.nvb;
end
for indc=1:length(mode.vic)
    itrans=[itrans,mode.viv+(mode.vic(indc)-1)*mesh.nvb];
end
mode.ntrans=itrans;
%
G=zeros(lP,lL,lT);
Rsp=zeros(lP,lT);
Es=G;
Dep=G;
%
eDensityv=NQW*Densityv;
hDensityv=NQW*Densityv;
%==========================================================================
% First temperature loop: computing complex refractive index
%==========================================================================
tic



%eval(FOR)
for indT=1:length(Tvet) % temperature loop
    
    [GT,EsT,DepT,RspT,eD,hD,EFc,EFv,Ren]=pf_functionFit(eDensityv,hDensityv,lambdavet,Tvet,indT,Ban,mesh,mode,Pargain);
    G(:,:,indT)=GT/NQW;
    Dep(:,:,indT)=DepT/NQW;
    Es(:,:,indT)=EsT/NQW;
    Rsp(:,indT)=RspT/NQW;
    
    eD1(:,indT)=eD/NQW;
    hD1(:,indT)=hD/NQW;
    EFcv(:,indT)=EFc;
    EFvv(:,indT)=EFv;
    
    disp(['Loop 1, temperature ',num2str(indT),' of ',num2str(length(Tvet))])
    
end% Temp

lambda=lambdavet*1e9;
N=Densityv;

 h=figure;
 set(h,'pos',[120         415        1796         525])
 subplot(131)
 plot(lambda,squeeze(Es(:,:,1)))
 xlabel('Wavelenght nm')
 ylabel('Es 1/s')

 subplot(133)
 plot(N,Rsp)
  xlabel('Carrier Density') 
  ylabel('R_{sp} 1/(cm^3 s)')
  
 subplot(132)
  plot(lambda,squeeze(G(:,:,1)))
   xlabel('Wavelenght nm')

    ylabel('Gain 1/cm ')
  
eval(['save ',mode.fileNameLUT,' lambda N Tvet G Rsp Es Dep'])