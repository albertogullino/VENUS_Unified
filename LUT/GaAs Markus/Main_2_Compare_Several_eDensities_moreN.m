%==========================================================================
% Initialization
%==========================================================================
clear
clear global
close all
clc
%==========================================================================
% Quantum well parameters
%==========================================================================
mesh.L=25e-9; % Total length of the domain, m
mesh.nn=251; % Number of spatial mesh nodes
mesh.Lz=7.7e-9; % quantum well width, m
mesh.xmol_barrier=0.286; % Al molar fraction (barrier)
mesh.xmol_well=0.0; % Al molar fraction (well)
mesh.Eg=1.412; % fitted from Gerlach PL measurements
mesh.DeltaEg=1.247; % DeltaEg
mesh.Qc=+0.62; % conduction band-offset percentage of DeltaEg
mesh.Delta=0.34; % spin-orbit coupling, eV
mesh.meffn=0.067; % electron conduction mass
mesh.ncinit=4; % number of conduction subbands to be computed
mesh.nvinit=6; % number of valence subbands to be computed
mesh.num_kvectors=282; % number of k points
mesh.max_k=0.14*2; % maximum k point, angstrom
mesh.parabolic=0; % if 1 the parabolic approximation is applied
%==========================================================================
% Optical response: settings
%==========================================================================
mode.iplot=0;
mode.LUT=0; % if 0 enables additional saves/plots
mode.flagLimitMaxDensity=0; % 1: limits max carrier density; 0: doesn't
mode.iline='lorentzian'; % 'nMark' ; 'lorentzian'; 'Landsberg'; Landsberg_extended' and Landsberg_modified  models 
mode.Expgammak=3; % exponent for carrier-carrier scattering corrections
mode.iren=1; % enable gap renormalization computation
mode.ieh_equal=1; % 1: e- and h- densities are assumed equal; 0: they aren't
mode.ifit=1; % 1: use a larger k grid by fitting it with spline; 0: don't
mode.vic=[]; % if empty, all transitions are computed
mode.viv=[]; % if empty, all transitions are computed
% % to compute the contributions to gain and spontaneous emission with
% % few transitions, use the following lines
% mode.vic=1;
% mode.viv=1;
%==========================================================================
% Optical response: parameters
%==========================================================================
Tvet=[300];
% lambdavet=linspace(700,900,101)*1e-9;
lambdavet=linspace(840,900,101)*1e-9;
%  Densityv=[ 0.01 5 10 20]*1e12;
% Densityv=[3 6 10]*1e12;
Densityv=[4 4.01]*1e12;
%Densityv=[10]*1e12;
%==========================================================================
%
%
%
%
%
%==========================================================================
% Computing or loading subbands
%==========================================================================
mesh.fileName=['subbands_WQW=',num2str(mesh.Lz*1e9),'nm_',num2str(mesh.num_kvectors),'k_',num2str(mesh.max_k),'kMax.mat'];
%
if(exist(mesh.fileName,'file'))
    disp('Loading existing subbands file')
    load(mesh.fileName)
else
    [Ban,mesh]=f_ComputeQWSubbands(mesh,mode);
end
%==========================================================================
% Loading constants
%==========================================================================
s_LoadConstants
s_GaAsConstants
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
if mode.ieh_equal~=1
    G=zeros(lP,lP,lL,lT);
    Rsp=zeros(lP,lP,lT);
else
    G=zeros(lP,lL,lT);
    Rsp=zeros(lP,lT);
end
Es=G;
Dep=G;
%
eDensityv=Densityv;
hDensityv=Densityv;
%==========================================================================
% First temperature loop: computing complex refractive index
%==========================================================================
dGdN=zeros(length(Tvet),length(lambdavet));
tic
for indT=1:length(Tvet) % temperature loop
    
    [GT,EsT,DepT,RspT,eD,hD,EFc,EFv,Ren]=pf_functionFit(eDensityv,hDensityv,lambdavet,Tvet,indT,Ban,mesh,mode,Pargain);
    
    if mode.ieh_equal~=1
        G(:,:,:,indT)=GT;
        Dep(:,:,:,indT)=DepT;
        Es(:,:,:,indT)=EsT;
        Rsp(:,:,indT)=RspT;
    else
        G(:,:,indT)=GT;
        Dep(:,:,indT)=DepT;
        Es(:,:,indT)=EsT;
        Rsp(:,indT)=RspT;
    end
    
    eD1(:,indT)=eD;
    hD1(:,indT)=hD;
    EFcv(:,indT)=EFc;
    EFvv(:,indT)=EFv;
    
    disp(['Loop 1, temperature ',num2str(indT),' of ',num2str(length(Tvet))])
    
%     dGdN(indT,:)=diff(GT)./diff(eDensityv);
end% Temp
FirstLoopTime=toc




col='bgmky';
%col='b'

% for k=1:length(eDensityv), plot(lambdavet*1e9,squeeze(G(k,:)),'linewidth',2), end
figure,plot(lambdavet*1e9,squeeze(G))
grid on

dGdN=diff(G)./diff(eDensityv);
figure,plot(lambdavet*1e9,dGdN)

return



