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
mesh.NQW=1;
mesh.L=25e-9; % Total length of the domain, m (25 nm, 1 QW) 
mesh.nn=251; % Number of spatial mesh nodes (251 pts, 1 QW)
mesh.Lz=7e-9; % quantum well width, m

mesh.xmol_barIn=0.33;   % Al molar fraction (inner barrier) - Al(x)GaAs(y)P
mesh.ymol_barIn=0.90;   % Al molar fraction (inner barrier) - Al(x)GaAs(y)P
mesh.xmol_barOut=0.45;  % Al molar fraction (outer barrier) - Al(x)GaAs
mesh.ymol_barOut=[];    % [] molar fraction (outer barrier) - Al(x)GaAs
mesh.xmol_well=0.767;   % Ga molar fraction (well) - Ga(x)InAs(y)P
mesh.ymol_well=0.604;   % As molar fraction (well) - Ga(x)InAs(y)P

mesh.Eg=1.412; % fitted from Gerlach PL measurements (GaAs)
mesh.DeltaEg=1.247; % DeltaEg = Eg(barrier) - Eg(well)
mesh.Qc=+0.62; % conduction band-offset percentage of DeltaEg
mesh.Delta=0.34; % spin-orbit coupling, eV
mesh.meffn=0.067; % electron conduction mass
mesh.ncinit=4*mesh.NQW; % number of conduction subbands to be computed
mesh.nvinit=6*mesh.NQW; % number of valence subbands to be computed
mesh.num_kvectors=282; % number of k points
mesh.max_k=0.14*2; % maximum k point, angstrom
mesh.parabolic=0; % if 1 the parabolic approximation is applied
%==========================================================================
% Optical response: settings
%==========================================================================
% mode.fileNameLUT=(['LUT4D_FBH_nMark_x=',num2str(mesh.xmol_well),'_y=',num2str(mesh.ymol_well)]); % filename for LUT
mode.fileNameLUT=(['LUT4D_FBH_nMark_x=',num2str(mesh.xmol_well)]); % filename for LUT
mode.LUT=1; % if 0 enables additional saves/plots
mode.DeltaDensity_Perc=1+0.005; % increment of carrier density for Jacobian derivatives
mode.iplot=1; % if 1, several intermediate plots are produced
mode.flagLimitMaxDensity=0; % 1: limits max carrier density; 0: doesn't
mode.iline='nMark'; % 'nMark' ; 'lorentzian'; 'Landsberg'; Landsberg_extended' and Landsberg_modified  models
mode.Expgammak=6; % exponent for carrier-carrier scattering corrections
% mesh.gammak = .90e13;     % 1/s
mesh.gammak=1.25e13;     % 1/s (Julian)
% mesh.tnm=70e-15;  %s  for non markovian effects
mesh.tnm=120e-15;  %s  for non markovian effects
%
mode.iren=0; % enable gap renormalization computation (0 faster, for first attempts)
mode.ieh_equal=0; % 1: e- and h- densities are assumed equal; 0: they aren't
mode.ifit=1; % 1: use a larger k grid by fitting it with spline; 0: don't
mode.vic=[]; % if empty, all transitions are computed
mode.viv=[]; % if empty, all transitions are computed
% to compute the contributions to gain and spontaneous emission with
% few transitions, use the following lines

% Full LUT parameters
Tvet=298;
Densityv=2*1e12;
lambdavet=794*1-9;%  % 780-825 for Al(0.04)GaAs, step 2 nm
%
%==========================================================================
% Region to investigate
%==========================================================================
vz_offset=0;    % 1 QW
% vz_offset=[-mesh.Lz,mesh.Lz]; % 2 QWs case
% vz_offset=[-(mesh.NQW-1)*mesh.Lz,0,(mesh.NQW-1)*mesh.Lz]; % 3 QWs case
%
mesh.vz_offset=vz_offset;
%
%==========================================================================
% Computing or loading subbands
%==========================================================================
%mesh.fileName=['subbands_WQW=',num2str(mesh.Lz*1e9),'nm_',num2str(mesh.num_kvectors),'k_',num2str(mesh.max_k),'kMax.mat'];
mesh.fileName=['subbands_WQW=',num2str(mesh.Lz*1e9),'nm_','x=',num2str(mesh.xmol_well),'y_=',num2str(mesh.xmol_well),'NQW=',num2str(mesh.NQW),'.mat'];
%
[Ban,mesh]=f_ComputeQWSubbands(mesh,mode);
%==========================================================================
% Loading constants
%==========================================================================
s_LoadConstants
s_AlGaAsConstants
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
if(isempty(mode.vic) && isempty(mode.vic))
    mode.vic=1:mesh.ncb;
    mode.viv=1:mesh.nvb;
end
for indc=1:length(mode.vic)
    itrans=[itrans,mode.viv+(mode.vic(indc)-1)*mesh.nvb];
end
mode.ntrans=itrans;
%
G=zeros(lP,lP,lL,lT);
Rsp=zeros(lP,lP,lT);
Es=G;
Dep=G;
%
eDensityv=Densityv;
hDensityv=Densityv;
%==========================================================================
% First temperature loop: computing complex refractive index
%==========================================================================

parfor indT=1:length(Tvet) % temperature loop
    
    [GT,EsT,DepT,RspT,eD,hD,EFc,EFv,Ren]=pf_functionFit(eDensityv,hDensityv,lambdavet,Tvet,indT,Ban,mesh,mode,Pargain);
    
    G(:,:,:,indT)=GT;
    Dep(:,:,:,indT)=DepT;
    Es(:,:,:,indT)=EsT;
    Rsp(:,:,indT)=RspT;
    
    eD1(:,indT)=eD;
    hD1(:,indT)=hD;
    EFcv(:,indT)=EFc;
    EFvv(:,indT)=EFv;
        
end


for indT=1:length(Tvet) % temperature loop
    EFcm=EFcv(ceil(end/2),indT);
    EFvm=EFvv(ceil(end/2),indT);
    [eDensity1,hDensity1,eDenV,hDenV] = f_charge_WFk(mesh,Ban,EFcm,EFvm,kBTev) ;
    
    Livh=Ban.SBV(:,1);
    
    cDS0=4*pi*m0*kBT/h^2*1e-4;
    for kh=1:length(Livh)
        DelE=(Livh(kh)-EFvm)/kBTev;
        DelEb=(-mesh.V0-EFvm)/kBTev;
        Lon=cDS0*(log(1+exp(-DelE))-log(1+exp(-DelEb)));
        mheff(kh)=1e-4*hDenV(kh)/Lon;
    end
    Mef_h(indT,1:length(Livh))=mheff;
end

eDensityv=Densityv*mode.DeltaDensity_Perc;
hDensityv=Densityv;

meshQW=mesh;
modeQW=mode;

port_2De=diag(eD1);
port_2Dh=diag(hD1);

h2m0=hbar^2/(2*m0);

clear MefH
for kh=1:Ban.nvh
    mdu=Mef_h(:,kh);
    finz=find(mdu>0);
    mval=mdu(finz);
    MefH(kh)=mean(mval);
end
mh=MefH';
SBVap=Ban.SBV(:,1)*ones(1,mesh.num_kvectors)+1./mh*(h2m0.*Ban.kgrid.^2/qel);

figure
plot(Ban.kgrid*1e-9,-Ban.SBV,'linewidth',1.5), chold,grid on
plot(Ban.kgrid*1e-9,-SBVap,'.','linewidth',2)
xlabel('k_{||}, nm^{-1}')
ylabel('Energy (VB)'),ylim([-0.2 0])


figure,grid on,box on,hold on
plot(mesh.x*1e7,abs(Ban.XVC./max(Ban.XVC)).^2,'.-')
xlabel('z, nm'),ylabel('Norm. Wavefunction (CB)'),ylim([0 1])

figure,grid on,box on,hold on
plot(mesh.x*1e7,abs(Ban.XVV./max(Ban.XVV)).^2,'.-')
xlabel('z, nm'),ylabel('Norm. Wavefunction (VB)'),ylim([0 1])

