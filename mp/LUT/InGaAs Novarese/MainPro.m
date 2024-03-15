%==========================================================================
% Initialization
%==========================================================================
clear
clear global
close all
clc
% Select quantum well 1-> InGas/GaAsP 0-> GaAs/AlGaAs
Qw_sel=1;
%==========================================================================
% Quantum well parameters
%==========================================================================



if Qw_sel==0
    mesh_ALGaAs;
end

if Qw_sel==1
  mesh_InGaAsP;
end
mesh.ncinit=4; % number of conduction subbands to be computed
mesh.nvinit=6; % number of valence subbands to be computed
mesh.num_kvectors=282; % number of k points
mesh.max_k=0.14*2; % maximum k point, angstrom
mesh.parabolic=0; % if 1 the parabolic approximation is applied
%==========================================================================
% Optical response: settings
%==========================================================================
mode.iplot=1;
mode.fgain=0; % if equal to 1 computes gain by summing the gain from Upper and Lower matrix of 4x4 obtained from 8x8 Hamiltonian.
setfg=0;
fileex=0;
% NON TOCCARE !!!!!!
mode.LUT=0; % if 0 enables additional saves/plots
mode.flagLimitMaxDensity=0; % 1: limits max carrier density; 0: doesn't
mode.iline='lorentzian'; % 'nMark' ; 'lorentzian'; 'Landsberg'; Landsberg_extended' and Landsberg_modified  models 
mode.iline='nMark'; % 'nMark' ; 'lorentzian'; 'Landsberg'; Landsberg_extended' and Landsberg_modified  models 
mode.Expgammak=3; % exponent for carrier-carrier scattering corrections
mode.ieh_equal=1; % 1: e- and h- densities are assumed equal; 0: they aren't
mode.ifit=1; % 1: use a larger k grid by fitting it with spline; 0: don't
% NON TOCCARE !!!!!! FINO A QUI!!!!
%
mode.iren=0; % enable gap renormalization computation
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
if Qw_sel==0
lambdavet=linspace(800,1100,201)*1e-9;
else
    lambdavet=linspace(800,900,201)*1e-9;
end
Densityv=[5 10 20]*1e12;
 Densityv=[0.01 1 3 10]*1e12; % for PL plot
 
 Densityv=[4.]*1e12;
% Densityv=[10]*1e12;
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
    fileex=1;
else
if mode.fgain==1
mode.fgain=0;
setfg=1;
end
    [Ban,mesh]=f_ComputeQWSubbands(mesh,mode,Qw_sel);
end

%==========================================================================
% Da qui in poi non dovrebbe fregartene più nulla.
%==========================================================================
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
    
end% Temp
FirstLoopTime=toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reapeat all previous caluclations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mode.fgain=0;
if setfg==1 && fileex==0
mode.fgain=1;
end
if mode.fgain==1 
Gu=G;
[Ban,mesh]=f_ComputeQWSubbands(mesh,mode,Qw_sel);

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
itrans=[];
if(isempty(mode.vic) & isempty(mode.vic))
    mode.vic=1:mesh.ncb;
    mode.viv=1:mesh.nvb;
end
for indc=1:length(mode.vic)
    itrans=[itrans,mode.viv+(mode.vic(indc)-1)*mesh.nvb];
end
mode.ntrans=itrans;
if mode.ieh_equal~=1
    G=zeros(lP,lP,lL,lT);
    Rsp=zeros(lP,lP,lT);
else
    G=zeros(lP,lL,lT);
    Rsp=zeros(lP,lT);
end
Es=G;
Dep=G;

eDensityv=Densityv;
hDensityv=Densityv;
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
end% Temp
Gl=G;

G=Gl/2+Gu/2; % final gain 

end

lambda = lambdavet*1e9;
N=Densityv*1e-12;

%save LUT G Es Rsp N T lambda


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end of the second iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=figure
set(h,'pos',[  420    58   550   887])
subplot(211)
hold on
for k=1:length(eDensityv), 
plot(lambdavet*1e9,squeeze(G(k,:)),'linewidth',2), end
% ylim([-500, 1500])
% grid, pausak
global EFsave ECVtotsave
% figure,hold on,plot((E Fnsave-EFpsave))
EF=EFsave;
EGAP=ECVtotsave;
ax = gca;
ax.ColorOrderIndex = 1;
lambdaF=hbar*2*pi*c_light./(EF*qel)*1e9
lambdaGAP=2*pi*c_light./EGAP*1e9;
%for k=1:length(eDensityv), plot(lambdaF(k),0,'o','linewidth',2), end
%plot(lambdaGAP,zeros(size(lambdaF)),'ko')
ylabel('Gain 1/cm')
grid
title([' Tamb =',num2str(T)])
xlim([lambdavet([1 end])*1e9])
ax = gca;
ax.ColorOrderIndex = 1;
figure(h)
subplot(212)
hold on
for k=1:length(eDensityv), P=squeeze(Es(k,:)); Pl=P/max(P); plot(lambdavet*1e9,Pl,'linewidth',2), end
%axis([1.424 1.524 -500 1500])
ylabel('Es arb.u.')
xlim([lambdavet([1 end])*1e9])
xlabel('Wavelengh nm')
grid 

legend(num2str(Densityv'*1e-12))
