%==========================================================================
% Initialization
%==========================================================================
clear
clear global
close all
colordef white
clc

global NQW
Qw_sel=1;
%==========================================================================
% Quantum well parameters
%==========================================================================

NQW=1;

BBoxSide=30;  % spessore barriera
BB=10;  % spessore barriera
WW=6;  % spessore QW

xW=0.27;

xW=0.318;
xB=0.1;

%==========================================================================
% Quantum well parameters
%==========================================================================



if Qw_sel==0
    mesh_ALGaAs;
end

if Qw_sel==1
  mesh_InGaAsP_GEN1;
end

%==========================================================================
% Optical response: settings
%==========================================================================
mode.fileNameLUT='Anders1060-6nmQw318In-Bar10nm_17Dic19SpRen'; % filename for LUT
%mode.fileNameLUT='Pippo'; % filename for LUT
%mode.fileNameLUT=['PRO',num2str(BB),'_Jun19']; % filename for LUT
mode.LUT=1; % if 0 enables additional saves/plots
mode.DeltaDensity_Perc=1+0.05; % increment of carrier density for Jacobian derivatives
mode.iplot=0; % if 1, several intermediate plots are produced
mode.fgain=0; % if equal to 1 computes gain by summing the gain from Upper and Lower matrix of 4x4 obtained from 8x8 Hamiltonian.
mode.flagLimitMaxDensity=0; % 1: limits max carrier density; 0: doesn't
mode.iline='nMark'; % 'nMark' ; 'lorentzian'; 'Landsberg'; Landsberg_extended' and Landsberg_modified  models 
%mode.iline='lorentzian'; % 'nMark' ; 'lorentzian'; 'Landsberg'; Landsberg_extended' and Landsberg_modified  models 
mode.Expgammak=6; % exponent for carrier-carrier scattering corrections
mode.iren=1; % enable gap renormalization computation
mode.ieh_equal=1; % 1: e- and h- densities are assumed equal; 0: they aren't
mode.ifit=1; % 1: use a larger k grid by fitting it with spline; 0: don't
mode.vic=[]; % if empty, all transitions are computed
mode.viv=[]; % if empty, all transitions are computed
% % to compute the contributions to gain and spontaneous emission with
mesh.ncinit=5; % number of conduction subbands to be computed
mesh.nvinit=18; % number of valence subbands to be computed
mesh.num_kvectors=92; % number of k points
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
 Tvet=[270:10:450];
 
 %Tvet=[270:100:450];
 lambdavet=linspace(800,900,101)*1e-9;
 Densityv=linspace(0,5,26); % for PL plot
 %Densityv=linspace(0,5,4); % for PL plot
 Densityv(1)=1e-4;
 Densityv=Densityv*1e12;
 %Densityv=[1e-5 1e-4 1e-3 0.01 0.1  0.3:.3:3 4:12]*1e12;
 
 Densityv=[.01 1:5];

 %Densityv=linspace(0,5,26); % for PL plot
 %Densityv=linspace(0,5,4); % for PL plot
 %Densityv(1)=1e-4;

%clear Densityv
% Densityv(1)=1e-4;

Dinc=1.1;




 Densityv=Densityv*1e12;
 Tvet=[300 330 360 390];
 
 Tvet=[300 350];
 Tvet=[350];
% Tvet=[290:2:330];
 %Tvet=[325];
% Tvet=[350];
%lambdavet=850*1e-9;

lambdavet=[1000:1:1100]*1e-9;% linspace(830,870,41)*1e-9;
%lambdavet=[1045:5:1065]*1e-9;% linspace(830,870,41)*1e-9;
lambdavet=[1045:2:1075]*1e-9;% linspace(830,870,41)*1e-9;

%lambdavet=[1060]*1e-9;% linspace(830,870,41)*1e-9;
%lambdavet=[1060]*1e-9;% linspace(830,870,41)*1e-9;

%==========================================================================
%
%
%
%
%
%==========================================================================
% Computing or loading subbands
%==========================================================================
mesh.fileName=['subbands_',num2str(NQW),'QW-w=',num2str(xW),'In%-',num2str(mesh.Lz*1e9),'-',num2str(BB),'Bar-',num2str(mesh.num_kvectors),'k_',num2str(mesh.max_k),'kMax.mat'];
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
'qui', keyboard
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

Gi=G;
Rspi=Rsp;
Esi=Es;
Depi=Dep;
%
eDensityv=NQW*Densityv;
hDensityv=NQW*Densityv;

eDensityvi=eDensityv*Dinc;
hDensityvi=hDensityv*Dinc;
%==========================================================================
% First temperature loop: computing complex refractive index
%==========================================================================
tic



%eval(FOR)
for indT=1:length(Tvet) % temperature loop
    
    [GT,EsT,DepT,RspT,eD,hD,EFc,EFv,Ren,EFsplit]=pf_functionFit(eDensityv,hDensityv,lambdavet,Tvet,indT,Ban,mesh,mode,Pargain);
    G(:,:,indT)=GT/NQW;
    Dep(:,:,indT)=DepT/NQW;
    Es(:,:,indT)=EsT/NQW;
    Rsp(:,indT)=RspT/NQW;
    
    eD1(:,indT)=eD/NQW;
    hD1(:,indT)=hD/NQW;
    EFcv(:,indT)=EFc;
    EFvv(:,indT)=EFv;
    EFspV(:,indT)=EFsplit;

    [GT,EsT,DepT,RspT,eD,hD,EFc,EFv,Ren,EFsplit]=pf_functionFit(eDensityvi,hDensityvi,lambdavet,Tvet,indT,Ban,mesh,mode,Pargain);
    Gi(:,:,indT)=GT/NQW;
    Depi(:,:,indT)=DepT/NQW;
    Esi(:,:,indT)=EsT/NQW;
    Rspi(:,indT)=RspT/NQW;
    
    eD1i(:,indT)=eD/NQW;
    hD1i(:,indT)=hD/NQW;
    EFcvi(:,indT)=EFc;
    EFvvi(:,indT)=EFv;
    EFspVi(:,indT)=EFsplit;
    
    
    disp(['Loop 1, temperature ',num2str(indT),' of ',num2str(length(Tvet))])
    
end% Temp

lambda=lambdavet*1e9;
N=Densityv;


 

 h=figure;

eval(['save ',mode.fileNameLUT,' lambda N Tvet G Gi Rsp Es Dep Depi EFspV'])
lambda=Densityv;
if length(size(G))==3
  plot(lambda,G,'.-')
   if length(lambdavet)>length(Densityv)
    xlabel('Wavelenght nm')
   else
    xlabel('Carrier Dens')
   end
else   
  plot(lambda,G)
  xlabel('Carrier Dens')

end
    ylabel('Gain 1/cm ')
if length(Tvet)>1    
legend(num2str(Tvet'))
else
legend(num2str(1e9*lambdavet'))
end

G=squeeze(G);
Gi=squeeze(Gi);
Dep=squeeze(Dep);
Depi=squeeze(Depi);
 
LaMat=ones(length(N),1)*(.02*pi./lambdavet);    
gg=(G-Gi)./LaMat;  
nn=squeeze(Dep);
nn=squeeze(Dep-Depi);

al=nn./gg;

figure, plot(1e9*lambdavet,al)
%ylim([-10 0])
ylabel('Henry factor')
xlabel('Wavelength (nm)')  
if length(Tvet)>1    
legend(num2str(Tvet'))
else
legend(num2str(1e-12*Densityv'))
end  
ylim([-2 0.1])
  
