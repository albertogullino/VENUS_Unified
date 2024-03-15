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

NQW=5;

BBoxSide=30;  % spessore barriera
BB=6;  % spessore barriera
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
mode.fileNameLUT='LUT5QW-Bar6nm_Jun19'; % filename for LUT
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
 Tvet=[290:10:520];
 lambdavet=linspace(800,900,101)*1e-9;
 Densityv=linspace(0,20,31); % for PL plot
 Densityv(1)=.001;
 Densityv=Densityv*1e12;

%lambdavet=1060*1e-9;

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


figure, plot(lambdavet*1e9,G), pausak

return

hh=figure
set(hh,'pos',[  420    58   550   887])
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
figure(hh)
subplot(212)
hold on
for k=1:length(eDensityv), P=squeeze(Es(k,:)); Pl=P/max(P); plot(lambdavet*1e9,Pl,'linewidth',2), end
%axis([1.424 1.524 -500 1500])
ylabel('Es arb.u.')
xlim([lambdavet([1 end])*1e9])
xlabel('Wavelengh nm')
grid 

legend(num2str(Densityv'*1e-12))
pausak



if mode.iren==0
openfig('g1qw0.fig')
else
%openfig('g1qw.fig')
%openfig('g1qwSN.fig')
%openfig('g1qwSNew.fig')
openfig('g1qwS4.fig')
end
hold on
plot(lambdavet*1e9,squeeze(G(1,1,:)),'g.','markersize',15), 
xlabel('Wavelength, nm')
ylabel('Gain, 1/cm @4e12/cm^2')
pausak

return

FirstLoopTime=toc

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


GE=zeros(lP,lP,lL,lT);
RspE=zeros(lP,lP,lT);
EsE=GE;
DepE=GE;   

eDensityv=NQW*Densityv*mode.DeltaDensity_Perc;
hDensityv=NQW*Densityv;
   
tic
% Calcolo G(El,El)  e G(El,Ho); termine diagonale e Lacune fuori diagonale   
%eval(FOR)
for indT=1:length(Tvet) % temperature loop
    [GT,EsT,DepT,RspT,eD,hD,EFc,EFv]=pf_functionFit(eDensityv,hDensityv,lambdavet,Tvet,indT,Ban,mesh,mode,Pargain);
    
    GE(:,:,:,indT)=GT/NQW;
    DepE(:,:,:,indT)=DepT/NQW;
    EsE(:,:,:,indT)=EsT/NQW;
    RspE(:,:,indT)=RspT/NQW;
    
    eD1E(:,indT)=eD/NQW;
    hD1E(:,indT)=hD/NQW;
    EFcvE(:,indT)=EFc;
    EFvvE(:,indT)=EFv;    
    
    disp(['Loop 2, temperature ',num2str(indT),' of ',num2str(length(Tvet))])

end% Temp
Loop_2_Time=toc

GH=zeros(lP,lP,lL,lT);
RspH=zeros(lP,lP,lT);
EsH=GH;
DepH=GH;  

eDensityv=NQW*Densityv;
hDensityv=NQW*Densityv*mode.DeltaDensity_Perc;
   
tic
% Calcolo G(El,El)  e G(El,Ho); termine diagonale e Lacune fuori diagonale   
%eval(FOR)

for indT=1:length(Tvet) % temperature loop

    [GT,EsT,DepT,RspT,eD,hD,EFc,EFv]=pf_functionFit(eDensityv,hDensityv,lambdavet,Tvet,indT,Ban,mesh,mode,Pargain);
    
    GH(:,:,:,indT)=GT/NQW;
    DepH(:,:,:,indT)=DepT/NQW;
    EsH(:,:,:,indT)=EsT/NQW;
    RspH(:,:,indT)=RspT/NQW;

    eD1H(:,indT)=eD/NQW;
    hD1H(:,indT)=hD/NQW;
    EFcvH(:,indT)=EFc;
    EFvvH(:,indT)=EFv;    
    
    disp(['Loop 3, temperature ',num2str(indT),' of ',num2str(length(Tvet))])

end% Temp
Loop_3_Time=toc

DeltaN = Dep/2*Nb; % questo ? il DeltaN assoluto
DeltaNE = DepE/2*Nb; % questo ? il DeltaN assoluto
DeltaNH = DepH/2*Nb; % questo ? il DeltaN assoluto
lav=lambdavet*1e6;
Tv=Tvet;
port=Densityv/(mesh.Lz*100)*1e-18;

Gtot=G*c_light/Nb*1e-9*100;
GtotE=GE*c_light/Nb*1e-9*100;
GtotH=GH*c_light/Nb*1e-9*100;

meshQW=mesh;
modeQW=mode;

port_2De=diag(eD1);
port_2Dh=diag(hD1);

%'port', keyboard

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
 figure, plot(Ban.kgrid,Ban.SBV,'linewidth',1.5), hold on    
 plot(Ban.kgrid,SBVap,'.','linewidth',2), %ylim([0 .2])
 pausak

port_2D=Densityv;

DeltaN_Perc=mode.DeltaDensity_Perc;

var=[' port_2D port_2De port_2Dh eD1 hD1 lav Tv G Es Rsp DeltaN '];
varder='GH GE EsE EsH RspE RspH DeltaNE DeltaNH DeltaN_Perc ';
var_more=['MefH Ban meshQW modeQW'];
eval(['save ',mode.fileNameLUT,'_more.mat ',var_more])

eval(['save ',mode.fileNameLUT,'.mat ',var, varder])




%var=[' port_2D port_2De port_2Dh eDm hDm eD1 hD1 lav Tv G Es Rsp DeltaN '];
%varder='GE EsE EsH RspE RspH DeltaNE DelatNH DeltaN_Perc ';

for indTemp=1:length(Tv)

por=eD1(:,indTemp)'*1e-4;
porInc=eD1E(:,indTemp)'*1e-4;
porE=por;
porH=hD1(:,indTemp)'*1e-4;
porHInc=hD1H(:,indTemp)'*1e-4;
iDele=diag(1./(porInc-por));
iDlac=diag(1./(porHInc-porH));
     Rspt=squeeze(Rsp(:,:,indTemp));
     RsptdE=squeeze(RspE(:,:,indTemp));
     RsptdH=squeeze(RspH(:,:,indTemp));
     dRE=(RsptdE-Rspt)*iDele;
     dRH=iDlac*(RsptdH-Rspt);
     dRicdE(:,:,indTemp)=dRE;  
     dRicdH(:,:,indTemp)=dRH;     
%     'rixc', keyboard
    for indlambda=1:length(lav)
        % indlambda=1; indTemp=1; % debug purposes
        
     Gtemp=squeeze(G(:,:,indlambda,indTemp));
     GtempdE=squeeze(GE(:,:,indlambda,indTemp));
     GtempdH=squeeze(GH(:,:,indlambda,indTemp));
     dGE=(GtempdE-Gtemp)*iDele;
     dGH=iDlac*(GtempdH-Gtemp);
     dGdE(:,:,indlambda,indTemp)=dGE;  
     dGdH(:,:,indlambda,indTemp)=dGH;  
     
     Etemp=squeeze(Es(:,:,indlambda,indTemp));
     EtempdE=squeeze(EsE(:,:,indlambda,indTemp)); 
     EtempdH=squeeze(EsH(:,:,indlambda,indTemp));    
     dEE=(EtempdE-Etemp)*iDele;
     dEH=iDlac*(EtempdH-Etemp);
     dEdE(:,:,indlambda,indTemp)=dEE;  
     dEdH(:,:,indlambda,indTemp)=dEH;       
    end     
end

%'qui', keyboard

porE_Rsp=[0,porE];
porH_Rsp=[0,porH];

Rsp0=zeros(size(Rsp(:,:,1))+1);
for kT=1:length(Tv)
 Rdu=Rsp(:,:,kT);
 Rsp0(2:end,2:end)=Rdu;
 Ric(:,:,kT)=Rsp0;
 Rdu=dRicdE(:,:,kT);
 Rsp0(2:end,2:end)=Rdu;
 dRspdEF(:,:,kT)=Rsp0; 
 Rdu=dRicdH(:,:,kT);
 Rsp0(2:end,2:end)=Rdu;
 dRspdHF(:,:,kT)=Rsp0;  
end

dRicdE=dRspdEF;
dRicdH=dRspdHF;
por_E=porE;
por_H=porH;

fileSav=[mode.fileNameLUT,'_Der.mat '];

var=' G dGdE dGdH lav por_E por_H porE_Rsp porH_Rsp Ric dRicdE dRicdH Es dEdE dEdH Tv DeltaN';

%'dentro process', keyboard
eval([' save ',fileSav,var])

