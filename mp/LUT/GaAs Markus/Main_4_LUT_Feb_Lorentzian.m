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
mesh.Lz=7.9e-9; % quantum well width, m
mesh.xmol_barrier=0.270; % Al molar fraction (barrier)
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
mode.fileNameLUT='LUT4D_Feb_Lorentzian_noRen'; % filename for LUT
mode.LUT=1; % if 0 enables additional saves/plots
mode.DeltaDensity_Perc=1+0.005; % increment of carrier density for Jacobian derivatives
mode.iplot=0; % if 1, several intermediate plots are produced
mode.flagLimitMaxDensity=0; % 1: limits max carrier density; 0: doesn't
mode.iline='lorentzian'; % 'nMark' ; 'lorentzian'; 'Landsberg'; Landsberg_extended' and Landsberg_modified  models 
mode.Expgammak=6; % exponent for carrier-carrier scattering corrections
mode.iren=0; % enable gap renormalization computation
mode.ieh_equal=0; % 1: e- and h- densities are assumed equal; 0: they aren't
mode.ifit=1; % 1: use a larger k grid by fitting it with spline; 0: don't
mode.vic=[]; % if empty, all transitions are computed
mode.viv=[]; % if empty, all transitions are computed
% % to compute the contributions to gain and spontaneous emission with
% % few transitions, use the following lines
% % mode.vic=1;
% % mode.viv=1:2;
%==========================================================================
% Optical response: parameters
%==========================================================================
% Tvet=[300 350];
% lambdavet=linspace(800,870,101)*1e-9;
% Densityv=[ 0.01 5 10]*1e12;
Tvet=[290:10:520];
Densityv=[0.0001 0.001 0.01 0.1 1 2 2.5:.5:8 9:20]*1e12;
lambdavet=linspace(830,870,41)*1e-9;
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
tic
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
    
    disp(['Loop 1, temperature ',num2str(indT),' of ',num2str(length(Tvet))])
    
end% Temp
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

eDensityv=Densityv*mode.DeltaDensity_Perc;
hDensityv=Densityv;
   
tic
% Calcolo G(El,El)  e G(El,Ho); termine diagonale e Lacune fuori diagonale   
parfor indT=1:length(Tvet) % temperature loop
    [GT,EsT,DepT,RspT,eD,hD,EFc,EFv]=pf_functionFit(eDensityv,hDensityv,lambdavet,Tvet,indT,Ban,mesh,mode,Pargain);
    
    GE(:,:,:,indT)=GT;
    DepE(:,:,:,indT)=DepT;
    EsE(:,:,:,indT)=EsT;
    RspE(:,:,indT)=RspT;
    
    disp(['Loop 2, temperature ',num2str(indT),' of ',num2str(length(Tvet))])

end% Temp
Loop_2_Time=toc

GH=zeros(lP,lP,lL,lT);
RspH=zeros(lP,lP,lT);
EsH=GH;
DepH=GH;  

eDensityv=Densityv;
hDensityv=Densityv*mode.DeltaDensity_Perc;
   
tic
% Calcolo G(El,El)  e G(El,Ho); termine diagonale e Lacune fuori diagonale   
parfor indT=1:length(Tvet) % temperature loop

    [GT,EsT,DepT,RspT,eD,hD,EFc,EFv]=pf_functionFit(eDensityv,hDensityv,lambdavet,Tvet,indT,Ban,mesh,mode,Pargain);
    
    GH(:,:,:,indT)=GT;
    DepH(:,:,:,indT)=DepT;
    EsH(:,:,:,indT)=EsT;
    RspH(:,:,indT)=RspT;
    
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
 plot(Ban.kgrid,SBVap,'.','linewidth',2), ylim([0 .2])
 pausak

port_2D=Densityv;

DeltaN_Perc=mode.DeltaDensity_Perc;

var=[' port_2D port_2De port_2Dh eD1 hD1 lav Tv G Es Rsp DeltaN '];
varder='GH GE EsE EsH RspE RspH DeltaNE DeltaNH DeltaN_Perc ';
var_more=['MefH Ban meshQW modeQW'];
eval(['save ',mode.fileNameLUT,'_more.mat ',var_more])

eval(['save ',mode.fileNameLUT,'.mat ',var, varder])


nsum=1e6*[1:3:300];
por=eD1(:,1)'*1e-4;
di=find(diff(por)==0)+1;
por(di)=por(di)+nsum(di);
porE=por;
porH=hD1(:,1)'*1e-4;
iDele=diag(1./((mode.DeltaDensity_Perc-1)*porE));
iDlac=diag(1./((mode.DeltaDensity_Perc-1)*porH));

%var=[' port_2D port_2De port_2Dh eDm hDm eD1 hD1 lav Tv G Es Rsp DeltaN '];
%varder='GE EsE EsH RspE RspH DeltaNE DelatNH DeltaN_Perc ';

for indTemp=1:length(Tv)
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

