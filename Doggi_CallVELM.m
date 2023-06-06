close all
clear
colordef black
% dbstop if error

addpathVENUS



% load LW_Stephan_GRok_.mat
% % load LW_Julian_TJ_2AR_pDBR_fake_.mat
% load LW_Julian_TJ_2AR_pDBR_fake_new.mat

% load LW_MarkusN_BTJetch_DD_ROSSO_LgDBR_.mat

% load('out\LW_MarkusN_FINALE_LgDBR_oldFit_FCAindox_FatPerd=2.4_Tindox_Texp=1.2_noDoubleTFCA_eMobQW=0.2_NX=2.5.mat');

% load('out\LW_MarkusN_TJ_oxAbove_LgDBR_fixed_oldFit_FCAindox_FatPerd=2.4_Tindox_Texp=1.2_noDoubleTFCA_eMobQW=0.2_NX=2.5.mat');
% 
% load('out\LW_MarkusN_TJ_oxBelow_LgDBR_fixed_oldFit_FCAindox_FatPerd=2.4_Tindox_Texp=1.2_noDoubleTFCA_eMobQW=0.2_NX=2.5.mat');

% load LW_MarkusN_FINALE_LgDBR_FatPerd0=2.8_Texp=0_eQWmob=.4_dndT=2.3_betaT=1.mat

% load LW_MarkusN_FINALE_LgDBR_110_FatPerd0=2.9_Texp=0_eQWmob=.4_dndT=2.37_betaT=1.1_tauSRH=10ns_fatRad=0.05_Pmin=0.1_fCondTer=0.9_FatPerCoefTemp=0.002_FattoreZ=1.0015.mat

load('C:\Users\albig\Politecnico\Dottorato\3b_VENUS\VENUS_Unified\out\LW_MarkusN_TJ_oxBelow_LgDBR_20C_Pelec-PTJ_noHeatTJ.mat')
% iold=input('  old VELM? [ENTER]: Yes; [Any key]: No  ');
% iold=[];
% if isempty(iold)
%     rmpath('OtticoBar\new23OpticaGR')
%     rmpath('OtticoBar\new22Optica')
% %     load LW_MarkusN_FINALEdisFitto_testOLDvelm_abs.mat
% %     load LW_MarkusN_FINALE_LgDBR_oldVELM.mat
%     %     load LW_MarkusN_FINALEdisFitto_oldVELM.mat
% else
%     rmpath('OtticoBar\new22Optica')
% %     load LW_MarkusN_FINALE_LgDBR_testNEWvelm.mat
%     %     load LW_MarkusN_FINALEdisFitto_newVELM.mat
% end

ParVet=MODEplot{1}.ParVet;
VelmOptions=MODEplot{1}.VelmOptions;

structureName=[nomeSR,mode.strName];
fis= strfind(structureName,'\');
strName=structureName(fis(end)+1:end);

% fil_str=[structureName,'.str'];
fil_str=['dati\',strName,'.str'];
    
verVE=0;    % verbVELM
ian=1;      % guiding
MulT=1;     % Temperature pre-factor

nqw=3;      % Number of QW

% mode.DT0=200;
% NP_k=20;

% load velmSa

if nqw==1
    fat_gain=2.6;
else
    fat_gain=1;
end

s_LoadConstants
vph=Clight./mode.nindexQW;


vind=[];    % DD indexes of VELM called (indVoltage)
veind=[];   % Progressive vector of VELM calling (indVELM)
for k=1:length(VELMInput)-1
    vv=VELMInput(k).indVoltage;
    vev=VELMInput(k).indVELM;
    if length(vv)==1
        vind=[vind vv];
        veind=[veind vev];
        fPold(k,:)=VELMInfo(k).fPdif;
        Laold(k,:)=VELMInfo(k).vlambda;
        Ga(k,:)=VELMInfo(k).Lm/vph;
        Tvet(k)=VELMInfo(k).DeltaTmax;
        E2v(:,:,k)=VELMInfo(k).E2';
    end
end
modeold=mode;

%'corrente', keyboard
Pcor=input(' Corrente = ')
% Pcor=.6

Cor=1000*modeold.ii_dd(vind);
% Cor=1000*modeold.ii_dd;         % Computed current from DD
[~,fim]=min(abs(Pcor-Cor));     % find closest current index to Cor (DD)
pu=fim;                         % for plot purposes

% ves= find(vind==fim);           % index of VELM calling (1st, 2nd,...)
ves=fim;           % index of VELM calling (1st, 2nd,...)
vindm=vind(1:end-1);

IPLOT=1

if IPLOT==1
    figure(1), hold on, grid on
    plot(1000*mode.ii_dd,modeold.Gmod,'g','LineWidth',2)
    plot(1000*mode.ii_dd,modeold.Lmod,'r--','LineWidth',2)
    if(mode.nmodes>1)
        plot(1000*modeold.ii_dd(vind),modeold.Gmod(:,vind)','go','LineWidth',2)
        plot(1000*modeold.ii_dd(vind),modeold.Lmod(:,vind)','ro','LineWidth',2)
    else
        plot(1000*modeold.ii_dd(vindm),modeold.Gmod(vindm),'go','LineWidth',2)
        plot(1000*modeold.ii_dd(vind),modeold.Lmod(vind),'ro','LineWidth',2)
    end
    xlabel('Current, mA')
    ylabel('Modal gain vs losses, cm^{-1}')
    ylim([0 100])
    pausak
    
    I=1000*modeold.ii_dd(vind);     % Current when VELM is called
    I(ves)
    I0=Pcor;                        % target current
    
    figure(2), hold on, grid on
    plot(I,Tvet(veind),'LineWidth',2)   % VELM current vs VELM Temperature
    plot(I(ves),Tvet(ves),'mo','LineWidth',2),
    xlabel('Corrente')
    ylabel(' Tmax ')
    pausak
    
    figure(3), hold on, grid on
    plot(1e4*mesh.xgrid,squeeze(E2v(:,:,ves)),'LineWidth',1.5),
    xlabel(' \rho, \mum ')
    ylabel(' Field intensity ')
    xlim([0 10])
    pausak
    
    
    figure(4), hold on, grid on
    plot(I,squeeze(E2v(10,:,veind)),'LineWidth',2),
    if size(E2v,2)>1
        hold on, plot(I(ves),squeeze(E2v(10,:,ves)),'mo','LineWidth',2),
    end
    xlabel('Corrente, mA')
    ylabel(['Fields intensity at x(10)=',num2str(mesh.xgrid(10)*1e4),' \mum'])
    pausak
    
    figure(5), hold on, grid on
    plot(I,fPold(veind,:),'LineWidth',2),
    plot(I(ves),fPold(ves,:),'mo','LineWidth',2),
    xlabel('Corrente')
    ylabel(' fPdif ')
    pausak
    
    figure(6), hold on, grid on
    plot(I,Laold(veind,:),'LineWidth',2),
    plot(I(ves),Laold(ves,:),'mo','LineWidth',2),
    xlabel('Corrente')
    ylabel(' Lambda ')
    pausak
end
%return

%'corrente', keyboard
%Pcor=input(' Corrente = ')

Cor=1000*modeold.ii_dd(vind);
[~,fim]=min(abs(Pcor-Cor));
veind=fim;

ico=0;
for indVELM=veind
    ico=ico+1;
    
    mode.verbVELM=verVE;
    %     VelmOptions.ianti_gui=ian;  % 0 per LP
    DL=VelmOptions.Dlam;  % 0 per LP
    %DL(5)=.8;
    %VelmOptions.Dlam=DL;  % 0 per LP
    %     VelmOptions.gain_gui=ian;  % 0 per LP
    %     VelmOptions.NP_k=NP_k;  % 0 per LP
%     VelmOptions.itutmir=0; % 1: caso termico completo; 0: caso ridotto (+ veloce)
%         VelmOptions.Pf.nmasce=-3;
    
    mode.matgain=VELMInput(indVELM).matgain;
    mode.DeltaN=VELMInput(indVELM).DeltaN;
    mode.Deltan=VELMInput(indVELM).Deltan;
    mode.efield_y=VELMInput(indVELM).efield_y;
    mode.efield_x=VELMInput(indVELM).efield_x;
    mode.vlambda=VELMInput(indVELM).vlambda;
    mode.alpha=VELMInput(indVELM).alpha;
    
    mode.elecABSvelm=VELMInput(indVELM).elecABS;
    mode.holeABSvelm=VELMInput(indVELM).holeABS;
    if indVELM>1
        mode1.TmVelm=VELMInput(indVELM).TmVelm;
        mode1.LamVelm=VELMInput(indVELM).LamVelm;
    else
        mode1.a=0;
    end
    mesh.DeltaTvelm=MulT*VELMInput(indVELM).DeltaTvelm;
    mesh.ygrid=VELMInput(indVELM).ygrid;
    mesh.xgrid=VELMInput(indVELM).xgrid;
    mesh.nnx=VELMInput(indVELM).nnx;
    mesh.nny=VELMInput(indVELM).nny;
    
%     VelmOptions.Dlam=[-5 +5 10 0 .4]; % ORIGINAL, -10???? (30°C)
%     VelmOptions.Dlam=[-.5 +7 10 0 .4]; % ORIGINAL, -10???? (30°C)
% %     VelmOptions.Dlam=[-6 -2 10 0 .4]; % ORIGINAL, -10???? (room T)
%     VelmOptions.Dlam=[-.5 4.5 10 0 .4];
%     VelmOptions.NP_k=30;
%     VelmOptions.krel_max=.12;
%     VelmOptions.ianti_gui=0;
%     
%     ParVet(3)=8;
%     mode.verbVELM=2;
    
    'ver prima di call', keyboard
    
    [velm] = f_CallVELM(mesh,mode,mode1,ParVet,VelmOptions,fil_str);
%     [velm] = f_CallVELM_old(mesh,mode,mode1,ParVet,VelmOptions,fil_str);
    indv=vind(ves);
    
    Lmod=velm.Lm/vph;
    lambda=velm.vlambda;
    fP=velm.fPdif;
    E2=velm.E2'; %dà errore qui (multimodo)
end

colordef white

indQW=2;
inQW = mesh.inMQW{indQW};
% xQW = mesh.node(1,inQW);
xQW=mesh.xgrid(1:mesh.nnxQW{1})*1e4;
% figure, plot(xQW,squeeze(E2v(1:mesh.nnxQW{1},:,pu)),'LineWidth',2),
figure, plot(xQW,E2(1:mesh.nnxQW{1},:),'LineWidth',2),
xlabel(' \rho, \mum ')
ylabel('Optical Field intensity, a.u.')
pausak

if IPLOT==1
    vindm=vind(1:end-1);
    figure(1)
    plot(mode.ii_dd(indv)*1e3,Lmod,'y+','LineWidth',2)
    pausak
    
    figure(3)
    plot(1e4*mesh.xgrid,E2,'--','LineWidth',2)
    xlim([0 5])
    pausak
        
    figure(4)
    plot(I(ves),squeeze(E2(10,:)),'g*')
    ylabel(['Fields intensity at x(10)=',num2str(mesh.xgrid(10)*1e4),' \mum'])
    pausak
    
    figure(5)
    plot(I(ves),fP,'g*')
    pausak
    
    figure(6)
    plot(I(ves),lambda,'g*')
    pausak
    %
end
