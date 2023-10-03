close all
% clear
clear global
colordef black
% dbstop if error

addpathVENUS

% load LW_MarkusN_BTJetch_DD_ROSSO_LgDBR_.mat

% load LW_MarkusN_TJ_oxAbove_LgDBR_FINAL.mat
% load LW_MarkusN_TJ_oxBelow_LgDBR_NUSOD2023.mat

% load out\LW_MarkusN_TJ_oxAbove_LgDBR_20C_lambda_radial.mat
% load out\LW_MarkusN_TJ_oxBelow_LgDBR_20C_lambda_radial.mat

% load out\LW_MarkusN_FINALE_LgDBR_new.mat % 1D
% load out\LW_MarkusN_FINALE_LgDBR_1D_PerCoefExT=2e-2.mat % 1D
load out\LW_MarkusN_FINALE_LgDBR_20C_radial.mat % 3D

ParVet=MODEplot{1}.ParVet;
VelmOptions=MODEplot{1}.VelmOptions;

structureName=[nomeSR,mode.strName];
fis= strfind(structureName,'\');
strName=structureName(fis(end)+1:end);

% fil_str=[structureName,'.str'];
fil_str=['dati\',strName,'.str'];
    
verVE=2   % verbVELM

IPLOT=1

% if mode.quasi1D==1
%     ian=0;      % NO anti-guiding
% else
%     ian=1;      % anti-guiding
% end
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
if length(Pcor)==0
 Pcor=8
end 

Cor=1000*modeold.ii_dd(vind);
% Cor=1000*modeold.ii_dd;         % Computed current from DD
[~,fim]=min(abs(Pcor-Cor));     % find closest current index to Cor (DD)
pu=fim;                         % for plot purposes

% ves= find(vind==fim);           % index of VELM calling (1st, 2nd,...)
ves=fim;           % index of VELM calling (1st, 2nd,...)
vindm=vind(1:end-1);


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
    
    if mode.quasi1D==0
    figure(4), hold on, grid on
    plot(I,squeeze(E2v(10,:,veind)),'LineWidth',2),
    if size(E2v,2)>1
        hold on, plot(I(ves),squeeze(E2v(10,:,ves)),'mo','LineWidth',2),
    end
    xlabel('Corrente, mA')
    ylabel(['Fields intensity at x(10)=',num2str(mesh.xgrid(10)*1e4),' \mum'])
    pausak
    end
    
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
    if exist('ian')
        VelmOptions.ianti_gui=ian;  % 0 per LP
    end
%     DL=VelmOptions.Dlam;  % 0 per LP
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
    
%     load fca
%     mode.elecABSvelm=fca.eleccentro';
%     mode.holeABSvelm=fca.holecentro';    
    mode.elecABSvelm=VELMInput(indVELM).elecABS;
    mode.holeABSvelm=VELMInput(indVELM).holeABS;

    if indVELM>1
        mode1.TmVelm=VELMInput(indVELM).TmVelm;
        mode1.LamVelm=VELMInput(indVELM).LamVelm;
    else
        mode1.a=0;
    end
    mesh.DeltaTvelm=MulT*VELMInput(indVELM).DeltaTvelm;
	iTaroccoT=0
    if iTaroccoT==1
        IPLOT=1;
        StrTT=MODEplot{1}.StrTT;
        'TERMICO TAROCCO', keyboard
        [deltaT,PTherm,T_Contributi,fattore_correttivo,condzTe] = f_ThermD1ANA(mesh,mode,StrTT,IPLOT);
        R= reshape(deltaT,1,prod(size(mesh.DeltaTvelm)));
        'Fine TERMICO TAROCCO', keyboard
        mesh.DeltaTvelm=R*1.4;
    end
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
% %     ParVet(3)=8;

%     mode.verbVELM=0;
    mode.verbVELM=verVE;
    
    'ver prima di call', keyboard
    
    if mode.quasi1D==0
        [velm] = f_CallVELM(mesh,mode,mode1,ParVet,VelmOptions,fil_str);
%         [velm] = f_CallVELM_old(mesh,mode,mode1,ParVet,VelmOptions,fil_str);
    else
        [velm] = f_CallVELM_1D(mesh,mode,mode1,ParVet,VelmOptions,fil_str,mode.quasi1D);
    end

    indv=vind(ves);
    
    Lmod=velm.Lm/vph;
    lambda=velm.vlambda;
    fP=velm.fPdif;
    if isfield(mode,'quasi1D') && mode.quasi1D==1
        velm.E2=ones(mode.nmodes,mesh.nnx)/(mode.AreaOx);
        velm.E2=ones(mode.nmodes,mesh.nnx)/(mode.AreaOpt);
        velm.SW=0;
    end
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
        
    if mode.quasi1D==0
        figure(4)
        plot(I(ves),squeeze(E2(10,:)),'g*')
        ylabel(['Fields intensity at x(10)=',num2str(mesh.xgrid(10)*1e4),' \mum'])
        pausak
    end
    
    figure(5)
    plot(I(ves),fP,'g*')
    pausak
    
    figure(6)
    plot(I(ves),lambda,'g*')
    pausak
    %
end
