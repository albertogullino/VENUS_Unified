% Look at presentation: "20230929_HeatSources_Stephan.pptx"

% close all
% clear
% clear global
% colordef black
% dbstop if error

addpathVENUS


% load LW_MarkusN_TJ_oxAbove_LgDBR_FINAL.mat
% load LW_MarkusN_TJ_oxBelow_LgDBR_NUSOD2023.mat

% load out\LW_MarkusN_TJ_oxAbove_LgDBR_20C_lambda_radial.mat
% load out\LW_MarkusN_TJ_oxBelow_LgDBR_20C_lambda_radial.mat

% load out\LW_MarkusN_FINALE_LgDBR_new.mat % 1D
% load out\LW_MarkusN_FINALE_LgDBR_1D_PerCoefExT=2e-2.mat % 1D
% load out\LW_MarkusN_FINALE_LgDBR_20C_radial.mat % 3D

% load out\LW_MarkusN_FINALE_LgDBR_20C_Stephan.mat
% load out\LW_MarkusN_FINALE_LgDBR_20C_Stephan_repelem.mat
% load out\LW_MarkusN_FINALE_LgDBR_xB=0.20_DeltaEgFIX_newFullLUT.mat
load out\LW_MarkusN_FINALE_LgDBR_xB=0.38_DeltaEgFIX_newFullLUT.mat

ParVet=MODEplot{1}.ParVet;
VelmOptions=MODEplot{1}.VelmOptions;
StrTT=MODEplot{1}.StrTT;

structureName=[nomeSR,mode.strName];
fis= strfind(structureName,'\');
strName=structureName(fis(end)+1:end);

% fil_str=[structureName,'.str'];
fil_str=['dati\',strName,'.str'];

IPLOT=1

s_LoadConstants

x=mesh.xgrid*1e4;
y=mesh.ygrid*1e4;

[~,irox]=min(abs(x-StrTT.Rox));

modet=mode;
%
Pcor=input(' Corrente = ')
if length(Pcor)==0
    Pcor=8
end

Cor=1e3*mode.ii_dd;
% Cor=1000*modeold.ii_dd;         % Computed current from DD
[~,iCurr]=min(abs(Pcor-Cor));     % find closest current index to Cor (DD)

%% Extract the quantities needed to call again the Thermal Solver
%
% Absorption coefficient (alpha)
ifFC=mode.ifFC;
if ifFC==1
    elecABS=squeeze(MODEplot{1}.elec(iCurr-1,:,:));
    holeABS=squeeze(MODEplot{1}.hole(iCurr-1,:,:));
    
    if mode.flgBTJ==1
        iRagTJ=mesh.nnxQW{1};
        for iTJ=1:iRagTJ
            elecABS(mesh.LBTJ(iTJ):mesh.RBTJ(iTJ))=elecABS(mesh.LBTJ(iTJ):mesh.RBTJ(iTJ))+mode.dop_dp(mesh.LBTJ(iTJ):mesh.RBTJ(iTJ));
            holeABS(mesh.LBTJ(iTJ):mesh.RBTJ(iTJ))=holeABS(mesh.LBTJ(iTJ):mesh.RBTJ(iTJ))+mode.dop_am(mesh.LBTJ(iTJ):mesh.RBTJ(iTJ));
        end
    end
elseif ifFC==0
    elecABS=mode.dop_dp;
    holeABS=mode.dop_am;
elseif ifFC==2
    elecABS=(mode.dop_dp+mode.elec)/2;
    holeABS=(mode.dop_am+mode.hole)/2;
end
elecABS=elecABS*1e-18*mode.CarrierNorm;
holeABS=holeABS*1e-18*mode.CarrierNorm;

T300=mode.T300;
Temp=(squeeze(MODEplot{1}.Temp(iCurr-1,:,:)))-T300;

ABS_Texp=mode.ABS_Texp+mode.PerCoefExT*Temp;

modet.elecABS=elecABS.*(1+Temp/T300).^ABS_Texp;
modet.holeABS=holeABS.*(1+Temp/T300).^ABS_Texp;

if isfield(mode,'Fat_PerCoefTemp')
    Fat_Perd=mode.Fat_Perd0+Temp*mode.Fat_PerCoefTemp;
    ABSh=mode.ABSh0.*Fat_Perd;
    ABSe=mode.ABSe0.*Fat_Perd;
else
    ABSh=mode.ABSh;
    ABSe=mode.ABSe;
end

modet.alpha=(modet.holeABS.*ABSh+modet.elecABS.*ABSe)*100;
modet.alpha(not(mesh.iq))=0;
%
% Heat sources (taken at the bias point iCurr-1, due to the lack of self
% consistency between electro-optical solver and thermal solver)
modet.HeatJoule=squeeze(MODEplot{1}.HeatJoule(iCurr-1,:,:));        % Joule effect
modet.HeatRec_Cap=squeeze(MODEplot{1}.HeatRec_Cap(iCurr-1,:,:));    % QW capture
modet.HeatRec_RAD=squeeze(MODEplot{1}.HeatRec_RAD(iCurr-1,:,:));    % Radiative GR
modet.HeatRec_13=squeeze(MODEplot{1}.HeatRec_13(iCurr-1,:,:));      % non-Rad GR
modet.HeatOptAbs=squeeze(MODEplot{1}.HeatOptAbs(iCurr-1,:,:));      % FCA losses
% 
modet.PDissPred(end)=mode.PDissPred(iCurr);
%
% Temperature at the previous step, required by Thermic solver
Tprec=squeeze(MODEplot{1}.Tprec(iCurr-1,:,:));
%
% DeltaT from VENUS simulation, at current index "iCurr" (for plot purpose)
T=squeeze(MODEplot{1}.Temp(iCurr,:,:))-T300;

colordef white

% Plot results from VENUS simulation
if IPLOT==1
    %% DeltaT
    figure
    hold on
    grid on
    box on
    plot(mode.ii_dd*1e3,mode.DeltaTmax,'k','LineWidth',2)
    plot(mode.ii_dd*1e3,mode.DeltaTmax_Joule,'c--','LineWidth',2)
    plot(mode.ii_dd*1e3,mode.DeltaTmax_OptAbs,'g--','LineWidth',2)
    plot(mode.ii_dd*1e3,mode.DeltaTmax_Ccap,'r--','LineWidth',2)
    plot(mode.ii_dd*1e3,mode.DeltaTmax_RAD,'b--','LineWidth',2)
    plot(mode.ii_dd*1e3,mode.DeltaTmax_srhAu,'y--','LineWidth',2)
    xline(Pcor,'--')
    legend('\DeltaT','Joule','Opt. abs.','Ccap','Rad','SRH/Aug','location','northwest')
    xlabel('Current, mA')
    ylabel('Temperature rise, K')
    set(gca,'FontSize',16,'FontName','Times','Box','on')

    % Compute the thermal power from each source, as a function of current 
%     PthermIV
    
    figure
    set(gcf,'pos',[631         454        1276         420])
    subplot(131)
    surf(x,y,T)
    shading interp
    view(2)
    ylim([0 120])
    colorbar
    xlabel('\rho, \mum'),ylabel('z, \mum')
    set(gca,'FontSize',16,'FontName','Times new roman')
    pausak
    axis([0 8 114 118.5])
    xline(0,'r--','linewidth',4)
    xline(2,'g--','linewidth',4)
    xline(5,'c--','linewidth',4)
    yline(115,'k--','linewidth',4)
    yline(118.5,'m--','linewidth',4)
    
    subplot(132)
    hold on,grid on
    xlabel('z, \mum')
    ylabel('\DeltaT, K')
    plot(y,T(:,1),'r','linewidth',2)
    % set(gca,'yscale','log')
    set(gca,'FontSize',16,'FontName','Times new roman')
    plot(y,T(:,14),'g','linewidth',2)
    plot(y,T(:,25),'c','linewidth',2)
    legend('\rho=0','\rho=2 \mum','\rho=5 \mum','location','northwest')
    
    subplot(133)
    xlabel('z, \mum')
%     ylabel('\DeltaT, K')
    set(gca,'FontSize',16,'FontName','Times new roman')
    hold on,grid on
    plot(x,T(118,:),'k','linewidth',2)
    plot(x,T(end,:),'m','linewidth',2)
    legend('z=115 \mum (QW)','z=118.5 \mum (top)')
    legend('\rho=0','\rho=2 \mum','\rho=5 \mum')

    %% HeatSources
    Joule=modet.HeatJoule;
    FCA=modet.HeatOptAbs;
    GR=modet.HeatRec_13+modet.HeatRec_RAD;
    % var=modet.HeatRec_RAD;
    Ccap=modet.HeatRec_Cap;
    
    figure
    set(gcf,'pos',[631         454        1276         420])
    subplot(121)
%     surf(x,y,real(log10(Joule))),title('Joule (log)')
    % surf(x,y,real(log10(Ccap))),title('Ccap')
    % surf(x,y,real(log10(GR))),title('GR')
    % surf(x,y,real(log10(FCA))),title('FCA')
    surf(x,y,real(log10(Joule+FCA+Ccap))),title('FCA+Ccap+Joule (log)')
    caxis([-10 -3])
    shading interp
    view(2)
    colorbar
    xlabel('\rho, \mum')
    ylabel('z, \mum')
    set(gca,'FontSize',16,'FontName','Times new roman')
    axis([0 8 114 118.5])
    
    subplot(122)
%     surf(x,y,Joule),title('Joule (lin)')
    % surf(x,y,real(log10(Ccap))),title('Ccap')
    % surf(x,y,real(log10(GR))),title('GR')
    % surf(x,y,real(log10(FCA))),title('FCA')
    surf(x,y,Joule+FCA+Ccap),title('FCA+Ccap+Joule (lin)')
    shading interp
    view(2)
    colorbar
    xlabel('\rho, \mum')
    ylabel('z, \mum')
    set(gca,'FontSize',16,'FontName','Times new roman')
    axis([0 8 114 118.5])
    
    Prag=input(' Radius = ')
    if length(Prag)==0
        Prag=0
    end
    [~,iRag]=min(abs(Prag-x));
    
    figure
    set(gcf,'pos',[631         454        1276         420])
    subplot(121)
    hold on,grid on
    xlabel('z, \mum')
    ylabel('Heat sources, W/cm^3')
    plot(y,Joule(:,iRag),'linewidth',2)
    plot(y,FCA(:,iRag),'linewidth',2)
    plot(y,GR(:,iRag),'linewidth',2)
    plot(y,Ccap(:,iRag),'linewidth',2)
    plot(y,FCA(:,iRag)+Joule(:,iRag)+GR(:,iRag)+Ccap(:,iRag),'k--','linewidth',1.5)
    legend('Joule','FCA','GR','Ccap','location','southeast')
    axis([114 117 1e-15 1e0])
    set(gca,'yscale','log','FontSize',16,'FontName','Times new roman')
    
    subplot(122)
    hold on,grid on
    xlabel('z, \mum')
    plot(y,Joule(:,iRag),'linewidth',2)
    plot(y,FCA(:,iRag),'linewidth',2)
    plot(y,GR(:,iRag),'linewidth',2)
    plot(y,Ccap(:,iRag),'linewidth',2)
    plot(y,FCA(:,iRag)+Joule(:,iRag)+GR(:,iRag)+Ccap(:,iRag),'k--','linewidth',1.5)
    set(gca,'FontSize',16,'FontName','Times new roman')
    xlim([114 117])

end

    
%% Fake Heat source
% 
iFakeFlg=input(' Fake Heat source? [Enter or 1: YES; Any other key: NO]\n');
if length(iFakeFlg)==0
    iFakeFlg=1
end

if iFakeFlg==1
    % thermal sources are set to a value close to 0 (NO ZERO!)
    modet.HeatJoule=1e-30*modet.HeatJoule;          % Joule effect
    modet.HeatRec_Cap=1e-30*modet.HeatRec_Cap;      % QW capture
    modet.HeatRec_RAD=1e-30*modet.HeatRec_RAD;      % Radiative GR
    modet.HeatRec_13=1e-30*modet.HeatRec_13;        % non-Rad GR
    modet.HeatOptAbs=1e-30*modet.HeatOptAbs;        % FCA losses
    
    % Stephan way: Q=(Pin-Pout)/(Area*z)
    Pin=mode.ii_dd(iCurr-1)*mode.vv_dd(iCurr-1);    % W
    Pout=sum(mode.Pst_dd(:,iCurr-1),1)*1e-3;        % W
    
    % z: zbot -> bottom point; ztop -> top point

    % Values used in "20230929_HeatSources_Stephan.pptx"
%     zbot=110;          % um                                     
%     ztop=118.3;        % um
%     zbot=114.5;          % um                                     
%     ztop=116;            % um
%     zbot=114;              % um                                     
%     ztop=117.5;            % um
    %
%     % it was probably WRONG
%     zbot=112.5;              % um                                     
%     ztop=116.5;            % um
    % FIXED (values used in "20231115_HeatSources_MHerper.pptx")
    zbot=114;              % um                                     
%     ztop=116.5;            % um
    
    ztop=118.33;            % um (to verify the decrease of T after QWs)
%     ztop=y(end);            % um (to verify the decrease of T after QWs)

    [~,izbot]=min(abs(zbot-y)); % extract the corresponding mesh index
    [~,iztop]=min(abs(ztop-y)); % extract the corresponding mesh index
    dz=y(iztop)-y(izbot);       % Length along z

    % rho: r -> radial point for the area where the source is inserted
    % Values used in "20230929_HeatSources_Stephan.pptx"
%     rSimpl=StrTT.Rox+1;          % um
    rSimpl=StrTT.Rox;          % um
%     rSimpl=StrTT.ro_max;
    
    [~,iR]=min(abs(x-rSimpl));   % extract the corresponding mesh index
    
    % Compute the simplified source (Stephan way)
    Qfake=(Pin-Pout)/(pi*rSimpl^2*dz);   % W/um^3

    modet.HeatJoule(izbot:iztop,1:iR)=Qfake;    % inserted in Joule
    
    if IPLOT==1
        %% Plot the simplified source
        Joule=modet.HeatJoule;
        FCA=modet.HeatOptAbs;
        GR=modet.HeatRec_13+modet.HeatRec_RAD;
        % var=modet.HeatRec_RAD;
        Ccap=modet.HeatRec_Cap;

        figure
        set(gcf,'pos',[631         454        1276         420])
        subplot(121)
        surf(x,y,Joule+FCA+Ccap),title('FCA+Ccap+Joule (lin)')
        shading interp
        view(2)
        colorbar
        xlabel('\rho, \mum')
        ylabel('z, \mum')
        set(gca,'FontSize',16,'FontName','Times new roman')
        axis([0 8 zbot-1 y(end)])
        
        subplot(122)
        hold on,grid on
        xlabel('z, \mum')
        plot(y,Joule(:,iRag),'linewidth',2)
        plot(y,FCA(:,iRag),'linewidth',2)
        plot(y,GR(:,iRag),'linewidth',2)
        plot(y,Ccap(:,iRag),'linewidth',2)
        plot(y,FCA(:,iRag)+Joule(:,iRag)+GR(:,iRag)+Ccap(:,iRag),'k--','linewidth',1.5)
        set(gca,'FontSize',16,'FontName','Times new roman')
        xlim([zbot-1 y(end)])
    end
end

%% Call Thermal solver
fprintf('therm prima di call\n'), keyboard
% modet.iplotTerm=0;
% modet.iTfig=1;  % inner plot of thermal solver
modet.iTfig=0;  % inner plot of thermal solver
% 
% Thermal conductivities scaling
if iFakeFlg==1
    %% Herper data
    StrTT.Bcond.CondZmet=StrTT.Bcond.Cond_air;
    StrTT.Bcond.Condpas=StrTT.Bcond.Cond_air;
    %
    % % Uniform thermal conductivity
    % StrTT.Bcond.CondZmet=StrTT.Bcond.CondZc;
    % StrTT.Bcond.Cond_air=StrTT.Bcond.CondZc;
    % StrTT.Bcond.Condpas=StrTT.Bcond.CondZc;
    % StrTT.Bcond.CondZb=StrTT.Bcond.CondZc;

    % % POSSIBLY DANGEROUS!
    Tprec=0              % Avoid the first rescaling of k(T)
    modet.ioldTemp=-1;   % k(T) is not recomputed & fattore_correttivo=1!!!
else
    %% Reporoduce VENUS results from a stored workspace
    % mesh.fCondTer=1; % in VENUS: 0.90
    mesh.fCondTer=0.90; % in VENUS: 0.90
    mesh.fCondTerZ=mesh.fCondTer*0.72;
    % mesh.fCondTerZ=mesh.fCondTer;
    % mesh.fCondTerZ=1; % in VENUS, 0.72*mesh.fCondTer
    % mesh.fCondTerZ=0.72*mesh.fCondTer;
    %
    % 
end

if isfield(mode,'Oflg')==0 && mode.oflg==1
    modet.Oflg=1;
end

[DeltaT,Tprec_new,PTherm_new,T_Contributi,fattore_correttivo]=f_ThermicFunMOD(Tprec,mesh,modet,StrTT,IPLOT);

DeltaTJoule=T_Contributi{1};
DeltaTRec_srhAu=T_Contributi{2};
DeltaTRec_Ccap=T_Contributi{3};
DeltaTRec_RAD=T_Contributi{4};
DeltaTOptAbs=T_Contributi{5};

% IPLOT=0
% if IPLOT==1
    [DTtot,iTz,iTr]=f_maxMatrix(DeltaT,x,y,IPLOT);
    DTsum=sum([DeltaTJoule(iTz,iTr),DeltaTRec_srhAu(iTz,iTr),DeltaTRec_RAD(iTz,iTr),DeltaTOptAbs(iTz,iTr),DeltaTRec_Ccap(iTz,iTr)]);
    [DTjoule,iTzJou,~]=f_maxMatrix(DeltaTJoule,x,y,IPLOT);
    [DTfca,iTzFca,~]=f_maxMatrix(DeltaTOptAbs,x,y,IPLOT);
    [DTrec,iTzRec,~]=f_maxMatrix(DeltaTRec_srhAu+DeltaTRec_RAD,x,y,IPLOT);
    [DTccap,iTzCap,~]=f_maxMatrix(DeltaTRec_Ccap,x,y,IPLOT);
% end

% 
% DeltaTRec_srhAu(iTz,iTr)
% DeltaTRec_RAD(iTz,iTr)
% DeltaTOptAbs(iTz,iTr)
% DeltaTRec_Ccap(iTz,iTr)

%% Control plots
% T comes from VENUS, DeltaT is the recomputed quantity
if IPLOT==1
    figure(33)
    set(gcf,'pos',[516    89     1397     420])
    subplot(131)
    hold on,grid on,box on
    plot(y,DeltaT(:,1),'linewidth',2)
    plot(y,T(:,1),'--','linewidth',2)
    ylabel('\DeltaT, K'),xlabel('z, \mum')
    xlim([y(1) y(end)])
%     xlim([100 y(end)])
    legend('Doggi','VENUS','location','northwest')
    title('\rho=0')
    
    subplot(132)
    hold on,grid on,box on
    plot(y,DeltaT(:,irox),'linewidth',2)
    plot(y,T(:,irox),'--','linewidth',2)
    xlabel('z, \mum')
    xlim([y(1) y(end)])
%     xlim([100 y(end)])
    title('\rho=\rho_{ox}')

    subplot(133)
    hold on,grid on,box on
    plot(x,DeltaT(mesh.inMQW{2}(1),:),'linewidth',2)
    plot(x,T(mesh.inMQW{2}(1),:),'--','linewidth',2)
    xlabel('\rho, \mum')
    title('z=z_{QW}\approx 115 \mum')

    figure,hold on,box on,grid on
    plot(y,DeltaT(:,1),'k','linewidth',2)
    xlim([110 y(end)])
    plot(y,DeltaTJoule(:,1),'c','linewidth',2)
    plot(y,DeltaTOptAbs(:,1),'g','linewidth',2)
    plot(y,DeltaTRec_Ccap(:,1),'r','linewidth',2)
    plot(y,DeltaTRec_RAD(:,1),'b','linewidth',2)
    plot(y,DeltaTRec_srhAu(:,1),'y','linewidth',2)
%     xline(y(mesh.inMQW{2}(1)),'--')
    xlabel('z, \mum'),ylabel('\DeltaT, K')
    legend('Tot','Joule','FCA','Ccap','Rad','nr','location','northwest')
    title('\DeltaT profiles at \rho=0')
    set(gca,'FontSize',16,'FontName','Times','Box','on')
end
