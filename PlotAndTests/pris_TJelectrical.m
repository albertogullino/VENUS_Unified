IPLOT=0

% 0: no experimental results plotted;
% 1: FBH curves on GaAs and AlGaAs (x=0.25)
% 2: Chalmers, from 2020Tibaldi_PRAPP
iplotEXP=0

if IPLOT==1
    figure(70)
    hold on,grid on,box on
    plot(mesh.ygrid*1e4,mesh.xmol(1:mesh.nny),'-','LineWidth',2)
    xlabel('z, \mum'),ylabel('Al molar fraction')
    set(gca,'FontSize',14,'FontName','Arial')

    figure,
    hold on,grid on,box on
    plot(mesh.ygrid*1e4,mesh.dop_d(1:mesh.nny),'.-')
    plot(mesh.ygrid*1e4,mesh.dop_a(1:mesh.nny),'.-')
    xlabel('z, \mum'),ylabel('Doping levels, cm^{-3}')
    legend('N_D','N_A','location','best')
    set(gca,'yscale','log','FontSize',14,'FontName','Arial')
end

figure(79)
hold on,grid on,box on
plot(mode.vv_dd,squeeze(mode.VTJ(:,:,1)),'.-','linewidth',2)
xlabel('Voltage, V'),ylabel('V_{TJ}, V')
xlim([0 max(mode.vv_dd)])
set(gca,'FontSize',14,'FontName','Arial')

figure(80)
% subplot(121)
% hold on,grid on,box on
% plot(mode.vv_dd,mode.ii_dd/mode.AreaOx,'.-','linewidth',2)
% xlim([0 1])
% xlabel('Voltage, V'),ylabel('Current, A/cm^2')
% set(gca,'FontSize',14,'FontName','Arial')

% subplot(122)
hold on,grid on,box on
if iplotEXP==1
    load('FBHdata.mat')
    hold on
    plot(GaAs_Te(:,1),flip(GaAs_Te(:,2)),'r','linewidth',2)
    plot(AlGaAs_Te(:,1),flip(AlGaAs_Te(:,2)),'k','linewidth',2)
    chold
elseif iplotEXP==2
    %% Load NEGF-DD results from original folder (no Rs inserted)
    load IV_noR_2020Tibaldi_PRAPP.mat

    dox = 17; % DIFFERENT FROM 17um OF THE EXPERIMENTAL DATAAAAAAAA
    Area = pi*(dox/2*1e-4)^2; % cm^2
    Rs = input('Rs value?\n'); % ohm*cm^2, series resistance to fit the experimental curve

    Ibias = Jbias*Area;
    Vtot = abs(Vbias) + Rs*Ibias;

    %% Chalmers experimental curve, used in 2020Tibaldi_PRAPP: Rox=8.5 um
    load p_on_n_dox=17.mat

    Vmeas=Vvet;
    Imeas=Ivet;

    %% Plot section

    % I(V) plot
    plot(Vtot,1000*Ibias,'b','linewidth',2)
%     plot(Vtot,1000*Jbias,'b','linewidth',2)
    plot(Vvet,Imeas*1000,'bo','linewidth',2)
    axis([0 2.5 0.1 90])
    legend('Simulations','Measurements','Location','Best')
end

if exist("Rs")==0
    Rs=0
end

if mode.Zmat>0
    % v0_dd, BE CAREFUL when:
    % - Current driving
    % - Optical simulation ("Bisez"/"gambero")
    if iplotEXP==2
        plot(mode.v0_dd,mode.ii_dd*1e3,'r','LineWidth',2)
        ylabel('Current, mA')
        ylim([0.1 90])
        set(gca,'FontSize',14,'FontName','Arial')
    else
        plot(mode.v0_dd,mode.ii_dd/mode.AreaOx,'LineWidth',2)
        axis([0 max(mode.vv_dd+mode.ii_dd*Rs) 1e2 1e5])
        ylabel('Current density, A/cm^2')
        set(gca,'yscale','log','FontSize',14,'FontName','Arial')
    end
else%if Rs>0
    if iplotEXP==2
        % Plot current (I) instead of current density (J)
        plot(mode.vv_dd+mode.ii_dd*Rs,mode.ii_dd*1e3,'k--','LineWidth',2)
        ylabel('Current, mA')
        ylim([0.1 90])
        set(gca,'FontSize',14,'FontName','Arial')
    else
        plot(mode.vv_dd+mode.ii_dd*Rs,mode.ii_dd/mode.AreaOx,'LineWidth',2)
        axis([0 max(mode.vv_dd+mode.ii_dd*Rs) 1e2 1e5])
        ylabel('Current density, A/cm^2')
        set(gca,'yscale','log','FontSize',14,'FontName','Arial')
    end
% else
%     plot(mode.vv_dd,mode.ii_dd/mode.AreaOx,'.-','linewidth',2)
end
xlim([0 max(mode.vv_dd+mode.ii_dd*Rs)])
xlabel('Voltage, V')
set(gcf,'Position',[680    51   557   420])

%% Band diagram
ecb=squeeze(MODEplot{1}.ecb(end,:,end));
evb=squeeze(MODEplot{1}.evb(end,:,end));
EFn=squeeze(MODEplot{1}.EFn(end,:,end));
EFp=squeeze(MODEplot{1}.EFp(end,:,end));
z=mesh.ygrid*1e4;
figure,hold on,grid on
plot(z,ecb,'b.-','linewidth',1)
plot(z,evb,'r.-','linewidth',1)
plot(z,EFn,'k--','linewidth',1)
plot(z,EFp,'k-.','linewidth',1)
xlabel('z, \mum'),ylabel('Energy band diagram, eV')
title(['V = ',num2str(mode.vv_dd(end)),' V'])
