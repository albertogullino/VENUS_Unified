clear
% % close all
clc

% addpathVENUS

%% Load NEGF-DD results from original folder (no Rs inserted)
load IV_noR_2020Tibaldi_PRAPP.mat

dox = 18; % DIFFERENT FROM 17um OF THE EXPERIMENTAL DATAAAAAAAA
Area = pi*(dox/2*1e-4)^2; % cm^2
R = 6; % ohm*cm^2, series resistance to fit the experimental curve

Ibias = Jbias*Area;
Vtot = abs(Vbias) + R*Ibias;

%% Chalmers experimental curve, used in 2020Tibaldi_PRAPP: Rox=8.5 um
load p_on_n_dox=17.mat

Vmeas=Vvet;
Imeas=Ivet;

%% Plot section

% I(V) plot
figure(43)
set(gcf,'Position',[1097 283 560 420])
hold on
box on
grid on
[AxisSim,HandleSimVI,HandleSimRI] = plotyy(1000*Ibias,Vtot,1000*Ibias,gradient(Vtot)./gradient(Ibias));
set(AxisSim(1),'XLim',[1,90])
set(AxisSim(1),'YLim',[0,2.5])
set(AxisSim(1),'YTick',[])
set(HandleSimVI,'LineWidth',1.5)
set(HandleSimVI,'Color',[0 0 1])
set(AxisSim(1),'YColor',[0 0 1])
%
set(AxisSim(2),'XLim',[1,90])
set(AxisSim(2),'YLim',[0,40])
set(AxisSim(2),'YTick',[])
set(HandleSimRI,'LineWidth',1.5)
set(HandleSimRI,'Color',[1 0 0])
set(AxisSim(2),'YColor',[1 0 0])
%
xlabel(AxisSim(1),'Current, mA','FontSize',14,'FontName','Times New Roman')
%
[AxisMeas,HandleMeasVI,HandleMeasRI] = plotyy(1000*Imeas,Vmeas,1000*Imeas,gradient(Vmeas)./gradient(Imeas));
legend('Simulations','Measurements','Location','Best')

set(AxisMeas(1),'XLim',[1,90])
set(AxisMeas(1),'YLim',[0,2.5])
set(AxisMeas(1),'YTick',[0:0.5:2.5],'FontSize',13,'FontName','Times New Roman')
set(HandleMeasVI,'LineStyle','none')
set(HandleMeasVI,'LineWidth',1.5)
set(HandleMeasVI,'Marker','o')
set(HandleMeasVI,'Color',[0 0 1])
set(AxisMeas(1),'YColor',[0 0 1])
ylabel(AxisMeas(1),'Voltage, V','FontSize',14,'FontName','Times New Roman')
%
set(AxisMeas(2),'XLim',[1,90])
set(AxisMeas(2),'YLim',[0,40])
set(AxisMeas(2),'YTick',[0:8:40],'FontSize',13,'FontName','Times New Roman')
set(HandleMeasRI,'LineStyle','none')
set(HandleMeasRI,'LineWidth',1.5)
set(HandleMeasRI,'Marker','o')
set(HandleMeasRI,'Color',[1 0 0])
set(AxisMeas(2),'YColor',[1 0 0])
ylabel(AxisMeas(2),'Differential resistance, \Omega','FontSize',14,'FontName','Times New Roman')