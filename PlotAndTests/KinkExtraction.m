% Kink voltage extraction
%
clear 

s_LoadConstants
h2m0=hbar^2/(2*m0);

Imin=0.3;
Imax=3;

%% Load the workspace 
% Barrier xmol (xB=0.1)
load('out\LW_MarkusN_FINALE_LgDBR_1D_newLUT_xB=0.1.mat')
[xB1,xQW1,Vk1,EgB1,EgQW1]=LUT_plotter(mode,geom,Imin,Imax);

% (xB=0.2)
load('out\LW_MarkusN_FINALE_LgDBR_1D_newLUT_xB=0.2.mat')
[xB2,xQW2,Vk2,EgB2,EgQW2]=LUT_plotter(mode,geom,Imin,Imax);

% (xB=0.286)
load('out\LW_MarkusN_FINALE_LgDBR_1D_newLUT.mat')
[xB3,xQW3,Vk3,EgB3,EgQW3]=LUT_plotter(mode,geom,Imin,Imax);

% load('out\LW_MarkusN_FINALE_LgDBR_1D_newLUT_bulk.mat')
% [xB,xQW,Vk,EgB,EgQW]=LUT_plotter(mode,geom,Imin,Imax);

% (xB=0.4)
load('out\LW_MarkusN_FINALE_LgDBR_1D_newLUT_xB=0.4.mat')
[xB4,xQW4,Vk4,EgB4,EgQW4]=LUT_plotter(mode,geom,Imin,Imax);

% QW xmol (xQW=0.1 - xB=0.4)
load('out\LW_MarkusN_FINALE_LgDBR_1D_newLUT_xB=0.4_xQW=0.1.mat')
[xB41,xQW41,Vk41,EgB41,EgQW41]=LUT_plotter(mode,geom,Imin,Imax);

% (xQW=0.2 - xB=0.4)
load('out\LW_MarkusN_FINALE_LgDBR_1D_newLUT_xB=0.4_xQW=0.2.mat')
[xB42,xQW42,Vk42,EgB42,EgQW42]=LUT_plotter(mode,geom,Imin,Imax);

pausak

% xB variation
xB=[xB1 xB2 xB3 xB4];
VkB=[Vk1 Vk2 Vk3 Vk4];

% xQW variation
xQW=[xQW4 xQW41 xQW42]
VkQW=[Vk4 Vk41 Vk42]

% Subbands bandgap = bottom CB - top VB
% Values extracted graphically from the plots "VB_Holger" and "CB_Holger":
% it is possible to compute them (and plot the same figures) by loading the
% results from the LUT (stored in the corresponding folder)
Eg2D=[EgQW1 EgQW2 EgQW3 EgQW4]; % Subbands bandgap (xB variation)
EgB=[EgB1 EgB2 EgB3 EgB4]; % Barriers bandgap
%
Eg2Dqw=[EgQW4 EgQW41 EgQW42]; % Subbands bandgap (xQW variation)

%% Figures
% xB variation
figure
yyaxis left
grid on
plot(xB,VkB,'o-')
ylabel('Kink voltage, V')
xlabel('Barrier molar fraction (xB)')
yyaxis right
hold on
plot(xB,Eg2D,'*-')
plot(xB,EgB,'*--')
legend('','subbands','barrier','location','southeast')
ylabel('Energy bandgap, eV')
set(gca,'FontSize',14,'FontName','Arial','box','on')

% xQW variation
figure
yyaxis left
grid on
plot(xQW,VkQW,'o-')
ylabel('Kink voltage, V')
xlabel('QW molar fraction (xQW)')
yyaxis right
hold on
plot(xQW,Eg2Dqw,'*-')
% plot(xB,EgB,'*--')
% legend('','subbands','barrier','location','southeast')
ylabel('Energy bandgap, eV')
set(gca,'FontSize',14,'FontName','Arial','box','on')

% as a function of Eg2D
figure
hold on,grid on
plot(Eg2D,VkB,'ko-','LineWidth',2)
plot(Eg2Dqw,VkQW,'m*-','LineWidth',2)
legend('xB variation','xQW variation','location','southeast')
set(gca,'FontSize',14,'FontName','Arial','box','on')
legend('xB=0.1->0.4 (xQW=0)','xQW=0->0.2 (xQW=0.4)','location','southeast')
xlabel('E_g^{QW}, eV'),ylabel('Kink voltage, V')
