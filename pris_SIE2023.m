clear
close all
dbstop if error

if exist('CORRENTI')==0
    CORRENTI=[1 6 10];
end

addpathVENUS


% P=load('out\LW_MarkusN_FINALE_LgDBR_NUSOD2023.mat');
% A=load('out\LW_MarkusN_TJ_oxAbove_LgDBR_FINAL.mat');
% B=load('out\LW_MarkusN_TJ_oxBelow_LgDBR_NUSOD2023.mat');

P=load('out\LW_MarkusN_FINALE_LgDBR_20C_radial.mat');
A=load('out\LW_MarkusN_TJ_oxAbove_LgDBR_20C_lambda_radial.mat');
B=load('out\LW_MarkusN_TJ_oxBelow_LgDBR_20C_lambda_radial.mat');


if length(P.MODEplot)>1
    iT=input('iT?\n');
    P.modeplot=P.MODEplot{iT};
    A.modeplot=A.MODEplot{iT};
    B.modeplot=B.MODEplot{iT};
else
    P.modeplot=P.MODEplot{end};
    A.modeplot=A.MODEplot{end};
    B.modeplot=B.MODEplot{end};
end


COL='rgbcm';
newcolors={'r','#FF8A00','b'};
COL=newcolors;

T0=P.modeplot.T0

load(['MarkusN_4_T',num2str(double(T0)),'.mat'])
pa=1:10:length(Imeas);

%% IV
figure(80)
hold on,grid on,box on
plot(P.modeplot.vv_dd,P.modeplot.ii_dd*1000,'color',COL{1},'LineWidth',2)
plot(A.modeplot.vv_dd,A.modeplot.ii_dd*1000,'color',COL{2},'LineWidth',2)
plot(B.modeplot.vv_dd,B.modeplot.ii_dd*1000,'color',COL{3},'LineWidth',2)
plot(Vmeas(pa),Imeas(pa),'r.','markersize',12)

set(gca,'FontSize',16,'FontName','Times new roman')
axis([1.4 3.2 0 15])
xticks([1.4:.3:3.2])
xlabel('Voltage, V'),ylabel('Current, mA')
legend('\it pin','TJ - Ox. Above','TJ - Ox. Below','location','northwest')

%% LI
figure(81)
hold on,grid on,box on
plot(P.modeplot.ii_dd*1000,sum(P.modeplot.Pst_dd,1),'color',COL{1},'LineWidth',2)
plot(A.modeplot.ii_dd*1000,sum(A.modeplot.Pst_dd,1),'color',COL{2},'LineWidth',2)
plot(B.modeplot.ii_dd*1000,sum(B.modeplot.Pst_dd,1),'color',COL{3},'LineWidth',2)
plot(Imeas(pa),Lmeas(pa),'r.','LineWidth',2)

set(gca,'FontSize',16,'FontName','Times new roman')
axis([0 16 0 6.5])
xticks([0:2:16])
ylabel('Output optical power, mW'),xlabel('Current, mA')
legend('\it pin','TJ - Ox. Above','TJ - Ox. Below','location','northwest')

%% Gain-field product
nnQWP=P.mesh.nnxQW{1};
xqwP=P.mesh.xgrid(1:nnQWP)*1e4;
iQWP=P.mesh.inMQW{2}(1);
yP=P.mesh.ygrid*1e4;
nnQWA=A.mesh.nnxQW{1};
xqwA=A.mesh.xgrid(1:nnQWA)*1e4;
iQWA=A.mesh.inMQW{2}(1);
yA=A.mesh.ygrid*1e4;
nnQWB=B.mesh.nnxQW{1};
xqwB=B.mesh.xgrid(1:nnQWB)*1e4;
iQWB=B.mesh.inMQW{2}(1);
yB=B.mesh.ygrid*1e4;

iCurrP=CurrIndex(CORRENTI,P.modeplot.ii_dd*1e3);
iCurrA=CurrIndex(CORRENTI,A.modeplot.ii_dd*1e3);
iCurrB=CurrIndex(CORRENTI,B.modeplot.ii_dd*1e3);

figure(1185+iCurrP(1))

lines{1}='-'; lines{2}='-.'; lines{3}='--';  lines{4}='.';  lines{5}='s';  lines{6}='o';  lines{7}='d'; lines{8}='*';
LW=[1 1.8 2.5];     % linewidth vector
% LW=[1 2.5];     % linewidth vector
% COL='rgbcmyk';

for iI=[1 3]
    
    hold on,box on
        plot(xqwP,P.modeplot.nQW{iCurrP(iI)}{2},lines{iI},'color',COL{1},'linewidth',LW(iI)) % indicates the central QW!
        plot(xqwA,A.modeplot.nQW{iCurrA(iI)}{2},lines{iI},'color',COL{2},'linewidth',LW(iI)) % indicates the central QW!
        plot(xqwB,B.modeplot.nQW{iCurrB(iI)}{2},lines{iI},'color',COL{3},'linewidth',LW(iI)) % indicates the central QW!
%     ResetColor
        plot(xqwP,P.modeplot.pQW{iCurrP(iI)}{2},lines{iI+4},'color',COL{1},'linewidth',LW(iI))
        plot(xqwA,A.modeplot.pQW{iCurrA(iI)}{2},lines{iI+4},'color',COL{2},'linewidth',LW(iI))
        plot(xqwB,B.modeplot.pQW{iCurrB(iI)}{2},lines{iI+4},'color',COL{3},'linewidth',LW(iI))
    ylabel(' Sheet carrier  density, 1/cm^2 ')
    xlabel(' \rho, \mum')
    
        xlim([0 5]),grid on
        set(gca,'FontSize',16,'FontName','Times new roman')
end
legend('n^{2D} (1 mA)','','','p^{2D}','','','n^{2D} (10 mA)','','','p^{2D}')



