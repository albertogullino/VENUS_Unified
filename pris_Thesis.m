clear
close all
dbstop if error

if exist('CORRENTI')==0
    CORRENTI=[2 8 10];
    CORRENTI=[1 4 8];
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

T0=P.modeplot.T0

load(['MarkusN_4_T',num2str(double(T0)),'.mat'])
pa=1:10:length(Imeas);

%% IV
figure(80)
hold on,grid on,box on
plot(P.modeplot.vv_dd,P.modeplot.ii_dd*1000,'r.-','LineWidth',2)
plot(A.modeplot.vv_dd,A.modeplot.ii_dd*1000,'g.-','LineWidth',2)
plot(B.modeplot.vv_dd,B.modeplot.ii_dd*1000,'b.-','LineWidth',2)
plot(Vmeas(pa),Imeas(pa),'r.','markersize',12)

set(gca,'FontSize',14,'FontName','Times new roman')
axis([1.4 3.2 0 15])
xticks([1.4:.3:3.2])
xlabel('Voltage, V'),ylabel('Current, mA')
legend('\it pin','TJ - Ox. Above','TJ - Ox. Below','location','northwest')

%% LI
figure(81)
hold on,grid on,box on
plot(P.modeplot.ii_dd*1000,sum(P.modeplot.Pst_dd,1),'r.-','LineWidth',2)
plot(A.modeplot.ii_dd*1000,sum(A.modeplot.Pst_dd,1),'g.-','LineWidth',2)
plot(B.modeplot.ii_dd*1000,sum(B.modeplot.Pst_dd,1),'b.-','LineWidth',2)
plot(Imeas(pa),Lmeas(pa),'r.','LineWidth',2)

plot(P.modeplot.ii_dd*1000,P.modeplot.Pst_dd(1,:),'ro','LineWidth',1,'markersize',4)
plot(P.modeplot.ii_dd*1000,P.modeplot.Pst_dd(2,:),'rv','LineWidth',1,'markersize',4)
plot(P.modeplot.ii_dd*1000,P.modeplot.Pst_dd(3,:),'r*','LineWidth',.5,'markersize',6)

plot(A.modeplot.ii_dd*1000,A.modeplot.Pst_dd(1,:),'go','LineWidth',1,'markersize',4)
plot(A.modeplot.ii_dd*1000,A.modeplot.Pst_dd(2,:),'gv','LineWidth',1,'markersize',4)
plot(A.modeplot.ii_dd*1000,A.modeplot.Pst_dd(3,:),'g*','LineWidth',.5,'markersize',6)

plot(B.modeplot.ii_dd*1000,B.modeplot.Pst_dd(1,:),'bo','LineWidth',1,'markersize',4)
plot(B.modeplot.ii_dd*1000,B.modeplot.Pst_dd(2,:),'bv','LineWidth',1,'markersize',4)
plot(B.modeplot.ii_dd*1000,B.modeplot.Pst_dd(3,:),'b*','LineWidth',.5,'markersize',6)

set(gca,'FontSize',16,'FontName','Times new roman')
axis([0 16 0 6.5])
xticks([0:2:16])
ylabel('Output optical power, mW'),xlabel('Current, mA')
legend('\it pin','TJ - Ox. Above','TJ - Ox. Below','location','northwest')

%% WPE
figure(82)
hold on,grid on,box on
plot(P.modeplot.ii_dd*1e3,sum(P.modeplot.Pst_dd,1)./(P.modeplot.ii_dd*1000.*P.modeplot.vv_dd)*100,'r.-','LineWidth',2)
plot(A.modeplot.ii_dd*1e3,sum(A.modeplot.Pst_dd,1)./(A.modeplot.ii_dd*1000.*A.modeplot.vv_dd)*100,'g.-','LineWidth',2)
plot(B.modeplot.ii_dd*1e3,sum(B.modeplot.Pst_dd,1)./(B.modeplot.ii_dd*1000.*B.modeplot.vv_dd)*100,'b.-','LineWidth',2)
plot(Imeas(pa),Lmeas(pa)./(Vmeas(pa).*Imeas(pa))*100,'r.','markersize',12)

set(gca,'FontSize',14,'FontName','Times new roman')
axis([0 14 0 30])
ylabel('\eta_{WP}, %'),xlabel('Current, mA')
legend('\it pin','TJ - Ox. Above','TJ - Ox. Below','location','south')

%% Wavelength
figure(83)
hold on,grid on,box on
plot(P.modeplot.ii_dd*1000,P.modeplot.lambda(1,:),'r','LineWidth',1.5)
plot(A.modeplot.ii_dd*1000,A.modeplot.lambda(1,:),'g','LineWidth',1.5)
plot(B.modeplot.ii_dd*1000,B.modeplot.lambda(1,:),'b','LineWidth',1.5)
plot(Cur,LAM,'ro','LineWidth',1,'markersize',4)

thP11=find(P.modeplot.ii_dd*1e3>3);
thA11=find(A.modeplot.ii_dd*1e3>3);
thB11=find(B.modeplot.ii_dd*1e3>3);
plot(P.modeplot.ii_dd(thP11:end)*1000,P.modeplot.lambda(2,thP11:end),'r--','LineWidth',1.5)
plot(A.modeplot.ii_dd(thA11:end)*1000,A.modeplot.lambda(2,thA11:end),'g--','LineWidth',1.5)
plot(B.modeplot.ii_dd(thB11:end)*1000,B.modeplot.lambda(2,thB11:end),'b--','LineWidth',1.5)

thP02=find(P.modeplot.ii_dd*1e3>8);
thA02=find(A.modeplot.ii_dd*1e3>8);
thB02=find(B.modeplot.ii_dd*1e3>8);
plot(P.modeplot.ii_dd(thP02:end)*1000,P.modeplot.lambda(3,thP02:end),'r:','LineWidth',1.5)
plot(A.modeplot.ii_dd(thA02:end)*1000,A.modeplot.lambda(3,thA02:end),'g:','LineWidth',1.5)
plot(B.modeplot.ii_dd(thB02:end)*1000,B.modeplot.lambda(3,thB02:end),'b:','LineWidth',1.5)

set(gca,'FontSize',14,'FontName','Times new roman')
axis([0.5 14 847 852])

xlabel('Current, mA'),ylabel('Wavelength, nm')
legend('\it pin','TJ - Ox. Above','TJ - Ox. Below','location','northwest')


%% Leakage current
figure(84)
hold on,grid on,box on
plot(P.modeplot.ii_dd*1000,P.modeplot.ii_dd.*1000.*P.modeplot.Fleak,'r.-','LineWidth',2)
plot(A.modeplot.ii_dd*1000,A.modeplot.ii_dd.*1000.*A.modeplot.Fleak,'g.-','LineWidth',2)
plot(B.modeplot.ii_dd*1000,B.modeplot.ii_dd.*1000.*B.modeplot.Fleak,'b.-','LineWidth',2)

Curr=10; % current value where plot a black dot in the leakage plot
tolCurr=1e-4; % tolerance of the find in mode.ii_dd
indIP=find(abs(P.modeplot.ii_dd*1e3-Curr)<tolCurr);
indIA=find(abs(A.modeplot.ii_dd*1e3-Curr)<tolCurr);
indIB=find(abs(B.modeplot.ii_dd*1e3-Curr)<tolCurr);

plot(P.modeplot.ii_dd(indIP)*1000,P.modeplot.ii_dd(indIP).*1000.*P.modeplot.Fleak(indIP),'kd','LineWidth',1)
plot(A.modeplot.ii_dd(indIA)*1000,A.modeplot.ii_dd(indIA).*1000.*A.modeplot.Fleak(indIA),'kd','LineWidth',1)
plot(B.modeplot.ii_dd(indIB)*1000,B.modeplot.ii_dd(indIB).*1000.*B.modeplot.Fleak(indIB),'kd','LineWidth',1)

axis([0 14 0 3])
legend('\it pin','TJ - Ox. Above','TJ - Ox. Below','location','northwest')
xlabel('Total current, mA'),ylabel('Leakage current, mA')
set(gca,'FontSize',14,'FontName','Times new roman')
set(gcf,'Position',[268   441   695   496])
% Inset with leak cur percentage
ax2 = axes('Box','on');
set(gca,'position',[.22 .3 .4 .4])
ylabel('Leak. cur. %')
hold on, box on, grid on
plot(P.modeplot.ii_dd*1000,P.modeplot.Fleak*100,'r.-','LineWidth',2)
plot(A.modeplot.ii_dd*1000,A.modeplot.Fleak*100,'g.-','LineWidth',2)
plot(B.modeplot.ii_dd*1000,B.modeplot.Fleak*100,'b.-','LineWidth',2)

plot(P.modeplot.ii_dd(indIP)*1000,P.modeplot.Fleak(indIP)*100,'kd','LineWidth',1)
plot(A.modeplot.ii_dd(indIA)*1000,A.modeplot.Fleak(indIA)*100,'kd','LineWidth',1)
plot(B.modeplot.ii_dd(indIB)*1000,B.modeplot.Fleak(indIB)*100,'kd','LineWidth',1)

set(gca,'FontSize',10,'FontName','Times new roman')
axis([1 14 0 25])

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
COL='rgbcmyk';

for iI=[1 3]
%     coloriWhite
%     figure(1185+iCurrP(iI))
    set(gcf,'pos',[175          77        1453         902])
    subplot(221)
    hold on, box on
%     plot(xqwP,squeeze(P.MODEplot{1}.E2(iCurrP(iI),:,1:nnQWP)),'r.-')  
%     plot(xqwA,squeeze(A.MODEplot{1}.E2(iCurrA(iI),:,1:nnQWA)),'g.-')
%     plot(xqwB,squeeze(B.MODEplot{1}.E2(iCurrB(iI),:,1:nnQWB)),'b.-')
    plot(xqwP,squeeze(P.modeplot.E2(iCurrP(iI),:,1:nnQWP)),[COL(1),lines{iI}],'linewidth',LW(iI))  
%     plot(xqwA,squeeze(A.modeplot.E2(iCurrA(iI),:,1:nnQWA)),[COL(2),lines{iI}],'linewidth',LW(iI))
    plot(xqwB,squeeze(B.modeplot.E2(iCurrB(iI),:,1:nnQWB)),[COL(3),lines{iI}],'linewidth',LW(iI))
    xlim([0 5]),grid on
    xlabel('\rho, \mum'),ylabel('Optical field intensity, a.u.')
    pausak
%     colororder(newcolors)
%     title('E2')
    
%     subplot(221)
%     hold on,box on
%     plot(xqwP,P.MODEplot{1}.matgain(iCurrP(iI),:).*squeeze(P.MODEplot{1}.E2(iCurrP(iI),:,1:nnQWP)),'r','linewidth',1.5)
%     plot(xqwA,A.MODEplot{1}.matgain(iCurrA(iI),:).*squeeze(A.MODEplot{1}.E2(iCurrA(iI),:,1:nnQWA)),'g','linewidth',1.5)
%     plot(xqwB,B.MODEplot{1}.matgain(iCurrB(iI),:).*squeeze(B.MODEplot{1}.E2(iCurrB(iI),:,1:nnQWB)),'b','linewidth',1.5)
%     xlim([0 5]),grid on
%     xlabel('\rho, \mum'),ylabel('Gain-field product')
%     title('matgain \cdot E2')
    title(['I=',num2str(CORRENTI(iI)),' mA'])
    set(gca,'FontSize',16,'FontName','Times new roman')
    
    subplot(222)
    hold on,box on
    plot(xqwP,P.modeplot.matgain(iCurrP(iI),:),[COL(1),lines{iI}],'linewidth',LW(iI))
    plot(xqwA,A.modeplot.matgain(iCurrA(iI),:),[COL(2),lines{iI}],'linewidth',LW(iI))
    plot(xqwB,B.modeplot.matgain(iCurrB(iI),:),[COL(3),lines{iI}],'linewidth',LW(iI))
    xlim([0 5]),grid on
    xlabel('\rho, \mum')
    ylabel('QW gain, cm^{-1}')
%     title('matgain')
    set(gca,'FontSize',16,'FontName','Times new roman')
    
    subplot(223)
    hold on,box on
    plot(xqwP,squeeze(P.modeplot.Temp(iCurrP(iI),iQWP,1:nnQWP))-(T0+273),[COL(1),lines{iI}],'linewidth',LW(iI))
    plot(xqwA,squeeze(A.modeplot.Temp(iCurrA(iI),iQWA,1:nnQWA))-(T0+273),[COL(2),lines{iI}],'linewidth',LW(iI))
    plot(xqwB,squeeze(B.modeplot.Temp(iCurrB(iI),iQWB,1:nnQWB))-(T0+273),[COL(3),lines{iI}],'linewidth',LW(iI))
    xlim([0 5]),grid on
    xlabel('\rho, \mum'), ylabel('\DeltaT radial profile, °C')
    set(gca,'FontSize',16,'FontName','Times new roman')
    
        ax2=axes;
        set(gca,'position',[0.3008 0.3237 0.1942 0.1840])
        hold on
        plot(yP,squeeze(P.modeplot.Temp(iCurrP(iI),:,1))-(T0+273),[COL(1),lines{iI}],'linewidth',LW(iI))
        plot(yA,squeeze(A.modeplot.Temp(iCurrA(iI),:,1))-(T0+273),[COL(2),lines{iI}],'linewidth',LW(iI))
        plot(yB,squeeze(B.modeplot.Temp(iCurrB(iI),:,1))-(T0+273),[COL(3),lines{iI}],'linewidth',LW(iI))
        xlabel('z, \mum'), xlim([110 118])
        set(gca,'FontSize',12,'FontName','Times new roman')

    
    subplot(224)
    hold on,box on
        plot(xqwP,P.modeplot.nQW{iCurrP(iI)}{2},[COL(1),lines{iI}],'linewidth',LW(iI)) % indicates the central QW!
        plot(xqwA,A.modeplot.nQW{iCurrA(iI)}{2},[COL(2),lines{iI}],'linewidth',LW(iI)) % indicates the central QW!
        plot(xqwB,B.modeplot.nQW{iCurrB(iI)}{2},[COL(3),lines{iI}],'linewidth',LW(iI)) % indicates the central QW!
%     ResetColor
        plot(xqwP,P.modeplot.pQW{iCurrP(iI)}{2},[COL(1),lines{iI+4}],'linewidth',LW(iI))
        plot(xqwA,A.modeplot.pQW{iCurrA(iI)}{2},[COL(2),lines{iI+4}],'linewidth',LW(iI))
        plot(xqwB,B.modeplot.pQW{iCurrB(iI)}{2},[COL(3),lines{iI+4}],'linewidth',LW(iI))
    ylabel(' Sheet carrier  density, 1/cm^2 ')
    xlabel(' \rho, \mum')
    
        xlim([0 5]),grid on
        set(gca,'FontSize',16,'FontName','Times new roman')
end
subplot(224)
legend('n^{2D} (1 mA)','','','p^{2D}','','','n^{2D} (10 mA)','','','p^{2D}')

% figure
% hold on,box on
% plot(xqwP,squeeze(P.MODEplot{1}.Temp(iCurrP,iQWP,1:nnQWP)),COL(1),'linewidth',2)
% plot(xqwA,squeeze(A.MODEplot{1}.Temp(iCurrA,iQWA,1:nnQWP)),COL(2),'linewidth',2)
% plot(xqwB,squeeze(B.MODEplot{1}.Temp(iCurrB,iQWB,1:nnQWP)),COL(3),'linewidth',2)
% xlim([0 5]),grid on
% xlabel('\rho, \mum'), ylabel('Radial T profile, K')

RdiffP = [NaN, diff(P.mode.vv_dd)./diff(P.mode.ii_dd)];
RdiffA = [NaN, diff(A.mode.vv_dd)./diff(A.mode.ii_dd)];
RdiffB = [NaN, diff(B.mode.vv_dd)./diff(B.mode.ii_dd)];

pa=1:3:length(Imeas)-300;
Rmeas=diff(Vmeas)./diff(Imeas)*1e3;
Rmeas=[NaN; Rmeas];
figure
hold on,grid on,box on
plot(Imeas(pa),Rmeas(pa),'r.','LineWidth',1.5,'MarkerSize',12)
plot(P.mode.ii_dd*1e3,RdiffP,'r.-','LineWidth',2)
plot(A.mode.ii_dd*1e3,RdiffA,'g.-','LineWidth',2)
plot(B.mode.ii_dd*1e3,RdiffB,'b.-','LineWidth',2)
axis([0 15 20 300])
xlabel('Current, mA')

% plot(Vmeas(pa),Rmeas(pa),'r.','LineWidth',1.5,'MarkerSize',12)
% plot(P.mode.vv_dd,RdiffP,'r.-','LineWidth',2)
% plot(A.mode.vv_dd,RdiffA,'g.-','LineWidth',2)
% plot(B.mode.vv_dd,RdiffB,'b.-','LineWidth',2)
% axis([1.4 3.2 20 300])
% xticks([1.4:.3:3.2])
% xlabel('Voltage, V')

ylabel('R_{diff}, \Omega')
set(gca,'FontSize',18,'FontName','Times new roman')
legend('','\it pin','TJ - Ox. Above','TJ - Ox. Below','location','northeast')

iRagCheck

figure(80)
plot(A.modeplot.vv_dd-VrhoA',A.modeplot.ii_dd*1000,'g--','LineWidth',2)
plot(B.modeplot.vv_dd-VrhoB'+.1,B.modeplot.ii_dd*1000,'b--','LineWidth',2)

HeatSourcePlot(P.mesh,P.mode,P.MODEplot{end})
HeatSourcePlot(A.mesh,A.mode,A.MODEplot{end})
HeatSourcePlot(B.mesh,B.mode,B.MODEplot{end})

