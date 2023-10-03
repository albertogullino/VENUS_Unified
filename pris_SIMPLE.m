addpathVENUS
colordef white 

clear Imeas Lmeas Rmeas Vmeas T0

T0=mode.T0-273;

% if mode.flgBTJ==0
    load(['MarkusN_4_T',num2str(double(T0)),'.mat'])
    pa=1:10:length(Imeas);
% else
%     load VELM3-TJ
%     pa=1:length(Imeas);
% end

vind=[];
for indplot=1:length(VELMInfo)
    vind=[vind,VELMInfo(indplot).indVoltage];
end

figure(80)
hold on
grid on
box on
plot(Vmeas(pa),Imeas(pa),'ro','markersize',4,'linewidth',1)

plot(mode.vv_dd,mode.ii_dd*1000,'b.-','LineWidth',2,'markersize',8)
% plot(mode.vv_dd(vind),mode.ii_dd(vind)*1000,'bo','markersize',4)
axis([1.4 mode.vv_dd(end)+.2 0 1000*mode.ii_dd(end)+1])
xlabel('Voltage, V')
ylabel('Current, mA')

figure(81)
hold on
grid on
box on
plot(Imeas(pa),Lmeas(pa),'ro','markersize',4,'linewidth',1)
chold
plot(Cur,10.^(P_dB/10),'s--','LineWidth',2)
PPst=sum(mode.Pst_dd,1)+mode.Psp_dd;

plot(mode.ii_dd*1000,PPst,'b.-','LineWidth',2,'markersize',8)
chold
plot(mode.ii_dd*1000,mode.Pst_dd+mode.Psp_dd,'.','markersize',6)

axis([0 1000*mode.ii_dd(end)+1 0 max(sum(mode.Pst_dd,1))*1.1+.1])
xlabel('Current, mA')
ylabel('Optical power, mW')

figure(82)
hold on,grid on,box on
plot(Cur,LAM,'ro','linewidth',1,'markersize',4)
chold
plot(mode.ii_dd*1e3,mode.lambda,'b.-','linewidth',1.5)
xlabel('Current, mA')
ylabel('Emission wavelength, nm')
xlim([0 mode.ii_dd(end)*1e3+0.01])
title('VCSEL thermometer')

figure
                set(gcf,'Position',[  634    73   579   365])
                hold on
                grid on
                box on
                plot(mode.ii_dd*1e3,mode.DeltaTmax,'k','LineWidth',2)
                plot(mode.ii_dd*1e3,mode.DeltaTmax_Joule,'c--','LineWidth',2)
                plot(mode.ii_dd*1e3,mode.DeltaTmax_OptAbs,'g--','LineWidth',2)
                plot(mode.ii_dd*1e3,mode.DeltaTmax_Ccap,'r--','LineWidth',2)
                plot(mode.ii_dd*1e3,mode.DeltaTmax_RAD,'b--','LineWidth',2)
                plot(mode.ii_dd*1e3,mode.DeltaTmax_srhAu,'y--','LineWidth',2)

                legend('\DeltaT','Joule','Opt. abs.','Ccap','Rad','SRH/Aug','location','best')
                xlabel('Current, mA')
                ylabel('Temperature rise, K')
				
				
            figure(120)
            hold on
            grid on
            box on
            plot(mode.ii_dd*1e3,mode.Gmod,'k','LineWidth',2)
            plot(mode.ii_dd*1e3,mode.Lmod,'r--','LineWidth',2)
            if(mode.nmodes>1)
                plot(mode.ii_dd(vind)*1e3,mode.Gmod(:,vind)','ko','LineWidth',2)
                plot(mode.ii_dd(vind)*1e3,mode.Lmod(:,vind)','ro','LineWidth',2)
            else
                plot(mode.ii_dd(vind)*1e3,mode.Gmod(vind),'ko','LineWidth',2)
                plot(mode.ii_dd(vind)*1e3,mode.Lmod(vind),'ro','LineWidth',2)
            end
            xlabel('Current, mA')
            ylabel('Modal gain vs losses, cm^{-1}')		
ylim([min(mode.Lmod(1,:)) max(mode.Lmod(1,:)) ])			