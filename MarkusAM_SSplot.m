close all
% clear

if mode.SmallSignal == 1
    wafer = '01103_448C';
    device = 'C134R233';    % available: C134R233, C134R235, C134R237
    
    Tvet = [20];            % available: 20, 30, 40, 50, 60
    %     Ivet = [1:6];           % available: 1, 2, 3, 4, 5, 6
    Ivet = [1 2 3];           % available: 1, 2, 3, 4, 5, 6
    Ivet = [1 2];           % available: 1, 2, 3, 4, 5, 6
    
    Z0 = 50;
    colorMatrix = [[0 0.4470 0.7410],
        [0.8500 0.3250 0.0980],
        [0.9290 0.6940 0.1250],
        [0.4940 0.1840 0.5560],
        [0.4660 0.6740 0.1880],
        [0.3010 0.7450 0.9330],
        [0.6350 0.0780 0.1840]];
    %
    for indT = 1:length(Tvet)
        for indI = 1:length(Ivet)
            fileName = ['C:\Users\albig\Politecnico\Dottorato\3b_VENUS\SmallSignalJul19\',...
                wafer,'_',num2str(Tvet(indT)),'\',wafer,'_',device,'_I_',num2str(Ivet(indI)),'mA_',num2str(Tvet(indT)),'°C.s2p'];
            Sstr = sparameters(fileName);
            vFreq = Sstr.Frequencies;
            
            
            % Scattering parameters extraction
            S21_exp = squeeze(Sstr.Parameters(2,1,:));
            S11_exp = squeeze(Sstr.Parameters(1,1,:));
            
            Z_exp = Z0*(1+S11_exp)./(1-S11_exp);
            
            figure(21)%,subplot(1,3,3)
            hold on
%             semilogx(vFreq/1e9,10*log10(abs(S21_exp/S21_exp(1)).^2),'-.','color',colorMatrix(indI,:),'linewidth',2)
            semilogx(vFreq/1e9,10*log10(abs(S21_exp/S21_exp(1)).^2),'o','linewidth',2)
%             semilogx(vFreq/1e9,10*log10(abs(S21_exp/S21_exp(1))),'o','linewidth',2)
            xlim([1 100])
        end
    end
    %
    
    grid on,box on
    % xlim([0 20])
    xlabel('Frequency, GHz')
    ylabel('AM response, normalized')
    
    %% 1D simulation
%     strFileName='MarkusN_TopLAqw';
    strFileName='NUSOD_D1ANA';
%     load(['results_',strFileName,'.mat'])
    %
    CutoffLine = -3*ones(1,length(mode.fvet));
    
    f = mode.fvet;
    Rm = 4;
    Cp = 15.5e-12;
    fRC = 1/(2*pi*Rm*Cp)*1e-9
    
%     Rm = 0;
%     Cp = 0;
    figure(21)
    set(gcf,'Position',[1034 513 560 420])
    
    p = 1;
    for i = length(mode_ss)-20:4:length(mode_ss)
        Iplot(p) = mode.I_dd(i)*mode.Area*1e3;
        
        Ypad = 1./((1./mode_ss(i).Y)+Rm);
        
        Zpar = Rm-1i./(2*pi.*f.*Cp);
        Ypar=1./Zpar;
        
        Ytot = Ypad + 1i*2*pi*f*Cp;
        AM = abs(mode_ss(i).Pst./Ytot).^2;
%         AM = abs(mode_ss(i).Pst./mode_ss(i).Y).^2;
        AMy= abs(1./mode_ss(i).Y).^2;
        AMypar=abs(1./Ypar).^2;
        AMytot= abs(1./Ytot).^2;
        AMp= abs(mode_ss(i).Pst).^2;
        figure(21)
        semilogx(f/1e9,10*log10(AM./AM(1)),'LineWidth',2)
% %         semilogx(f/1e9,10*log10(AMy./AMy(1)),'b','LineWidth',2)
% %         semilogx(f/1e9,10*log10(AMypar./AMypar(1)),'LineWidth',2)
%         semilogx(f/1e9,10*log10(AMytot./AMytot(1)),'LineWidth',2)
%         semilogx(f/1e9,10*log10(AMp./AMp(1)),'g--','LineWidth',2)
        
        p = p+1;
    end
    semilogx(mode_ss(end).fvet/1e9,CutoffLine,'k--','linewidth',1.2)
    lgd2 = legend([num2str(Ivet(1),' %0.1f '),' (exp)'],[num2str(Ivet(2),' %0.1f '),' (exp)'],...
        [num2str(Iplot(1),' %0.1f '),' (1D)'],[num2str(Iplot(2),' %0.1f '),' (1D)'],...
        [num2str(Iplot(3),' %0.1f '),' (1D)'],[num2str(Iplot(4),' %0.1f '),' (1D)'],...
        [num2str(Iplot(5),' %0.1f '),' (1D)'],[num2str(Iplot(6),' %0.1f '),' (1D)']);
    
    xlabel('Frequency (GHz)')
    ylabel('Normalized AM response (dB)')
    format bank
%     title(lgd2,'Bias current (mA)')
    str = '-3 dB line';
    t = annotation('textbox','String',str,'EdgeColor','none','FitBoxToText','on');
    % xlim([0 1000])
    set(gca,'xscale','log')
    xlim([1e-2 20]),ylim([-50 10])
end

%% 1D LI experimental curve
eval(['load MarkusN_4_T20.mat'])
Rmeas=diff(Vmeas)./diff(Imeas/1000);
Rmeas=[Rmeas(1:end-2)];
Imeas_res=Imeas(1:end-3);
Imax=max(Imeas);
Imax=15;
Pmax=max(Lmeas)*1.2;

pa=1:5:length(Imeas);

iR = find(mode.Vbias >= 1.61,1);
R_static = diff(mode.Vbias(iR-1:end))./diff(mode.I_dd(iR-1:end));
if mode.SS == 1
    for i = iR:length(mode_ss)
        R_dynamic(i-iR+1) = mode_ss(i).R(1);
    end
end

figure(19)
% Resistence=diff(mode.Vbias)./diff(mode.I_dd);
%             Resistence=[Resistence(1),Resistence];
            % plot(mode.ii_dd*1000,Resistence,Imeasinterp,Resmeasinterp)
%            plot(Imeasinterp,Resmeasinterp,'r-','LineWidth',2)
             pu=1:4:length(Rmeas);
              semilogy(Imeas_res(pu),Rmeas(pu),'ro','Markersize',3)            
              hold on
%                  lres=1:length(Resistence);
            %semilogy(mode.ii_dd(lres)*1000,Resistence(lres),[colo(kpar),'.-'],'Markersize',10)            
%             semilogy(mode.ii_dd(lres)*1000,Resistence(lres),'r.-','linewidth',2)            
            ylim([30,500])
            grid on
            xlim([0.4 5])
  %            xlim([0.5,1000*mode.ii_dd(end)+.1])
            xlabel('Current, mA')
            ylabel('Differential resistence, \Omega')
            hold on, grid on, box on
% plot(mode.I_dd(iR:end)*mode.Area*1e3,R_static/mode.Area,'ro-','linewidth',2)
% plot(mode.I_dd(iR:end)*mode.Area*1e3,R_dynamic/mode.Area,'b+','linewidth',1.5)
semilogy(mode.I_dd(iR:end)*mode.Area*1e3,R_static/mode.Area,'k--','linewidth',2)
if mode.SS == 1
    semilogy(mode.I_dd(iR:end)*mode.Area*1e3,R_dynamic/mode.Area,'g+','linewidth',1.5)
end
% xlabel('Current, mA'),ylabel('Differential Resistance, \Omega\cdotcm^{-2}')
% set(gca,'yscale','log')
legend('static','dynamic - 1Hz')  
            
figure(22)
hold on, grid on, box on
plot(mode.I_dd*mode.Area*1e3,'bo','linewidth',1)
ylabel('Current, mA')
            
figure(25)
hold on, grid on, box on
plot(mode.Pst,'bo','linewidth',1)
ylabel('Optical power, mW')

figure(23)
hold on, grid on, box on
plot(Vmeas(pa),Imeas(pa),['r.'],'markersize',10)
plot(mode.Vbias,mode.I_dd*mode.Area*1e3,'b','linewidth',2)
xlim([1.4 2])
xlabel('Bias voltage, V')
ylabel('Current, mA')

figure(24)
hold on, grid on, box on
plot(Imeas(pa),Lmeas(pa),'ro','markersize',3,'linewidth',2)
plot(mode.I_dd*mode.Area*1e3,mode.Pst,'b','linewidth',2)
xlim([0 2.1])
xlabel('Current, mA')
ylabel('Optical power, mW')
legend('Experimental','1D simulation','location','best')
