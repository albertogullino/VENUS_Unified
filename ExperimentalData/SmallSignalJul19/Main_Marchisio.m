% Static + S21 and S11
clear
close all
clc

wafer = '01103_448C';
device = 'C134R233';    % available: C134R233, C134R235, C134R237

z0 = 50;

Tvet = [20 50 80 110];         % STATIC available: 20, 50, 80, 110°C
% Small-signal: 20 (S21 published in 2020Gullino_NUSOD), 50°C
Ivet = [1 2 6];         % available: 1, 2, 3, 4, 5, 6

indLeg=0;
indLegS=0;

figure(1)
set(gcf,'Position',[150 150 1600 688])

figure(2)
set(gcf,'Position',[577 121 852 688])
hold on
box on

figure(3)
set(gcf,'Position',[577 121 852 688])
hold on
box on

figure(4)
set(gcf,'Position',[150 150 1600 688])
hold on
box on

for indT = 1:length(Tvet)
    %% Static characteristics
    load(['Static/Marchisio_d4_T',num2str(Tvet(indT)),'.mat'])
    figure(4)
    hold on,box on,grid on
    plot(Imeas,Lmeas,'.-','linewidth',2)
    xlabel('Current, mA'), ylabel('Optical power, mW')
    indLegS = indLegS+1;
    legStrS{indLegS} = ['T = ',num2str(Tvet(indT)),' °C'];

    if Tvet(indT)<60
        for indI = 1:length(Ivet)
            fileName = [wafer,'_',num2str(Tvet(indT)),'\',wafer,'_',device,'_I_',num2str(Ivet(indI)),'mA_',num2str(Tvet(indT)),'°C.s2p'];
            Sstr = sparameters(fileName);
            vFreq = Sstr.Frequencies;

            %% S21
            S21 = squeeze(Sstr.Parameters(2,1,:));

            figure(1)
            subplot(121)
            semilogx(vFreq/1e9,10*log10(abs(S21/S21(1)).^2),'LineWidth',2)
            grid on
            hold on
            box on

            xlabel('Frequency, GHz')
            ylabel('|S_{21}|^2')
            indLeg = indLeg+1;
            legStr{indLeg} = ['T = ',num2str(Tvet(indT)),' °C, I = ',num2str(Ivet(indI)),' mA'];

            %% S11
            S11 = squeeze(Sstr.Parameters(1,1,:));
            figure(1)
            subplot(122)
            %pausak
            %        plot(vFreq/1e9,10*log10(abs(S11).^2),'LineWidth',2)
            P=polar(angle(S11),abs(S11));
            set(P,'Linewidth',2)
            hold on
            %ylabel('S_{11}')
            %         legStr1{indLeg} = ['T = ',num2str(Tvet(indT)),' °C, I = ',num2str(Ivet(indI)),' mA'];
            %pausak

            %% Impedance (ABBIAMO PROVATO, NON SIAMO SICURI)
            figure(2)
            %         subplot(1,2,1)

            z=(1+S11)./(1-S11);
            z=z*z0;
            yyaxis left
            semilogx(vFreq/1e9,real(z))
            xlim([0.21 20])
            xlabel('Frequency, GHz')
            ylabel('Re\{Z\}, \Omega')
            yyaxis right
            %         plot(vFreq/1e9,imag(z))
            C=-1./(2*pi*vFreq.*imag(z));
            semilogx(vFreq/1e9,C)
            hold on
            box on
            xlim([0.21 20])
            ylabel('Im\{Z\}, \Omega')

            legStr1{indLeg} = ['T = ',num2str(Tvet(indT)),' °C, I = ',num2str(Ivet(indI)),' mA'];

            figure(3)
            %         subplot(1,2,2)
            hold on
            grid on
            box on
            xlim([0 20])
            plot(vFreq/1e9,abs(z),'linewidth',2)
            xlabel('Frequency, GHz'),ylabel('|Z|, \Omega')

        end
    end
end

figure(1)
subplot(121)
xlim([0,max(vFreq)/1e9])
legend(legStr,'Location','Best')
subplot(122)
title('S_{11}')
legend(legStr1,'Location','Best')

figure(2)
title('VCSEL impedance')
% subplot(121)
% legend(legStr1,'Location','North')
% subplot(122)
legend(legStr1,'Location','Best')

figure(3)
title('VCSEL impedance')
% subplot(121)
% legend(legStr1,'Location','North')
% subplot(122)
legend(legStr1,'Location','Best')

figure(4)
title('Static characteristics')
legend(legStrS,'Location','Best')
