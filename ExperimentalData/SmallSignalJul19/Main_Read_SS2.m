% S12 and S22
clear
close all
clc

wafer = '01103_448C';
device = 'C134R233'; % available: C134R233, C134R235, C134R237
Tvet = [20:10:60]; % available: 20, 30, 40, 50, 60
Ivet = [5]; % available: 1, 2, 3, 4, 5, 6

% Tvet = [20]; % available: 20, 30, 40, 50, 60
% Ivet = [1:6]; % available: 1, 2, 3, 4, 5, 6

indLeg=0;
figure(1)
set(gcf,'Position',[100 301 1600 688])


%figure(2)
%set(gcf,'Position',[577 121 852 688])
%hold on
for indT = 1:length(Tvet)
    for indI = 1:length(Ivet)
        fileName = [wafer,'_',num2str(Tvet(indT)),'\',wafer,'_',device,'_I_',num2str(Ivet(indI)),'mA_',num2str(Tvet(indT)),'°C.s2p'];
        Sstr = sparameters(fileName);
        vFreq = Sstr.Frequencies;
%         S21 = squeeze(Sstr.Parameters(1,2,:));
        S12 = squeeze(Sstr.Parameters(1,2,:));
        figure(1)
        subplot(121)
        grid on
        hold on
        box on
%         plot(vFreq/1e9,10*log10(abs(S21/S21(1)).^2),'LineWidth',2)
        plot(vFreq/1e9,10*log10(abs(S12/S12(1)).^2),'LineWidth',2)
        xlabel('Frequency, GHz')
%         ylabel('|S_{21}|^2')
        ylabel('|S_{12}|^2')
        indLeg = indLeg+1;
        legStr{indLeg} = ['T = ',num2str(Tvet(indT)),' °C, I = ',num2str(Ivet(indI)),' mA'];        
        
%         S11 = squeeze(Sstr.Parameters(2,2,:));
        S22 = squeeze(Sstr.Parameters(2,2,:));
        figure(1)
        subplot(122)

        %pausak
%        plot(vFreq/1e9,10*log10(abs(S11).^2),'LineWidth',2)
%         P=polar(angle(S11),abs(S11));
        P=polar(angle(S22),abs(S22));
        set(P,'Linewidth',2)
        hold on
        %ylabel('S_{11}')
        legStr1{indLeg} = ['T = ',num2str(Tvet(indT)),' °C, I = ',num2str(Ivet(indI)),' mA'];        
        %pausak
    end
end

figure(1)
subplot(121)
xlim([0,max(vFreq)/1e9])
legend(legStr,'Location','Best')


figure(1)
subplot(122)
% title('S_{11}')
title('S_{22}')
%legend(legStr1,'Location','Best')