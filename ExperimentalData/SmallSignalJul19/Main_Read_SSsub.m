
clear
close all
clc

wafer = '01103_448C';
device = 'C134R233'; % available: C134R233, C134R235, C134R237
% Tvet = [20:10:60]; % available: 20, 30, 40, 50, 60
% Ivet = [1]; % available: 1, 2, 3, 4, 5, 6

Tvet = [20:10:60]; % available: 20, 30, 40, 50, 60
Ivet = [1:6]; % available: 1, 2, 3, 4, 5, 6

%% S11

% h=figure;
% set(h,'pos',[156          23        1656         960]);

figure(1)
set(gcf,'Position',[577 121 852 688])
hold on
box on

R=1;
C=1;

for indT = 1:length(Tvet)
    subplot(2,3,indT)
    indLeg=0;
    for indI = 1:length(Ivet)
        fileName = [wafer,'_',num2str(Tvet(indT)),'\',wafer,'_',device,'_I_',num2str(Ivet(indI)),'mA_',num2str(Tvet(indT)),'°C.s2p'];
        Sstr = sparameters(fileName);
        vFreq = Sstr.Frequencies;
        S11 = squeeze(Sstr.Parameters(R,C,:));
        
        plot(vFreq/1e9,10*log10(abs(S11/S11(1)).^2),'LineWidth',2)
        grid on
        grid minor
        hold on
        xlabel('Frequency, GHz')
        ylabel(['|S_{',num2str(R),num2str(C),'}|^2'])
        
        indLeg = indLeg+1;
        legStr{indLeg} = ['T = ',num2str(Tvet(indT)),' °C, I = ',num2str(Ivet(indI)),' mA'];
    end
    xlim([0,max(vFreq)/1e9])
    legend(legStr,'Location','Best')
    %    pausak
end


%% S21

figure(2)
set(gcf,'Position',[577 121 852 688])
hold on
box on

R=2;
C=1;

for indT = 1:length(Tvet)
    subplot(2,3,indT)
    indLeg=0;
    for indI = 1:length(Ivet)
        fileName = [wafer,'_',num2str(Tvet(indT)),'\',wafer,'_',device,'_I_',num2str(Ivet(indI)),'mA_',num2str(Tvet(indT)),'°C.s2p'];
        Sstr = sparameters(fileName);
        vFreq = Sstr.Frequencies;
        S21 = squeeze(Sstr.Parameters(R,C,:));
        
        plot(vFreq/1e9,10*log10(abs(S21/S21(1)).^2),'LineWidth',2)
        hold on
        grid on
        grid minor
        
        xlabel('Frequency, GHz')
        ylabel(['|S_{',num2str(R),num2str(C),'}|^2'])
        
        indLeg = indLeg+1;
        legStr{indLeg} = ['T = ',num2str(Tvet(indT)),' °C, I = ',num2str(Ivet(indI)),' mA'];
    end
    xlim([0,max(vFreq)/1e9])
    legend(legStr,'Location','Best')
    %    pausak
end



