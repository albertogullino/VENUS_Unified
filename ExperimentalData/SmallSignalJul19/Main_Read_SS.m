
clear
close all
clc

wafer = '01103_448C';
device = 'C134R233'; % available: C134R233, C134R235, C134R237
Tvet = [20:10:60]; % available: 20, 30, 40, 50, 60
% Ivet = [3]; % available: 1, 2, 3, 4, 5, 6

%Tvet = [20]; % available: 20, 30, 40, 50, 60
Ivet = [1:6]; % available: 1, 2, 3, 4, 5, 6

indLeg=0;
figure
set(gcf,'Position',[577 121 852 688])
% grid on
% hold on
% box on

for indT = 1:length(Tvet)
    for indI = 1:length(Ivet)
        fileName = [wafer,'_',num2str(Tvet(indT)),'\',wafer,'_',device,'_I_',num2str(Ivet(indI)),'mA_',num2str(Tvet(indT)),'°C.s2p'];
        Sstr = sparameters(fileName);
        vFreq = Sstr.Frequencies;
        S21 = squeeze(Sstr.Parameters(2,1,:));

        semilogx(vFreq/1e9,10*log10(abs(S21/S21(1)).^2),'LineWidth',2)
%                 plot(vFreq/1e9,abs(S21),'LineWidth',2)  % original

        grid on
        hold on
        box on
        xlabel('Frequency, GHz')
        ylabel('|S_{21}|^2')
        
        indLeg = indLeg+1;
        legStr{indLeg} = ['T = ',num2str(Tvet(indT)),' °C, I = ',num2str(Ivet(indI)),' mA'];
    end
end

xlim([0,max(vFreq)/1e9])
legend(legStr,'Location','Best')