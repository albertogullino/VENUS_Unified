
clear
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is aimed at plotting data from Anders Larsson measurements,
% received on 2019-06-10 from Pierluigi Debernardi.
%
% By cross-checking the various reports it can be understood that 
% 
% - structure 3TJ6 is the p+ on n+ (Gullino is working on this at present)
% - structure 4TJ3 is the n+ on p+
%
% moreover, 0501 (11um), 0601 (13um), 0701 (15um), 0801 (17um)
%           0502 (10um), 0602 (12um), 0702 (14um), 0802 (16um)
% are codes of measurements performed for structures with different oxide
% radii; 0501-0801 refer to the p+ on n+ structure, while 0502-0802 refer
% on the n+ on p+ device
%
% The 5 um oxide aperture cases seem to be peculiar, so it is not reported
% here. Moreover, wider apertures are more interesting, in view of 
% performing comparisons with 1D analyses (less crowding effects at the 
% oxide
%
% Alberto Tibaldi, 2019/11/22

% strType = 'n_on_p'; % structure 4TJ3
strType = 'p_on_n'; % structure 3TJ6

% d_ox = [10 12 14 16]; % diameter of the oxide aperture, um
d_ox = [11 13 15 17]; % diameter of the oxide aperture, um

for indox = 1:length(d_ox)
    
    load([strType,'_dox=',num2str(d_ox(indox)),'.mat'])

    Area = pi*(d_ox(indox)/2*1e-6)^2; % oxide aperture area, m^2
    figure(1)
    set(gcf,'Position',[311 249 1667 699])

    subplot(1,2,1)
    hold on
    grid on
    box on
    plot(Vvet,Ivet,'LineWidth',1.5)
    xlabel('Voltage, V')
    ylabel('Current, A')

    subplot(1,2,2)
    hold on
    grid on
    box on
    plot(Vvet,gradient(Vvet)./gradient(Ivet),'LineWidth',1.5)
    xlabel('Voltage, V')
    ylabel('Differential resistance, \Omega')
        
    legStr{indox} = ['d_{ox} = ',num2str(d_ox(indox)),' um'];
    
end

subplot(1,2,1)
legend(legStr,'Location','Best')
set(gca,'FontSize',14,'FontName','Times New Roman')

subplot(1,2,2)
legend(legStr,'Location','Best')
set(gca,'FontSize',14,'FontName','Times New Roman')

