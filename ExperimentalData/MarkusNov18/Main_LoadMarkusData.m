
clear
% close all
clc

Dact = 4; % available 2, 3, 4, 5, 6, 7, 8

% Different currents are available for each temperature.
cser='rgbmc';
Temperature = 20; % available 20, 50, 80, 110
Currents = {' 2.000',' 3.000',' 6.000',' 9.000','12.000'};
%Currents = {' 2.000',' 3.000',' 6.000',' 9.000','12.000','14.000','20.000'};
%Currents = {'2.000','3.000','5.000','8.000','10.000','14.000'};
% currents, mA, available 0.500, 1.000, 1.500, 2.000, 2.500, 3.000, 3.500, 4.000.
                      % Don't omit the .000 !!!!! 

deviceDir=[pwd,'\Dact_',num2str(Dact),'um'];

% DO NOT TOUCH: associating Dact to a device ID
switch Dact
    case 2
        deviceID='130229';
    case 3
        deviceID='130231';
    case 4
        deviceID='130233';
    case 5
        deviceID='130235';
    case 6
        deviceID='130237';
    case 7
        deviceID='130239';
    case 8
        deviceID='130241';
end
strDact=[num2str(Dact),'.0'];
LIVfileName=[deviceDir,'\LIV_',strDact,'um_Temp',num2str(Temperature),'_',deviceID,'.txt'];

[I,V,P] = f_ImportMarkusLIV(LIVfileName);

figure(1)
set(gcf,'Position',[295 518 1173 429])
subplot(1,2,1)
hold on
grid on
box on
plot(V,I)
xlabel('Voltage, V')
ylabel('Current, mA')
subplot(1,2,2)
hold on
grid on
box on
plot(I,P)
xlabel('Current, mA')
ylabel('Optical power, mW')

keyboard
figure(2)
set(gcf,'Position',[78 54 560 420])
hold on
grid on
box on
xlabel('Wavelength, nm')
ylabel('Intensity, A.U.')

for indCurrent=1:length(Currents)
    Spectrum_label{indCurrent}=[Currents{indCurrent},'mA'];
    SpectrumfileName=[deviceDir,'\OSI_',strDact,'um_Temp',num2str(Temperature),'_',Currents{indCurrent},'mA_',deviceID,'.txt'];
    [lambdavet,spectrum] = f_ImportMarkusSpectrum(SpectrumfileName);
    plot(lambdavet,spectrum)
    title(Spectrum_label{indCurrent})
    pausak
end

%legend(Spectrum_label,'Location','Best')

PeaksfileName=[deviceDir,'\OSI_Peaks_',strDact,'um_Temp',num2str(Temperature),'_',deviceID,'.txt'];
[Currentpeaks,lambdapeaks] = f_ImportMarkusSpectrumPeaks(PeaksfileName);

figure(3)
set(gcf,'Position',[1012 69 560 420])
hold on
grid on
box on
plot(Currentpeaks,lambdapeaks,'o')
xlabel('Current, mA')
ylabel('Wavelength, nm')