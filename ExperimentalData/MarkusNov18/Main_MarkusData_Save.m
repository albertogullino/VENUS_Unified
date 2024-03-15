
clear
close all
clc

Dact = input(' Device = (2 - 8 ) '); % available 2 - 8

Temperature = input(' Temp = (20, 50, 80, 110) '); % available 4, 5, 6
%Temperature = 20; % temperature, K, available 20, 30, 40, 50, 60, 70

%Dact = 4; % available 2, 3, 4, 5, 6, 7, 8

%Temperature = 20; % available 20, 50, 80, 110
Currents = {'2.000','3.000','4.000','5.000','6.000','7.000'};
% currents, mA, available 0.500, 1.000, 1.500, 2.000, 2.500, 3.000, 3.500, 4.000.
                      % Don't omit the .000 !!!!! 
cser='rgbymc';
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


PeaksfileName=[deviceDir,'\OSI_Peaks_',strDact,'um_Temp',num2str(Temperature),'_',deviceID,'.txt'];
[Currentpeaks,lambdapeaks] = f_ImportMarkusSpectrumPeaks(PeaksfileName);


figure(2)
set(gcf,'Position',[78 54 560 420])
hold on
grid on
box on
xlabel('Wavelength, nm')
ylabel('Intensity, A.U.')

for indCurrent=1:length(Currents)
    CurVal=str2num(Currents{indCurrent});
    Spectrum_label{indCurrent}=[Currents{indCurrent},'mA'];
    if str2num(Currents{indCurrent})<10
     SpectrumfileName=[deviceDir,'\OSI_',strDact,'um_Temp',num2str(Temperature),'_ ',Currents{indCurrent},'mA_',deviceID,'.txt'];
    else
     SpectrumfileName=[deviceDir,'\OSI_',strDact,'um_Temp',num2str(Temperature),'_',Currents{indCurrent},'mA_',deviceID,'.txt'];
    end 
    [lambdavet,spectrum] = f_ImportMarkusSpectrum(SpectrumfileName);
    
    iC=find(CurVal==Currentpeaks);
    LAV=lambdapeaks(iC,:);
    fi=find(isnan(LAV)==0);
    Po=NaN;
    Lo=NaN;
    if length(fi)>0
    for kk=fi
     li=LAV(kk);
     [du,fil]=min(abs(lambdavet-li));
     Po(kk)=spectrum(fil);
     Lo(kk)=li;
    end
    end
    x(indCurrent,1:length(Po))=Lo;
    p(indCurrent,1:length(Po))=Po;
    cu(indCurrent)=CurVal;
    plot(lambdavet,spectrum,cser(indCurrent)), hold on, plot(Lo,Po,['o',cser(indCurrent)])
      title(Spectrum_label{indCurrent})
    pausak
end

legend(Spectrum_label,'Location','Best')

x0=x;
p0=p;
fiV=find(isnan(x(:,1))==0);
x=x(fiV,:);
p=p(fiV,:);
cu=cu(fiV)';

si=size(x,2);
if si>1
 dx=diff(x(:,1))./diff(cu);
 fi=find(dx<.6*mean(dx));
 if length(fi)>0
  fTr=fi(1)+1:length(x(:,1));
  x(fTr,2:si)= x(fTr,1:si-1);
  x(fTr,1)= NaN;
  p(fTr,2:si)= p(fTr,1:si-1);
  p(fTr,1)= NaN; 
 end 
end 
fi=find(x==0);
x(fi)=NaN;
p(fi)=NaN;

figure(3)
set(gcf,'Position',[1012 69 560 420])
hold on
grid on
box on
plot(cu,x,'o')
xlabel('Current, mA')
ylabel('Wavelength, nm')
pausak



Vmeas=V;
Imeas=I;
Lmeas=P;

%x=ginput;

LAM=x;
P_dB=p;
Cur=cu;
P_lin=10.^(P_dB/10);
fiN=find(isnan(P_lin)==1);
P_lin(fiN)=0;
Ps=sum(P_lin,2);
for k=1:length(Ps)
 [du,fip]=min(abs(I-Cur(k)));
 Pl=P(fip);
 Cp(k)=Ps(k)/Pl;
end
P_ver=diag(1./Cp)*P_lin;

figure, plot(I,P,Cur,sum(P_ver,2),'o')
pausak
P_dB=10*log10(P_ver);
%save L1 Cur LAM P_dB
figure,
subplot(121), plot(Cur,P_ver)
subplot(122), plot(Cur,P_dB)
pausak


%'prima di salvare', keyboard
eval(['save Marchisio', num2str(Dact),'_T', num2str(Temperature),' Cur LAM P_dB Lmeas Vmeas Imeas'])