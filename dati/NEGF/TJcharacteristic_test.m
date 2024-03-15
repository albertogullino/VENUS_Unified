clear

% load('out\LW_FBH-testTJ-Julian_forward_Full_GaAs_Nd=3e19.mat')
% load('out\LW_FBH-testTJ-Julian_forward_Full_AlGaAs_Nd=3e19.mat')

% load('out\LW_FBH-testTJ-Julian_forward_Inc_GaAs_Nd=3e19.mat')
load('out\LW_FBH-testTJ-Julian_GaAs_Nd=3e19_NEGFinc_DDinc.mat')
% mode.strName='Julian';
% mode.strName=[];

% mode.iLtunn=1 % Jfit=VTJ*J0*exp(-lmin./L0);

LUT_GBTBT   % where NEGF results are stored 

iT=1;       % temperature index in NEGF results

cT=mode.cBTBT(iT,:);
pDeg=size(mode.cBTBT,2)-1;

% Vint=0:0.01:1;
Vint=linspace(-0.5,1,100); % Negative = Forward bias in this code! (and in DD)
VV=Vint'.^(pDeg:-1:0);

% Extrapolation from NEGF results
if iLtunn==1
    I_interp=J0*Vint'.*exp(-VV*cT'./L0);
else
    I_interp=10.^(VV*cT')-mode.I0_NEGF(iT);
end


rOx=2.175 % um, JSTQE oxide aperture radius

Area=(rOx*1e-4)^2*pi;

colordef white

fprintf('TJ characteristics from NEGF\n')

figure(4)
set(gcf,'Position',[489   336   905   538])
subplot(121)
hold on,grid on,box on
plot(Vint,I_interp,'.-')%*Area*1e3,'.-')
xlabel('V_{TJ}, V'),ylabel('J_{TJ}, A/cm^2')
subplot(122)
semilogy(Vint,abs(I_interp),'.-')
hold on,grid on,box on
ylim([1e-5 1e5])
xlabel('V_{TJ}, V'),ylabel('J_{TJ}, A/cm^2')
