clear

LUT_GBTBT   % where NEGF results are stored 


Vint=[0:0.01:1];

figure,hold on
plot(T_vet-273,'o')
ylabel('Temperature, °C')

Tin=input('\nWhich Temperature?     ') + 273;   % in K
[~,iT]=min(abs(Tin-T_vet));     % find closest current index to Cor (DD)
% iT=1;       % temperature index in NEGF results
plot(iT,Tin-273,'*')

cT=mode.cBTBT(iT,:);
pDeg=size(mode.cBTBT,2)-1;

% Extrapolation from NEGF results
for iVint=1:length(Vint)
    VV= Vint(iVint).^[pDeg:-1:0];
    I_interp(iVint)=10.^(VV*cT')-mode.I0_NEGF(iT);
end
%     II=1/LTJ*(-I_interp/qel);

rOx=2.175 % um, JSTQE oxide aperture radius

Area=(rOx*1e-4)^2*pi;

colordef white

fprintf('TJ characteristics from NEGF\n')

figure
set(gcf,'Position',[489   336   905   538])
subplot(121)
hold on,grid on,box on
plot(Vint,I_interp*Area*1e3,'.-')
xlabel('V_{TJ}, V'),ylabel('I_{TJ}, mA')
subplot(122)
semilogy(Vint,I_interp,'.-')
hold on,grid on,box on
xlabel('V_{TJ}, V'),ylabel('J_{TJ}, A/cm^2')
