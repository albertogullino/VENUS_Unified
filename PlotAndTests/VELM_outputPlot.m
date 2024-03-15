I=mode.ii_dd*1e3;
Ivelm=I([VELMInfo.indVoltage]);
fpdif=[VELMInput.fPdif];

figure(999)
hold on,grid on,box on
plot(Ivelm,fpdif(1:3:end),'LineWidth',2)
xlabel('Current, mA'),ylabel('fPdif')