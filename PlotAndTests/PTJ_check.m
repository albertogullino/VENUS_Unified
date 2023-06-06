iRagTJ=mesh.nnxQW{1};
iRagTJ=18;

rho=mesh.xgrid;
drho=diff([0 rho]);
xdxN=rho.*drho;

Area=pi*rho(iRagTJ)^2;
I=mode.ii_dd*1e3;

VTJ=squeeze(mode.VTJ(1,:,1:iRagTJ));
Vrho=VTJ(:,1:iRagTJ)*xdxN(1:iRagTJ)'*2*pi/Area;


figure
hold on,grid on,box on
plot(rho*1e4,squeeze(mode.VTJ(1,:,:)))
plot(rho(iRagTJ)*1e4,squeeze(mode.VTJ(1,:,iRagTJ)),'k*')
xlabel('\rho, \mum'),ylabel('V_{TJ}, V')
% xlim([rho(1) rho(iRagTJ)*1e4])

figure
hold on,grid on,box on
plot(I,squeeze(mode.VTJ(1,:,:)))
plot(I,Vrho,'k','linewidth',2)
xlabel('Current, mA'),ylabel('V_{TJ}, V')