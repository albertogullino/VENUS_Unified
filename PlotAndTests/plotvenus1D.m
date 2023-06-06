% puntatore=mesh.puntatore;
puntatore=2:mesh.nny-1;

if mode.oflg==0
    figure
    hold on,grid on,box on
    plot(mesh.node(2,puntatore),mode.elec(puntatore)*mode.CarrierNorm,'.-')
    plot(mesh.node(2,puntatore),mode.hole(puntatore)*mode.CarrierNorm,'.-')
    legend('n','p')
else
    figure
    hold on,grid on,box on
    plot(mesh.node(2,puntatore),mode.elec(puntatore)*mode.CarrierNorm+mode.N2D(puntatore)*mode.CarrierNorm2D,'.-')
    plot(mesh.node(2,puntatore),mode.hole(puntatore)*mode.CarrierNorm+mode.P2D(puntatore)*mode.CarrierNorm2D,'.-')
end


figure,hold on
plot(mesh.node(2,puntatore)*1e4,mode.ecb(puntatore),'.-','linewidth',2)
plot(mesh.node(2,puntatore)*1e4,mode.evb(puntatore),'.-','linewidth',2)
plot(mesh.node(2,puntatore)*1e4,mode.EFn(puntatore),'k-.')
plot(mesh.node(2,puntatore)*1e4,mode.EFp(puntatore),'k--')

figure(80),hold on,plot(mode.vv_dd,mode.ii_dd*1e3*mode.CarrierNorm)

% QW population
% for ii=1:length(mode.vv_dd)
%     n2D(ii)=mode.nQW{ii}{1}(1);
%     p2D(ii)=mode.pQW{ii}{1}(1);
% end
% 
% figure
% grid on,hold on,box on
% plot(mode.ii_dd*1e3*mode.CarrierNorm,n2D*mode.CarrierNorm/mesh.vWMQW{1},'.-')
% plot(mode.ii_dd*1e3*mode.CarrierNorm,p2D*mode.CarrierNorm/mesh.vWMQW{1},'.-')
% chold
% plot(mode.ii_dd*1e3,mode.n3MaxVet,'--')
% plot(mode.ii_dd*1e3,mode.p3MaxVet,'--')
% xlabel('Current, mA'),ylabel('nQW, pQW, cm^{-3}')
% legend('n2D','p2D','n3D','p3D','location','northwest')
% set(gca,'FontSize',12)

% Lmod - Gmod
figure
hold on,grid on
plot(mode.ii_dd*1e3,mode.Gmod)
hold on,plot(mode.ii_dd*1e3,mode.Lmod)
% legend('Gmod 1','Gmod 2','Lmod 1','Lmod 2')
set(gca,'FontSize',12)
xlabel('Current, mA')

%%
iV=25;
iV=60;
% Bands
figure
grid on,hold on,box on
plot(mesh.node(2,1:mesh.nny)*1e4,modes.Ec(iV,:),'.-','linewidth',2)
plot(mesh.node(2,1:mesh.nny)*1e4,modes.Ev(iV,:),'.-','linewidth',2)
plot(mesh.node(2,1:mesh.nny)*1e4,modes.EFc(iV,:),'k-.')
plot(mesh.node(2,1:mesh.nny)*1e4,modes.EFv(iV,:),'k--')
% xlim([110 118])
title(['I = ',num2str(mode.ii_dd(iV)*1e3),' mA'])

% Current densities
node='off'; triangle='off'; color='on'; vpath='off'; arrows='off';
sd=1:geom.nd; scale=[]; cmap=[]; grido='off'; cbar='off';

figure
plot_tri_mesh(geom,mesh,[],sd,color,triangle,node,grido,cbar,cmap,scale,vpath,arrows)
xlabel('z, \mum')
ylabel('\rho, \mum')
axis normal
 
Jnx=squeeze(modes.JXn(iV,:,:));
Jny=squeeze(modes.JYn(iV,:,:));
Jpx=squeeze(modes.JXp(iV,:,:));
Jpy=squeeze(modes.JYp(iV,:,:));

 hold on,quiver(mesh.xgrid*1e4,mesh.ygrid*1e4,Jnx,Jny,.5)
 hold on,quiver(mesh.xgrid*1e4,mesh.ygrid*1e4,Jpx,Jpy,.5)
title(['I = ',num2str(mode.ii_dd(iV)*1e3),' mA'])

% Ccap
iV=25;
iV=27;
iV=18;

figure
grid on,hold on,box on
title(['I = ',num2str(mode.ii_dd(iV)*1e3),' mA'])
plot(mesh.xgrid(1:mesh.nnxQW{1}(1))*1e4,mode.Ccapn{2,iV},'.-')
plot(mesh.xgrid(1:mesh.nnxQW{1}(1))*1e4,mode.Ccapp{2,iV},'.-')
xlabel('\rho, \mum')
legend('Ccapn','Ccapp')

Ccapn=cell2mat(mode.Ccapn);
Ccapp=cell2mat(mode.Ccapp);

figure
grid on,hold on,box on
plot(mode.ii_dd*1e3,Ccapn(2,1:mesh.nnx:end),'.-')
plot(mode.ii_dd*1e3,Ccapp(2,1:mesh.nnx:end),'.-')
xlabel('Current, mA')
legend('Ccapn','Ccapp')

EFn2D=cell2mat(mode.EFn2D);
EFp2D=cell2mat(mode.EFp2D);

figure
grid on,hold on,box on
% title(['I = ',num2str(mode.ii_dd(iV)*1e3),' mA'])
plot(mode.ii_dd*1e3,EFn2D(2,1:mesh.nnx:end),'.-')
% plot(mode.ii_dd*1e3,mode.EFp2D{2,:}(1),'.-')
plot(mode.ii_dd*1e3,modes.EFc(:,mesh.inMQW{2}(1)),'.-')
% plot(mode.ii_dd*1e3,modes.EFv{2,iV},'.-')
xlabel('Current, mA')
legend('EFn2D','EFn3D')

figure
grid on,hold on,box on
plot(mode.ii_dd*1e3,(EFn2D(2,1:mesh.nnx:end)'-modes.EFc(:,mesh.inMQW{2}(1)))./mode.Vt(:,mesh.inMQW{2}(1)),'.-')
xlabel('Current, mA')