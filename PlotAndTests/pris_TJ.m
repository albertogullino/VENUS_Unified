
JN_X=reshape(mode.Jn_x,mesh.nny,mesh.nnx);
JN_Y=reshape(mode.Jn_y,mesh.nny,mesh.nnx);
JP_X=reshape(mode.Jp_x,mesh.nny,mesh.nnx);
JP_Y=reshape(mode.Jp_y,mesh.nny,mesh.nnx);

J_X=JN_X+JP_X;
J_Y=JN_Y+JP_Y;

figure
surf(mesh.xgrid*1e4,mesh.ygrid*1e4,log10(abs(J_Y)))
title('J_z'),shading interp,xlabel('\rho'),ylabel('z')
view(2)
colorbar
caxis([0 4])

figure
hold on,grid on,
plot(mesh.xgrid*1e4,abs(J_Y(mesh.inMQW{1}(1),:)*mode.CarrierNorm))
xlabel('\rho, \mum')
ylabel('J_z, A/cm^2')

figure(80)
hold on
plot(mode.vv_dd,mode.ii_dd*1e3)
xlabel('Voltage, V')
ylabel('Current, mA')
box on,grid on

figure(81)
hold on
plot(abs(mode.ii_dd)*1e3,sum(mode.Pst_dd,1),'.-')
ylabel('Optical power, mW')
xlabel('Current, mA')
box on,grid on

clear n2D p2D
% QW population
for ii=1:length(mode.vv_dd)
    n2D(ii)=mode.nQW{ii}{1}(1);
    p2D(ii)=mode.pQW{ii}{1}(1);
end

figure
grid on,hold on,box on
plot(mode.ii_dd*1e3,n2D*mode.CarrierNorm/mesh.vWMQW{1},'.-')
plot(mode.ii_dd*1e3,p2D*mode.CarrierNorm/mesh.vWMQW{1},'.-')
chold
if length(mode.n3MaxVet)>length(mode.ii_dd)
    mode.n3MaxVet=mode.n3MaxVet(1:end-1);
    mode.p3MaxVet=mode.p3MaxVet(1:end-1);
end
plot(mode.ii_dd*1e3,mode.n3MaxVet,'--')
plot(mode.ii_dd*1e3,mode.p3MaxVet,'--')
xlabel('Current, mA'),ylabel('nQW, pQW, cm^{-3}')
legend('n2D','p2D','n3D','p3D','location','northwest')
set(gca,'FontSize',12)

% Lmod - Gmod
figure
hold on,grid on
plot(mode.ii_dd*1e3,mode.Gmod)
hold on,plot(mode.ii_dd*1e3,mode.Lmod)
set(gca,'FontSize',12)
xlabel('Current, mA')

%%
figure,plot(sum(mode.Pst_dd,1),'o-')
% keyboard
iV=25;
iV=input('iV?\n');
% iV=36;

if length(MODEplot)>1
    iT=input('iT?\n');
    modep=MODEplot{iT};
else
    modep=MODEplot{1};
end

% Bands
figure
grid on,hold on,box on
plot(mesh.node(2,1:mesh.nny)*1e4,modep.Ec(iV,:),'.-','linewidth',2)
plot(mesh.node(2,1:mesh.nny)*1e4,modep.Ev(iV,:),'.-','linewidth',2)
plot(mesh.node(2,1:mesh.nny)*1e4,modep.EFc(iV,:),'k-.')
plot(mesh.node(2,1:mesh.nny)*1e4,modep.EFv(iV,:),'k--')
xlim([mesh.ygrid(mesh.ICAV(1))*1e4 mesh.ygrid(mesh.ICAV(length(mesh.ICAV)/mesh.nnx))*1e4])
title(['I = ',num2str(mode.ii_dd(iV)*1e3),' mA'])

% Current densities
node='off'; triangle='off'; color='on'; vpath='off'; arrows='off';
sd=1:geom.nd; scale=[]; cmap=[]; grido='off'; cbar='off';

% figure
% plot_tri_mesh(geom,mesh,[],sd,color,triangle,node,grido,cbar,cmap,scale,vpath,arrows)
% xlabel('z, \mum')
% ylabel('\rho, \mum')
% axis normal
% %  
% Jnx=squeeze(modep.JXn(iV,:,:));
% Jny=squeeze(modep.JYn(iV,:,:));
% Jpx=squeeze(modep.JXp(iV,:,:));
% Jpy=squeeze(modep.JYp(iV,:,:));
% 
% hold on,quiver(mesh.xgrid*1e4,mesh.ygrid*1e4,Jnx,Jny,.5)
% hold on,quiver(mesh.xgrid*1e4,mesh.ygrid*1e4,Jpx,Jpy,.5)
% title(['I = ',num2str(mode.ii_dd(iV)*1e3),' mA'])

% Ccap
figure
grid on,hold on,box on
title(['I = ',num2str(mode.ii_dd(iV)*1e3),' mA'])
plot(mesh.xgrid(1:mesh.nnxQW{1}(1))*1e4,mode.Ccapn{2,iV},'.-')
plot(mesh.xgrid(1:mesh.nnxQW{1}(1))*1e4,mode.Ccapp{2,iV},'.-')
xlabel('\rho, \mum'),ylabel('Capture rates, 1/(cm^3\cdots)')
legend('Ccapn','Ccapp')

Ccapn=cell2mat(mode.Ccapn);
Ccapp=cell2mat(mode.Ccapp);

figure
grid on,hold on,box on
plot(mode.ii_dd*1e3,Ccapn(2,1:mesh.nnxQW{1}:end),'.-')
plot(mode.ii_dd*1e3,Ccapp(2,1:mesh.nnxQW{1}:end),'.-')
xlabel('Current, mA'),ylabel('Capture rates, 1/(cm^3\cdots)')
legend('Ccapn','Ccapp')

% Fermi levels
EFn2D=cell2mat(mode.EFn2D);
EFp2D=cell2mat(mode.EFp2D);

figure
set(gcf,'pos',[720   142   938   707])

subplot(221)
grid on,hold on,box on
plot(mode.ii_dd*1e3,EFn2D(2,1:mesh.nnxQW{1}:end),'.-')
plot(mode.ii_dd*1e3,modep.EFc(:,mesh.inMQW{2}(1)),'.-')
xlabel('Current, mA')
legend('EFn2D','EFn3D')

subplot(222)
grid on,hold on,box on
plot(mode.ii_dd*1e3,EFp2D(2,1:mesh.nnxQW{1}:end),'.-')
plot(mode.ii_dd*1e3,modep.EFv(:,mesh.inMQW{2}(1)),'.-')
xlabel('Current, mA')
legend('EFp2D','EFp3D')

subplot(223)
grid on,hold on,box on
plot(mode.ii_dd*1e3,(EFn2D(2,1:mesh.nnxQW{1}:end)'-modep.EFc(:,mesh.inMQW{2}(1)))./mode.Vt(:,mesh.inMQW{2}(1)),'.-')
xlabel('Current, mA'),ylabel('EFn difference')

subplot(224)
grid on,hold on,box on
plot(mode.ii_dd*1e3,(EFp2D(2,1:mesh.nnxQW{1}:end)'-modep.EFv(:,mesh.inMQW{2}(1)))./mode.Vt(:,mesh.inMQW{2}(1)),'.-')
xlabel('Current, mA'),ylabel('EFp difference')

fprintf('Quiver plots!\n'),pausak
quiverCorrenteNP_TJ

% Contain all the current densities at each bias point! (pcor=IV)
JY_cut=J_YN+J_YP;
JX_cut=J_XN+J_XP;

% indY=mesh.inMQW{1}(1);

figure(345)
hold on,box on,grid on
% plot(mesh.xgrid*1e4,abs(squeeze(JY_cut(pcor,indY,:))),'.-')
plot(mesh.xgrid*1e4,abs(squeeze(JY_cut(pcor,mesh.inMQW{3}(1),:))),'.-')
plot(mesh.xgrid*1e4,abs(squeeze(JY_cut(pcor,mesh.inMQW{2}(1),:))),'.-')
plot(mesh.xgrid*1e4,abs(squeeze(JY_cut(pcor,mesh.inMQW{1}(1),:))),'.-')
legend('QW 1','QW 2','QW 3','location','northeast')     
xlabel('\rho, \mum'),ylabel('J_y(\rho), A/cm^2')



HeatSourcePlot(mesh,mode,MODEplot)

