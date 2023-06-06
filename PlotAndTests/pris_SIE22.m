


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

% figure
% grid on,hold on,box on
% plot(mode.ii_dd*1e3,n2D*mode.CarrierNorm/mesh.vWMQW{1},'.-')
% plot(mode.ii_dd*1e3,p2D*mode.CarrierNorm/mesh.vWMQW{1},'.-')
% chold
% plot(mode.ii_dd*1e3,mode.n3MaxVet,'--')
% plot(mode.ii_dd*1e3,mode.p3MaxVet,'--')
% xlabel('Current, mA'),ylabel('nQW, pQW, cm^{-3}')
% legend('n2D','p2D','n3D','p3D','location','northwest')
% set(gca,'FontSize',12)

% Lmod - Gmod
% figure
% hold on,grid on
% plot(mode.ii_dd*1e3,mode.Gmod)
% hold on,plot(mode.ii_dd*1e3,mode.Lmod)
% set(gca,'FontSize',12)
% xlabel('Current, mA')

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

JN_X=squeeze(modep.JXn(iV,:,:));
JN_Y=squeeze(modep.JYn(iV,:,:));
JP_X=squeeze(modep.JXp(iV,:,:));
JP_Y=squeeze(modep.JYp(iV,:,:));

J_X=JN_X+JP_X;
J_Y=JN_Y+JP_Y;

figure(82)
hold on,grid on,
plot(mesh.xgrid*1e4,abs(J_Y(mesh.inMQW{2}(1),:)*mode.CarrierNorm))
xlabel('\rho, \mum')
ylabel('J_z, A/cm^2')

figure(83)
hold on,grid on,
plot(mesh.xgrid*1e4,abs(squeeze(modep.elec(iV,mesh.inMQW{2}(1),:))*mode.CarrierNorm))
plot(mesh.xgrid*1e4,abs(squeeze(modep.hole(iV,mesh.inMQW{2}(1),:))*mode.CarrierNorm))
xlabel('\rho, \mum')
ylabel('n_z, A/cm^2')

% Bands
figure
grid on,hold on,box on
plot(mesh.node(2,1:mesh.nny)*1e4,modep.Ec(iV,:),'.-','linewidth',2)
plot(mesh.node(2,1:mesh.nny)*1e4,modep.Ev(iV,:),'.-','linewidth',2)
plot(mesh.node(2,1:mesh.nny)*1e4,modep.EFc(iV,:),'k-.')
plot(mesh.node(2,1:mesh.nny)*1e4,modep.EFv(iV,:),'k--')
xlim([mesh.ygrid(mesh.ICAV(1))*1e4-0.5 mesh.ygrid(mesh.ICAV(1)+find(diff(mesh.ICAV)>1,1))*1e4+0.5])
title(['I = ',num2str(mode.ii_dd(iV)*1e3),' mA'])
xlabel('z, \mum'),ylabel('Energy, eV')

% Current densities
node='off'; triangle='off'; color='on'; vpath='off'; arrows='off';
sd=1:geom.nd; scale=[]; cmap=[]; grido='off'; cbar='off';

% figure
% plot_tri_mesh(geom,mesh,[],sd,color,triangle,node,grido,cbar,cmap,scale,vpath,arrows)
% xlabel('z, \mum')
% ylabel('\rho, \mum')
% axis normal
%  
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
plot(mesh.xgrid(1:mesh.nnxQW{1}(1))*1e4,mode.Ccapn{2,iV}*mode.CarrierNorm,'.-')
plot(mesh.xgrid(1:mesh.nnxQW{1}(1))*1e4,mode.Ccapp{2,iV}*mode.CarrierNorm,'.-')
xlabel('\rho, \mum'),ylabel('Capture rates, 1/(cm^3\cdots)')
legend('Ccapn','Ccapp')

% Ccapn=cell2mat(mode.Ccapn);
% Ccapp=cell2mat(mode.Ccapp);
% 
% figure
% grid on,hold on,box on
% plot(mode.ii_dd*1e3,Ccapn(2,1:mesh.nnxQW{1}:end),'.-')
% plot(mode.ii_dd*1e3,Ccapp(2,1:mesh.nnxQW{1}:end),'.-')
% xlabel('Current, mA'),ylabel('Capture rates, 1/(cm^3\cdots)')
% legend('Ccapn','Ccapp')


% Iniezione
Cn=modep.Cn;
Cn(:,1)=Cn(:,2);
Cp=modep.Cp;
Cp(:,1)=Cp(:,2);
Cpl=mode.ii_dd*1e3;

puc=iV;
pcor=puc;
SZ=.2;

quiverSTEP

SEGNO_elettroni=-1;  % -1 particelle
J_XN=SEGNO_elettroni*modep.JXn;
J_YN=SEGNO_elettroni*(modep.JYn);
J_XP=modep.JXp;
J_YP=(modep.JYp);

figure
set(gcf,'pos',[349   398   1525   525 ])
subplot(121)
plot(x,1e-15*Cn(pcor,:)*mode.CarrierNorm,'o-',x,1e-15*Cp(pcor,:)*mode.CarrierNorm,'.-')
xlim([0 x(mesh.nnxQW{1})])
legend('Electrons','Holes')
xlabel('\rho, \mum'),grid on
ylabel('2D Injection, 10^{15}/(cm^2\cdotns)')

subplot(122)
plot_tri_mesh(geom,mesh,[],sd,color,triangle,node,grido,cbar,cmap,scale,vpath,arrows)
yMQW=modep.yMQW{end}*1e4;
Xqw=[0 10 10 0];
thQ=((modep.yMQW{1}-modep.yMQW{end})+modep.vWMQW{1})*1e4;
Yqw=[yMQW yMQW  yMQW+thQ yMQW+thQ];
hhco=patch(Xqw,Yqw,'c'); set(hhco,'EdgeColor','none');
  
xlabel('z, \mum')
ylabel('\rho, \mum')
axis normal
hold on
% quiver(SX,SY,squeeze(J_XN(pcor,sY0,sX0)),squeeze(J_YN(pcor,sY0,sX0)),.2)
% quiver(SX,SY,squeeze(J_XP(pcor,sY0,sX0)),squeeze(J_YP(pcor,sY0,sX0)),.2)
quiver(mesh.xgrid*1e4,mesh.ygrid*1e4,squeeze(J_XN(pcor,:,:)),squeeze(J_YN(pcor,:,:)),.2)
quiver(mesh.xgrid*1e4,mesh.ygrid*1e4,squeeze(J_XP(pcor,:,:)),squeeze(J_YP(pcor,:,:)),.2)

Cur_2D=1;
xcm=x*1e-4; xQW=x(1:mesh.nnxQW{1}); WQW=mesh.vWMQW{1};
xDiv=(xcm+diff(xcm(1:2))/2)*2*pi;
jnQ_x=sum(squeeze(modep.JnQW(pcor,:,:)),2)*1000;
jpQ_x=sum(squeeze(modep.JpQW(pcor,:,:)),2)*1000;
jQn=SEGNO_elettroni*[jnQ_x' zeros(1,length(xcm)-length(xQW))]./xDiv;
jQp=[jpQ_x' zeros(1,length(xcm)-length(xQW))]./xDiv;

if Cur_2D==1
    yQWmedioE=modep.yMQW{3}*ones(size(sX0))*1e4;
    yQWmedioH=modep.yMQW{1}*ones(size(sX0))*1e4;
    quiver(x(sX0),yQWmedioE,jQn(sX0)*WQW,zeros(size(sX0)),.3,'b')
    quiver(x(sX0),yQWmedioH,jQp(sX0)*WQW,zeros(size(sX0)),.3,'r')
end
% ylim([mesh.ygrid(mesh.ICAV(1))*1e4-0.5 mesh.ygrid(mesh.ICAV(1)+find(diff(mesh.ICAV)>1,1))*1e4+0.5])

DensityCurrentPlot

prompt='Do you want detailed thermal results? (1: YES - Enter: NO) \n';
thPlot = str2num(input(prompt,'s'));

if thPlot==1
    HeatSourcePlot
end

% figure
% hold on,box on,grid on
% plot(mesh.ygrid*1e4,mode.Jn_y(1:mesh.nny),'.-')
% plot(mesh.ygrid*1e4,mode.Jp_y(1:mesh.nny),'.-')
% % xlim([114.4 115.6])
% ylabel('Current density, A/cm^2'),xlabel('z, \mum')
% legend('J_n','J_p','location','northwest')
% 
% figure
% hold on,box on,grid on
% plot(mesh.ygrid*1e4,mode.Jn_x(1:mesh.nny),'.-')
% plot(mesh.ygrid*1e4,mode.Jp_x(1:mesh.nny),'.-')
% % xlim([114.4 115.6])
% ylabel('Current density, A/cm^2'),xlabel('z, \mum')
% legend('J_n','J_p','location','northwest')