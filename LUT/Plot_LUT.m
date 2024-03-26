clear
% close all
%
s_LoadConstants


iLUT=1;
iLUT=input('LUT index?\n')

LUTname{1}='LUT4D_Jun_Markus_nMark_40';
LUTname{2}='LUT4D_Feb24_Markus_nMark';

LUTname{3}='LUT4D_Markus_nMark_xW=0_xB=0.2';
LUTname{4}='LUT4D_Markus_nMark_xW=0_xB=0.39';

LUTname{5}='LUT4D_Feb24_Markus_nMark_xW=0_xB=0.2';
LUTname{6}='LUT4D_Feb24_Markus_nMark_xW=0_xB=0.286';
LUTname{7}='LUT4D_Feb24_Markus_nMark_xW=0_xB=0.38';

LUTname{8}='LUT4D_Mar24_Julian795_nMark_xW=0.06_xB=0.33';

LUTname{iLUT}

LUTbands=[LUTname{iLUT},'_more.mat']; % results from Schrodinger equation
LUTgain=[LUTname{iLUT},'.mat'];       % the core of the LUT, with g, rsp, T, lamdba, carriers
LUTder=[LUTname{iLUT},'_Der.mat'];    % to compute the Jacobian matrix

load(LUTbands)
load(LUTgain)

lav=lav*1e3;

figure(1)
grid on,box on,hold on
grid minor
plot(Tv,'o')
title('Temperature')

inT=input('Temperature?\n');
plot(inT,Tv(inT),'k*')
pausak

figure(2)
grid on,box on,hold on
plot(port_2D,'o')
title('Carriers')

inPor=input('Carriers?\n');
plot(inPor,port_2D(inPor),'k*')
pausak
% port=Densityv/(mesh.Lz*100)*1e-18;

% dim(G)=(por_E,por_H,lav,Tv)
figure(10)
hold on,grid on,box on
% Assume one population constant!
plot(lav,squeeze(squeeze(G(inPor,inPor,:,inT))),'linewidth',2)

% load('LUT4D_Julian_nMark_xmol=0.06_T=350K.mat')
% Assume one population constant!
% chold
% plot(lav,squeeze(G(1,1,:,2)),'--','linewidth',2)

xlim([lav(1) lav(end)])
xlabel('Wavelength, nm')
ylabel('Gain, cm^{-1}')
set(gca,'FontSize',14,'FontName','Arial','Box','on')

% legend('Dens=0.1e12','Dens=1e12','Dens=2e12','Dens=3e12','Dens=4e12','numcolumns',1)
pausak

mesh=meshQW;
mode=modeQW;

port_2De=diag(eD1);
port_2Dh=diag(hD1);

h2m0=hbar^2/(2*m0);


mh=MefH';
SBVap=Ban.SBV(:,1)*ones(1,mesh.num_kvectors)+1./mh*(h2m0.*Ban.kgrid.^2/qel);
% figure
% grid on,box on,hold on
% plot(Ban.kgrid,Ban.SBV,'linewidth',1.5), chold
% plot(Ban.kgrid,SBVap,'.','linewidth',2), ylim([0 .2])
% chold
% plot(Ban.kgrid,-Ban.SBC,'linewidth',1.5),
pausak

figure(11)
set(gcf,'Position',[1060 512 560 420])
hold on
grid on
plot(Ban.kgrid,(-Ban.SBV),'LineWidth',2),chold
plot(Ban.kgrid,-SBVap,'.','linewidth',2), ylim([-0.2 0])
set(gca,'FontSize',14,'FontName','Arial','Box','on')
xlabel('k_{||}, nm^{-1}')
ylabel('Energy, eV')
title('Valence subbands')

figure(12)
set(gcf,'Position',[633 59 560 420])
hold on
grid on
plot(Ban.kgrid,(Ban.SBC),'LineWidth',2)
set(gca,'FontSize',14,'FontName','Arial','Box','on')
xlabel('k_{||}, nm^{-1}')
ylabel('Energy, eV')
title('Conduction subbands')

