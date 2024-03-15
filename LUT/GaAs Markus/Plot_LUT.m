clear
close all

s_LoadConstants

% inT=input('Temperature?\n');
% inPor=input('Carrier levels?\n');

% load('LUT4D_Julian_nMark_xmol=0.06.mat')
load('LUT4D_Jun_Markus_nMark_40_more.mat')
load('LUT4D_Jun_Markus_nMark_40.mat')
load('LUT4D_Jun_Markus_nMark_40_Der.mat')

load('LUT4D_Feb24_Markus_nMark_more.mat')
load('LUT4D_Feb24_Markus_nMark.mat')
load('LUT4D_Feb24_Markus_nMark_Der.mat')

lav=lav*1e3;

% port=Densityv/(mesh.Lz*100)*1e-18;

% dim(G)=(por_E,por_H,lav,Tv)
figure(1)
hold on,grid on,box on
% Assume one population constant!
plot(lav,squeeze(squeeze(G(end,end,:,1))),'linewidth',2)

% load('LUT4D_Julian_nMark_xmol=0.06_T=350K.mat')
% Assume one population constant!
chold
plot(lav,squeeze(G(1,1,:,2)),'--','linewidth',2)

xlim([lav(1) lav(end)])
xlabel('Wavelength, nm')
ylabel('Gain, cm^{-1}')
set(gca,'FontSize',14,'FontName','Arial','Box','on')

% legend('Dens=0.1e12','Dens=1e12','Dens=2e12','Dens=3e12','Dens=4e12','numcolumns',1)

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

figure(3)
    set(gcf,'Position',[1060 512 560 420])
    hold on
    grid on
    plot(Ban.kgrid,(-Ban.SBV),'LineWidth',2),chold
    plot(Ban.kgrid,-SBVap,'.','linewidth',2), ylim([-0.2 0])
    set(gca,'FontSize',14,'FontName','Arial','Box','on')
    xlabel('k_{||}, nm^{-1}')
    ylabel('Energy, eV')
    title('Valence subbands')
    
    figure(4)
    set(gcf,'Position',[633 59 560 420])
    hold on
    grid on
    plot(Ban.kgrid,(Ban.SBC),'LineWidth',2)
    set(gca,'FontSize',14,'FontName','Arial','Box','on')
    xlabel('k_{||}, nm^{-1}')
    ylabel('Energy, eV')
    title('Conduction subbands')

