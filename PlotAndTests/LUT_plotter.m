function [xB,xQW,Vk,EgB,EgQW]=LUT_plotter(mode,geom,Imin,Imax)

s_LoadConstants
h2m0=hbar^2/(2*m0);

figure(80)
hold on,grid on,box on
plot(mode.vv_dd,mode.ii_dd*1e3,'.-','LineWidth',1.5)
xlabel('Voltage, V'),ylabel('Current, mA')
set(gca,'FontSize',14,'FontName','Arial','box','on')

% Extract the kink voltage
indI=find(mode.ii_dd*1e3>Imin & mode.ii_dd*1e3<Imax);
p=polyfit(mode.vv_dd(indI),mode.ii_dd(indI)*1e3,1);
Vk=-p(2)/p(1)

% LUT quantities are plotted
load(geom.GLUTm)
xB=meshQW.xmol_barrier;
xQW=meshQW.xmol_well;
mh=MefH';
SBVap=Ban.SBV(:,1)*ones(1,meshQW.num_kvectors)+1./mh*(h2m0.*Ban.kgrid.^2/qel);
VB1=Ban.SBV(1,1); CB1=Ban.SBC(1,1);
EgQW=CB1+VB1;    % subbands gap
EgB=meshQW.Eg+meshQW.C0-meshQW.V0;

% xmol (QW + barriers)
figure(111),hold on,grid on,box on
xlabel('z, nm'),ylabel('Molar fration')
set(gca,'FontSize',14,'FontName','Arial','box','on')
title(['W_{QW}=',num2str(meshQW.Lz*1e9),' nm'])
plot(meshQW.xc*1e9,meshQW.xmol,'linewidth',2)

% Subbands
figure
subplot(121)
plot(Ban.kgrid*1e-9,-Ban.SBV,'linewidth',1.5), chold,grid on
plot(Ban.kgrid*1e-9,-SBVap,'.','linewidth',2)
xlabel('k_{||}, nm^{-1}')
ylabel('Energy (VB), eV'),ylim([-0.2 0])
set(gca,'FontSize',14,'FontName','Arial','box','on')
subplot(122)
plot(Ban.kgrid*1e-9,Ban.SBC,'linewidth',1.5), chold,grid on
ylabel('Energy (CB), eV'),xlabel('k_{||}, nm^{-1}')
set(gca,'FontSize',14,'FontName','Arial','box','on')
title(['x_{B}=',num2str(meshQW.xmol_barrier),'-x_{QW}=',num2str(meshQW.xmol_well)])