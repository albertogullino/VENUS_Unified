x=mesh.xgrid*1e4;
y=mesh.ygrid*1e4;

%
Pcor=input(' Corrente = ')
if length(Pcor)==0
 Pcor=1.5
end 

Cor=1e3*mode.ii_dd;
% Cor=1000*modeold.ii_dd;         % Computed current from DD
[~,iCurr]=min(abs(Pcor-Cor));     % find closest current index to Cor (DD)


Prag=input(' Radius = ')
if length(Prag)==0
 Prag=0
end

[~,iRag]=min(abs(Prag-x));     % find closest current index to Rag (x)

%
Pz=input(' Radius = ')
if length(Pz)==0
 Pz=0
end

[~,iz]=min(abs(Pz-y));     % find closest current index to z (y)

%% DeltaT

var=modep.Temp-mode.T0;

figure,surf(x,y,squeeze(var(iCurr,:,:)))
shading interp
view(2)
ylim([0 120])
colorbar
xlabel('\rho, \mum')
ylabel('z, \mum')
set(gca,'FontSize',16,'FontName','Times new roman')
% title('I_{bias}=1.5 mA')
% axis([0 15 100 120])
axis([0 10 110 120])
xline(0,'r--','linewidth',4)
xline(2,'g--','linewidth',4)
xline(5,'c--','linewidth',4)
yline(115,'k--','linewidth',4)
yline(118.5,'m--','linewidth',4)

figure,plot(y,squeeze(var(iCurr,:,1)),'r','linewidth',2)
hold on
xlabel('z, \mum')
ylabel('\DeltaT, K')
grid on
% set(gca,'yscale','log')
set(gca,'FontSize',16,'FontName','Times new roman')
plot(y,squeeze(var(iCurr,:,14)),'g','linewidth',2)
plot(y,squeeze(var(iCurr,:,25)),'c','linewidth',2)
legend('\rho=0','\rho=2 \mum','\rho=5 \mum','location','northwest')

figure,plot(x,squeeze(var(iCurr,118,:)),'k','linewidth',2)
xlabel('z, \mum')
ylabel('\DeltaT, K')
grid on
set(gca,'FontSize',16,'FontName','Times new roman')
hold on,plot(x,squeeze(var(iCurr,end,:)),'m','linewidth',2)
legend('z=115 \mum (QW)','z=118.5 \mum (top)')
legend('\rho=0','\rho=2 \mum (OX)','\rho=5 \mum (contact)')
legend('\rho=0','\rho=2 \mum','\rho=5 \mum')

%% HeatSources
Joule=squeeze(modep.HeatJoule(iCurr,:,:));
FCA=squeeze(modep.HeatOptAbs(iCurr,:,:));
GR=squeeze(modep.HeatRec_13(iCurr,:,:)+modep.HeatRec_RAD(iCurr,:,:));
% var=modep.HeatRec_RAD;
Ccap=squeeze(modep.HeatRec_Cap(iCurr,:,:));

figure
surf(x,y,real(log10(Joule))),title('Joule')
% surf(x,y,real(log10(Ccap))),title('Ccap')
% surf(x,y,real(log10(GR))),title('GR')
% surf(x,y,real(log10(FCA))),title('FCA')
shading interp
view(2)
colorbar
xlabel('\rho, \mum')
ylabel('z, \mum')
set(gca,'FontSize',16,'FontName','Times new roman')
% title('I_{bias}=1.5 mA')
% axis([0 15 100 120])
axis([0 10 112 119])
xline(0,'r--','linewidth',4)
xline(2,'g--','linewidth',4)
xline(5,'c--','linewidth',4)
yline(115,'k--','linewidth',4)
yline(118.5,'m--','linewidth',4)

Prag=input(' Radius = ')
if length(Prag)==0
 Prag=0
end
[~,iRag]=min(abs(Prag-x));

figure
hold on,grid on
xlabel('z, \mum')
ylabel('Heat sources, W/cm^3')
plot(y,Joule(:,iRag),'linewidth',2)
plot(y,FCA(:,iRag),'linewidth',2)
plot(y,GR(:,iRag),'linewidth',2)
plot(y,Ccap(:,iRag),'linewidth',2)
% plot(y,FCA(:,iRag)+Joule(:,iRag)+GR(:,iRag)+Ccap(:,iRag),'k--','linewidth',1.5)

axis([114 117 1e-15 1e0])
set(gca,'yscale','log','FontSize',16,'FontName','Times new roman')
plot(y,squeeze(var(iCurr,:,14)),'g','linewidth',2)
plot(y,squeeze(var(iCurr,:,25)),'c','linewidth',2)
plot(y,squeeze(var(iCurr,:,25)),'c','linewidth',2)
legend('\rho=0','\rho=2 \mum','\rho=5 \mum')

figure,plot(x,squeeze(var(iCurr,118,:)),'k','linewidth',2)
xlabel('z, \mum')
ylabel('\DeltaT, K')
grid on
set(gca,'FontSize',16,'FontName','Times new roman')
hold on,plot(x,squeeze(var(iCurr,end,:)),'m','linewidth',2)
legend('z=115 \mum (QW)','z=118.5 \mum (top)')
legend('\rho=0','\rho=2 \mum (OX)','\rho=5 \mum (contact)')
legend('\rho=0','\rho=2 \mum','\rho=5 \mum')
