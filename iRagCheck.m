% iRagHeatTJA=18       % radial index where the HeatTJ and PTJ mean are taken

xPerc=40;

% Here we want to compute the integral rho drho and divide it by the Area
for ii=1:length(A.mode.ii_dd)
    [maxQ,iQ]=max(squeeze(A.mode.HeatTJ(1,ii,:)));
    avgQ=maxQ*xPerc/100;
    ifQA(ii)=find(squeeze(A.mode.HeatTJ(1,ii,iQ:end))<avgQ,1)+iQ;
end

for ii=1:length(B.mode.ii_dd)
    [maxQ,iQ]=max(squeeze(B.mode.HeatTJ(1,ii,:)));
    avgQ=maxQ*xPerc/100;
    ifQB(ii)=find(squeeze(B.mode.HeatTJ(1,ii,iQ:end))<avgQ,1)+iQ;
end

%%
% Here we want to compute the integral rho drho and divide it by the Area
[maxA,iA]=max(squeeze(A.mode.HeatTJ(1,end,:)));
avgA=maxA*xPerc/100;
ifA=find(squeeze(A.mode.HeatTJ(1,end,iA:end))<avgA,1)+iA;
iRagHeatTJA=ifA

rhoA=A.mesh.xgrid;
drhoA=diff([0 rhoA]);
xdxNA=rhoA.*drhoA;

AreaA=pi*rhoA(iRagHeatTJA)^2;

% Mean are computed
VTJA=squeeze(A.mode.VTJ(1,:,1:iRagHeatTJA));
VrhoA=VTJA*xdxNA(1:iRagHeatTJA)'*2*pi/AreaA;

%%
% Here we want to compute the integral rho drho and divide it by the Area
[maxB,iB]=max(squeeze(B.mode.HeatTJ(1,end,:)));
avgB=maxB*xPerc/100;
ifB=find(squeeze(B.mode.HeatTJ(1,end,iB:end))<avgB,1)+iB;
iRagHeatTJB=ifB

rhoB=B.mesh.xgrid;
drhoB=diff([0 rhoB]);
xdxNB=rhoB.*drhoB;

AreaB=pi*rhoB(iRagHeatTJA)^2;

% Mean are computed
VTJB=squeeze(B.mode.VTJ(1,:,1:iRagHeatTJB));
VrhoB=VTJB*xdxNB(1:iRagHeatTJB)'*2*pi/AreaB;

xA=[1 cumsum(A.geom.div_x(1:end-1))+1];
figure,hold on,  grid on
% plot(A.mode.vv_dd(1:end-5),squeeze(A.mode.VTJ(1,1:end-5,xA)),'linewidth',2)
plot(A.mode.ii_dd*1e3,squeeze(A.mode.VTJ(1,:,xA)),'linewidth',2)
xlabel('Current, mA'),ylabel('V_{TJ}, V')
legendCell = cellstr([num2str(rhoA(xA)'*1e4, 'r=%.3f um')]);
legend(legendCell,'location','northwest')
set(gca,'FontSize',14,'FontName','Times new roman')
title('above')

xB=[1 cumsum(B.geom.div_x(1:end-1))+1];
figure,hold on, grid on
% chold,plot(B.mode.vv_dd,squeeze(B.mode.VTJ(1,:,xB)),'--','linewidth',2)
chold,plot(B.mode.ii_dd*1e3,squeeze(B.mode.VTJ(1,:,xB)),'--','linewidth',2)
xlabel('Current, mA'),ylabel('V_{TJ}, V')
legendCell = cellstr([num2str(rhoB(xB)'*1e4, 'r=%.3f um')]);
legend(legendCell,'location','northwest')
set(gca,'FontSize',14,'FontName','Times new roman')
title('below')

figure,hold on,box on,grid on
plot(A.mesh.xgrid*1e4,squeeze(A.mode.VTJ(1,iCurrA,:)),'linewidth',2)
chold
plot(B.mesh.xgrid*1e4,squeeze(B.mode.VTJ(1,iCurrB,:)),'--','linewidth',2)
xlim([0 13]),xlabel('\rho, \mum'),ylabel('V_{TJ}, V')
set(gca,'FontSize',14,'FontName','Times new roman')
legendCell = cellstr([num2str(round(B.mode.ii_dd(iCurrB)'*1e3), 'I=%d mA')]);
legend(legendCell,'location','northeast')

%%
% VTJ vs current
figure
hold on, grid on, box on
plot(B.mode.ii_dd*1e3,squeeze(B.mode.VTJ))
plot(B.mode.ii_dd*1e3,VrhoB,'k','linewidth',2)
xlabel('Current, mA'),ylabel('V_{TJ}')
title('below')

figure
hold on, grid on, box on
plot(A.mode.ii_dd*1e3,squeeze(A.mode.VTJ))
plot(A.mode.ii_dd*1e3,VrhoA,'k','linewidth',2)
xlabel('Current, mA'),ylabel('V_{TJ}')
title('above')

% VTJ vs xgrid
figure
hold on, grid on, box on
plot(B.mesh.xgrid*1e4,squeeze(B.mode.VTJ)')
plot(B.mesh.xgrid(iRagHeatTJB)*1e4,squeeze(B.mode.VTJ(1,:,iRagHeatTJB)),'ko')
xlim([0 13])
xlabel('\rho, \mum'),ylabel('V_{TJ}')
title('below')

figure
hold on, grid on, box on
plot(A.mesh.xgrid*1e4,squeeze(A.mode.VTJ))
plot(A.mesh.xgrid(iRagHeatTJA)*1e4,squeeze(A.mode.VTJ(1,:,iRagHeatTJA)),'ko')
xlim([0 13])
xlabel('\rho, \mum'),ylabel('V_{TJ}')
title('above')

% QTJ vs xgrid
figure
hold on, grid on, box on
plot(B.mesh.xgrid*1e4,squeeze(B.mode.HeatTJ)')
plot(B.mesh.xgrid(iRagHeatTJB)*1e4,squeeze(B.mode.HeatTJ(1,:,iRagHeatTJB)),'ko')
xlim([0 13])
xlabel('\rho, \mum'),ylabel('Q_{TJ}')
title('below')

figure
hold on, grid on, box on
plot(A.mesh.xgrid*1e4,squeeze(A.mode.HeatTJ))
plot(A.mesh.xgrid(iRagHeatTJA)*1e4,squeeze(A.mode.HeatTJ(1,:,iRagHeatTJA)),'ko')
xlim([0 13])
xlabel('\rho, \mum'),ylabel('Q_{TJ}')
title('above')




 
                    
                    