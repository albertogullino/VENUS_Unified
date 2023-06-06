% clear
%
% load('C:\Users\albig\Politecnico\Dottorato\3b_VENUS\20230213_VENUS_TJlitho\out\LW_MarkusN_FINALE_LgDBR_Thesis.mat')

if exist('modePlot')==0
    if length(MODEplot)>1
    iT=input('iT?\n');
    modep=MODEplot{iT};
else
    modep=MODEplot{1};
end
    % % xQW=P.mesh.xgrid(1:P.mesh.nnxQW{1})*1e4;
    xQW=mesh.xgrid(1:mesh.nnxQW{1})*1e4;
end
if exist('CORRENTI')==0
    CORRENTI=[1 4 8];
end
CorLav=CORRENTI;
corrente=modePlot.ii_dd*1000;
pu=[];
for kC=1:length(CorLav)
    [~,puk]=min(abs(corrente-CorLav(kC)));
    pu=[pu puk];
end

% pu=1:45;

s_LoadConstants
vph=Clight/mode.nindexQW;

Gamma_zMean=0.0099;                 % matgain è salvata solo su una QW ?
indQW=1;                            % qui si può mettere 1, tanto è uguale

nnQW=mesh.nnxQW{1};                 % number of lateral points in the QW (x-directed)
inQW=mesh.inMQW{indQW};
xQW=mesh.node(1,inQW);
LeQW=xQW(2:end)-xQW(1:end-1);       % edge length, cm
xcQW=(xQW(2:end)+xQW(1:end-1))/2;   % edge centers, cm
iiQW1=1:nnQW-1; iiQW2=2:nnQW;
Lp=zeros(1,nnQW);                   % box length, cm
Lp(1:(nnQW-1))=LeQW/2; Lp(2:nnQW)=Lp(2:nnQW)+LeQW/2;
Lp1=Lp(iiQW1)/2; Lp1(1)=2*Lp1(1);
Lp2=Lp(iiQW2)/2; Lp2(end)=2*Lp2(end);
% Cylindrical corrections
Lp1=2.*pi.*Lp1.*xcQW;
Lp2=2.*pi.*Lp2.*xcQW;

E=zeros(size(modePlot.E2(pu,:,1:nnQW)));
gE=E;

for iV=1:length(pu)
    E(iV,:,:)=modePlot.E2(pu(iV),:,1:nnQW);     % field from VELM
    E2=squeeze(E(iV,:,:));
    gE(iV,:,:)=modePlot.matgain(pu(iV),:).*E2.*Gamma_zMean;     % gain-field product
    
    GM(iV,:)=3*sum([Lp1.*squeeze(gE(iV,:,iiQW1)) Lp2.*squeeze(gE(iV,:,iiQW2))],2) ; % moltiplicato *3 (QW)
end

figure,grid on,box on
chold,plot(mode.ii_dd(pu)*1e3,mode.Lmod(:,pu)','*')
chold,plot(mode.ii_dd(pu)*1e3,GM(pu,:),'o')
xlabel('Current, mA'),ylabel('Modal gain vs. losses, 1/cm')
set(gca,'FontSize',14,'FontName','Times new roman')
legend('LP_{01}','LP_{11}','LP_{02}','location','best')
title('Stars: LOSSES - Open circles: GAIN')

x=mesh.xgrid*1e4;