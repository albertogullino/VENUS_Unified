% clear
% close all
dbstop if error

if exist('CORRENTI')==0
%     CORRENTI=[2 2 3];
end

addpathVENUS

% P=load('out\LW_MarkusN_FINALE_LgDBR_NUSOD2023.mat');
% A=load('out\LW_MarkusN_TJ_oxAbove_LgDBR_FINAL.mat');
% B=load('out\LW_MarkusN_TJ_oxBelow_LgDBR_NUSOD2023.mat');

% P=load('out\LW_MarkusN_FINALE_LgDBR_20C_radial.mat');
% A=load('out\LW_MarkusN_TJ_oxAbove_LgDBR_20C_lambda_radial.mat');
% B=load('out\LW_MarkusN_TJ_oxBelow_LgDBR_20C_lambda_radial.mat');

% P=load('out\LW_MarkusN_FINALE_LgDBR_1D_nomireq.mat');
% A=load('out\LW_MarkusN_TJ_oxAbove_LgDBR_1D_nomireq.mat');
% B=load('out\LW_MarkusN_TJ_oxBelow_LgDBR_1D_nomireq.mat');


if length(P.MODEplot)>1
    iT=input('iT?\n');
    P.modeplot=P.MODEplot{iT};
    A.modeplot=A.MODEplot{iT};
    B.modeplot=B.MODEplot{iT};
else
    P.modeplot=P.MODEplot{end};
    A.modeplot=A.MODEplot{end};
    B.modeplot=B.MODEplot{end};
end

iCurrP=CurrIndex(CORRENTI,P.modeplot.ii_dd*1e3);
iCurrA=CurrIndex(CORRENTI,A.modeplot.ii_dd*1e3);
iCurrB=CurrIndex(CORRENTI,B.modeplot.ii_dd*1e3);

% indA=find([A.VELMInput.indVoltage]==iCurrA(2));
% indB=find([B.VELMInput.indVoltage]==iCurrB(2));
% indP=find([P.VELMInput.indVoltage]==iCurrP(2));

yA=A.mesh.ygrid*1e4;
yB=B.mesh.ygrid*1e4;
yP=P.mesh.ygrid*1e4;

elecA=squeeze(A.modeplot.elec(iCurrA(1),:,:))*1e-18;
elecB=squeeze(B.modeplot.elec(iCurrB(1),:,:))*1e-18;
elecP=squeeze(P.modeplot.elec(iCurrP(1),:,:))*1e-18;
holeA=squeeze(A.modeplot.hole(iCurrA(1),:,:))*1e-18;
holeB=squeeze(B.modeplot.hole(iCurrB(1),:,:))*1e-18;
holeP=squeeze(P.modeplot.hole(iCurrP(1),:,:))*1e-18;

if isfield(A.mesh,'IBTJ')==1 % expand for MTJ! (each one has a different sigma)
    if A.mode.flgBTJ_lithographic>0
        [~,iRagTJ]=min(abs(A.mesh.xgrid*1e4-A.mode.rAperture));
    else
        iRagTJ=A.mesh.nnxQW{1};
    end
    for iTJ=1:iRagTJ
        elecA(A.mesh.LBTJ(iTJ):A.mesh.RBTJ(iTJ))=elecA(A.mesh.LBTJ(iTJ):A.mesh.RBTJ(iTJ))+A.mode.dop_dp(A.mesh.LBTJ(iTJ):A.mesh.RBTJ(iTJ))*1e-18;
        holeA(A.mesh.LBTJ(iTJ):A.mesh.RBTJ(iTJ))=holeA(A.mesh.LBTJ(iTJ):A.mesh.RBTJ(iTJ))+A.mode.dop_am(A.mesh.LBTJ(iTJ):A.mesh.RBTJ(iTJ))*1e-18;
    end
end

if isfield(B.mesh,'IBTJ')==1 % expand for MTJ! (each one has a different sigma)
    if B.mode.flgBTJ_lithographic>0
        [~,iRagTJ]=min(abs(B.mesh.xgrid*1e4-B.mode.rAperture));
    else
        iRagTJ=B.mesh.nnxQW{1};
    end
    for iTJ=1:iRagTJ
        elecB(B.mesh.LBTJ(iTJ):B.mesh.RBTJ(iTJ))=elecB(B.mesh.LBTJ(iTJ):B.mesh.RBTJ(iTJ))+B.mode.dop_dp(B.mesh.LBTJ(iTJ):B.mesh.RBTJ(iTJ))*1e-18;
        holeB(B.mesh.LBTJ(iTJ):B.mesh.RBTJ(iTJ))=holeB(B.mesh.LBTJ(iTJ):B.mesh.RBTJ(iTJ))+B.mode.dop_am(B.mesh.LBTJ(iTJ):B.mesh.RBTJ(iTJ))*1e-18;
    end
end

ABSh=P.mode.ABSh;
ABSe=P.mode.ABSe;
    
alphaA=elecA*ABSe+holeA*ABSh;
alphaB=elecB*ABSe+holeB*ABSh;
alphaP=elecP*ABSe+holeP*ABSh;

% alphaA=reshape([A.VELMInput(indA).alpha],A.mesh.nny,A.mesh.nnx)*1e-2;
% alphaB=reshape([B.VELMInput(indB).alpha],B.mesh.nny,B.mesh.nnx)*1e-2;
% alphaP=reshape([P.VELMInput(indP).alpha],P.mesh.nny,P.mesh.nnx)*1e-2;

figure,
grid on,hold on,box on
plot(yA,alphaA(:,1),'linewidth',2)
plot(yB,alphaB(:,1),'linewidth',2)
plot(yP,alphaP(:,1),'linewidth',2)
axis([114 116.5 4e-1 3e3])
xlabel('z, \mum'),ylabel('Absorption coefficient \alpha, cm^{-1}')
set(gca,'yscale','log','Fontname','Times New Roman','FontSize',16)
yticks([1e0 1e1 1e2 1e3])
