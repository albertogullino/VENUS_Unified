clear 
close all

load out\LW_MarkusN_FINALE_LgDBR_20C_Stephan.mat
% load out\LW_MarkusN_FINALE_LgDBR_20C_Stephan_repelem.mat

Doggi_CallTHERM

nnx=mesh.nnx;   % mesh points along rho (x)
nny=mesh.nny;   % mesh points along z (y)

z=y;    % mesh along z (y)
rho=x;  % mesh along rho (x)

Z=repmat(z,1,nnx)';     % Stacking of z mesh
R=repelem(rho,nny)';    % Stacking of rho mesh

% % Same mesh of thermal solver
% Z=repmat(zz,1,length(xx))';     % Stacking of z mesh
% R=repelem(xx,length(zz))';    % Stacking of rho mesh

% I=mode.ii_dd'*1e3;          % Current (mA)
V=mode.vv_dd';              % Voltage (V)
Pst=sum(mode.Pst_dd,1)';    % Output optical power (mW)

% Simplified source
Pin=mode.ii_dd(iCurr-1)*mode.vv_dd(iCurr-1);    % Input power (W)
Pout=Pst(iCurr-1)*1e-3;                                  % Output power (W)  

if iFakeFlg==1
    % FIXED (values used in "20231115_HeatSources_MHerper.pptx")
    zbot=114;              % um
    ztop=116.5;            % um
    [~,izbot]=min(abs(zbot-y)); % extract the corresponding mesh index
    [~,iztop]=min(abs(ztop-y)); % extract the corresponding mesh index
    
    dz=y(iztop)-y(izbot);       % Length along z
    rSimpl=StrTT.Rox;           % um
%     rSimpl=StrTT.ro_max;        % um (radial point)

    % Compute the simplified source (Stephan way)
    Qfake=(Pin-Pout)/(pi*rSimpl^2*(y(iztop)-y(izbot)))      % W/um^3

end


% Maximum T variation from each source 
% (commented the saved results from VENUS)
DeltaTmax=mode.DeltaTmax';
DeltaTmax_Joule=mode.DeltaTmax_Joule';
DeltaTmax_FCA=mode.DeltaTmax_OptAbs';
DeltaTmax_Ccap=mode.DeltaTmax_Ccap';
DeltaTmax_rad=mode.DeltaTmax_RAD';
DeltaTmax_nr=mode.DeltaTmax_srhAu';

% DeltaTmax=max(max(DeltaT));
% DeltaTmax_Joule=max(max(T_Contributi{1}));
% DeltaTmax_nr=max(max(T_Contributi{2}));
% DeltaTmax_Ccap=max(max(T_Contributi{3}));
% DeltaTmax_rad=max(max(T_Contributi{4}));
% DeltaTmax_FCA=max(max(T_Contributi{5}));

% Compute the thermal power from each source, as a function of current 
PthermIV

PlotThermalSources

% Thermal sources (normalized by the corresponding dissipated power)
Joule=reshape(Joule/PJoule(iCurr-1),1,mesh.nnx*mesh.nny)';
FCA=reshape(FCA/Pfca(iCurr-1),1,mesh.nnx*mesh.nny)';
Ccap=reshape(Ccap/Pcap(iCurr-1),1,mesh.nnx*mesh.nny)';
GR=reshape(GR/(Pnr(iCurr-1)+Prad(iCurr-1)),1,mesh.nnx*mesh.nny)';

DeltaT=reshape(DeltaT,1,mesh.nnx*mesh.nny)';

fcorr=mode.fattore_correttivo';


if mode.Pmin_Pfit>20
    iv=find(mode.vv_dd<mode.VIdrive);
    ii=find(mode.vv_dd>=mode.VIdrive);

    ibol=([0 diff(mode.vv_dd(iv)) diff(mode.ii_dd(ii))*1e3]>0);
    ibol(1)=true;

    I=I(ibol);
    V=V(ibol);
    Pst=Pst(ibol);
    PTherm=PTherm(ibol);
    PJoule=PJoule(ibol); Pfca=Pfca(ibol); Pnr=Pnr(ibol); Prad=Prad(ibol); Pcap=Pcap(ibol);
    DeltaTmax=DeltaTmax(ibol);
    DeltaTmax_Joule=DeltaTmax_Joule(ibol); DeltaTmax_FCA=DeltaTmax_FCA(ibol); 
    DeltaTmax_nr=DeltaTmax_nr(ibol); DeltaTmax_rad=DeltaTmax_rad(ibol); DeltaTmax_Ccap=DeltaTmax_Ccap(ibol);
    fcorr=fcorr(ibol);

    pausak
end


% save Data_VENUSThermalSources-2mA z rho Z R I V Pst Pin Pout DeltaT DeltaTmax DeltaTmax_Ccap DeltaTmax_FCA DeltaTmax_Joule DeltaTmax_nr DeltaTmax_rad Joule FCA GR Ccap PTherm PJoule Pfca Prad Pnr Pcap rSimpl dz

% filename='20231127_GullinoData_SimpleSource.xlsx';
filename='202400214_NonSelf_GullinoData.xlsx';

pausak 

NonMatrixData = [I, V, I.*V, Pst, I.*V-Pst, PTherm.*fcorr, fcorr, PTherm, PJoule, Pfca, Pnr+Prad, Pcap, DeltaTmax, DeltaTmax_Joule, DeltaTmax_FCA, DeltaTmax_nr + DeltaTmax_rad, DeltaTmax_Ccap];
NMD=array2table(NonMatrixData,'VariableNames',{'Current (mA)','Voltage (V)','Pin (mW)','Pout (mW)','Pdiss=Pin-Pout (mW)','Ptherm_corr (mW)','fcorr','PTherm (Diss. power) (mW)','Diss. power - Joule (mW)','Diss. power - FCA (mW)','Diss. power - GR (mW)','Diss. power - Ccap (mW)','DeltaT - total (K)','DeltaT - Joule (K)','DeltaT - FCA (K)','DeltaT - GR (K)','DeltaT - Ccap (K)'});
writetable(NMD,filename,'sheet','Current-dependent quantities')

pausak
% In case the fake source used: Add column referring to DeltaT from Qfake
if iFakeFlg==1
    SimpleSourceData = [Pin, Pout, dz, zbot, ztop, rSimpl, Qfake];
    SSD=array2table(SimpleSourceData,'VariableNames',{'Pin (W)','Pout (W)','dz (um)','zbot (um)','ztop (um)','rSimpl (um)','Qfake (W/um3)'});
    writetable(SSD,filename,'sheet',[num2str(Pcor),'mA_SimpleSourceData'])
    pausak 

    DeltaTsimple=DeltaT;
%     DeltaTsimple=reshape(Tprec_new,1,length(xx)*length(zz))';
    DTf=array2table([Z,R,DeltaTsimple],'VariableNames',{'z (um)','rho (um)','DeltaTsimple (K)'});
    writetable(DTf,filename,'sheet',[num2str(Pcor),'mA_SimpleSourceData'])
else
    MatrixData = [Z, R, DeltaT, Joule, FCA, GR, Ccap];
    MD=array2table(MatrixData,'VariableNames',{'z (um)','rho (um)','DeltaT (K)','Joule (W/um3)','FCA (W/um3)','GR (W/um3)','Ccap (W/um3)'});
    writetable(MD,filename,'sheet',[num2str(Pcor),'mA_MeshTvariationHeatSources'])
end
% xlswrite('Data_VENUSThermalSources-Curr.xlsx',NonMatrixData)
% % xlswrite('Data_VENUSThermalSources-Map.xlsx',MatrixData)
% xlswrite('aprova.xlsx',MatrixData)
% xlswrite('Data_VENUSThermalSources-SimpleQ.xlsx',SimpleSourceData)
