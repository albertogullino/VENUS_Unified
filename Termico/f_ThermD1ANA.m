function [deltaT,PTherm,T_Contributi,fattore_correttivo,condzTe] = f_ThermD1ANA(mesh,mode,StrTT,IPLOT)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of pointers to nodes and equations and related masks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nAir = 0; 
L_air = 0;%200e-7;
nn=mesh.nny-2+nAir;         % pointer to Heat equation
nodeAir = L_air/nAir*ones(1,nAir);
nodeAir = mesh.ygrid(end-1)+cumsum(nodeAir);
node = [mesh.ygrid(2:end-1), nodeAir];  % remove the contacts
%
% Defining indexes for assembling FEM and FD matrices
in=1:nn; % node index
in1=in(1:end-1); in2=in(2:end); % left, right nodes of each element
inr=[in1 in1 in2 in2]; inc=[in1 in2 in1 in2]; in12=[in1 in2];
maskr=true(1,nn); 
mask_in12=maskr(in12); in12=in12(mask_in12);
%
% Defining weighting areas to compute integrals when assembling matrices
mesh.Le = diff(node);
mesh.Lp = zeros(1,nn);
mesh.Lp(1:nn-1) = mesh.Le/2;
mesh.Lp(2:nn) = mesh.Lp(2:nn) + mesh.Le/2;
Lp1=mesh.Lp(in1)/2; Lp1(1)=2*Lp1(1);
Lp2=mesh.Lp(in2)/2; Lp2(end)=2*Lp2(end);
%
puntatore=mesh.nny+mesh.puntatore;
%
% Heat sources (w/out contacts)
Joule = mode.HeatJoule(puntatore)*1e+12;
Thomson = 0;%mode.HeatThomson(puntatore)*1e+12;
Rec_RAD = mode.HeatRec_RAD(puntatore)*1e+12; 
Rec_srhAu = mode.HeatRec_13(puntatore)*1e+12;
OptAbs = mode.HeatOptAbs(puntatore)*1e+12;
Rec_Cap = mode.HeatCap(puntatore)*1e+12;

TotalHeat = zeros(1,nn);
TotalHeat(1:nn) = (Joule + Thomson + Rec_Cap + Rec_srhAu + OptAbs + Rec_RAD); % 

% if isfield(mesh,'IBTJ')==1
% %     for iTJ=1:length(mesh.LBTJ)  
%         iTJ=1;
%         TotalHeat(mesh.LBTJ(iTJ):mesh.RBTJ(iTJ))=mode.HeatTJ(iTJ+1)/mode.CarrierNorm;
% %     end
% end

%
if mode.Oflg==1
    PTherm=TotalHeat(2:end)*diff(node')*1e3; % mW/cm2 VENUS
    fattore_correttivo=(mode.PDissPred(end))/PTherm/mode.AreaOx; % divided by area to convert PDiss in mW/cm^2
else
    fattore_correttivo=1;
end
fattore_correttivo
%
% Thermal conductivity
Bcondz = zeros(1,nn);
if abs(StrTT.Tbuf_dd-StrTT.Tbuf)<0.2
    Bcondz(mesh.ygrid(2:mesh.nny-1)*1e4<StrTT.Tbuf_dd)=StrTT.Bcond.CondZb*mesh.fCondTerZ*1e4*mode.FatQtot;
    Bcondz(mesh.ygrid(2:mesh.nny-1)*1e4>=StrTT.Tbuf_dd)=StrTT.Bcond.CondZm*1e4*mesh.fCondTerZ*mode.FatQtot/mode.FatV;
else
    Bcondz(mesh.ygrid(2:mesh.nny-1)*1e4<StrTT.Tbuf)=StrTT.Bcond.CondZb*1e4*mode.FatQtot;
    Bcondz(mesh.ygrid(2:mesh.nny-1)*1e4>=StrTT.Tbuf)=StrTT.Bcond.CondZm*1e4*mesh.fCondTerZ*mode.FatQtot/mode.FatV;
end
Bcondz_air = StrTT.Bcond.Cond_air*1e4;
%
% Transverse thermal conducibility
condz = [Bcondz, Bcondz_air*ones(1,nAir)];
%
% Temperature dependence
% betaT = -1; % -1.3 VENUS 
% betaT = 0;    % removed non linearity of thermal conducibility
betaT=mode.Exp_Temp0;

iterCount=0;
itermaxTherm=10;
% itermaxTherm=1; % removed non linearity of thermal conducibility
TolTherm=1e-3;
er=inf; 
deltaT_old=0;

T_cond=mesh.T(2:mesh.nny-1);

fprintf('Temperature - Thermal cond. self-consistency:\n')

while (iterCount<itermaxTherm && er>TolTherm)
    
    condzT = condz.*(T_cond./mode.T300).^betaT;
    %
    % From nodal to element quantities
    condzTe = (condzT(in1)+condzT(in2))/2;
    
    % System assembling
    A_1e1e=+condzTe./mesh.Le;
    A_1e2e=-condzTe./mesh.Le;
    A_2e1e=-condzTe./mesh.Le;
    A_2e2e=+condzTe./mesh.Le;
    MM=[A_1e1e, A_1e2e, A_2e1e, A_2e2e];
    
    Amat = sparse(inr,inc,MM,nn,nn);
    
    TotalHeat1 = TotalHeat(in1);
    TotalHeat2 = TotalHeat(in2);
    MM=[TotalHeat1.*Lp1, TotalHeat2.*Lp2];
    
    qvet=sparse(in12,1,MM(mask_in12),nn,1);
    %
    % Enforce the boundary condition
    qvet(1)= 0;
%     qvet(end)=100;
    %
    % Boundary conditions: Heat sink temperature (w.r.t. device)
    Amat(1,:) = 0;
    Amat(1,1) = 1; % line contact on the left
%     Amat(end,:) = 0;
%     Amat(end,end) = 1; % line contact on the right
    %
    % Solve the system
    deltaT = Amat\(fattore_correttivo*qvet);
    
    T_cond = mode.T0+deltaT';
        
    er=abs(1-max(deltaT_old)/max(deltaT));
    if isnan(er)
        er=0;
    end
    
    deltaT_old=deltaT;
    iterCount=iterCount+1;

    fprintf('  ')
    fprintf(['Iteration: ',num2str(iterCount),'  |  Error: ',num2str(double(er),'%10.5e'),' \n'])
    
end

T_Contributi{1}=f_ThermD1ANA_Contributi(mesh,condzTe,Joule,fattore_correttivo);
T_Contributi{2}=f_ThermD1ANA_Contributi(mesh,condzTe,Rec_srhAu,fattore_correttivo);
T_Contributi{3}=f_ThermD1ANA_Contributi(mesh,condzTe,Rec_Cap,fattore_correttivo);
T_Contributi{4}=f_ThermD1ANA_Contributi(mesh,condzTe,Rec_RAD,fattore_correttivo);
T_Contributi{5}=f_ThermD1ANA_Contributi(mesh,condzTe,OptAbs,fattore_correttivo);

deltaT=full(deltaT);
deltaT(end)=deltaT(end)*mode.FatQcontact;

deltaT=[0; deltaT; 0];
deltaT=repmat(deltaT,1,mesh.nnx);

% mode.condzTe=condzTe;
% mode.deltaT = deltaT(1:mesh.nn);
% mode.DeltaTmax(mode.indV)=max(deltaT);

if IPLOT==1
    figure(222),clf
    set(gcf,'position',[69 59.4 560 420])
    subplot(121)
    grid on,box on
    plot(mesh.ygrid(2:end-1)*1e4,deltaT(2:end-1,2),'b','linewidth',2)
    xlabel('z, \mum'),ylabel('\DeltaT(z), K')
    set(gca,'FontSize',12)

    subplot(122)
    grid on,box on
    plot(mesh.ygrid(2:end-1)*1e4,deltaT(2:end-1,2),'b','linewidth',2)
    xlim([StrTT.Tbuf mesh.ygrid(end)*1e4])
    xlabel('z, \mum'),ylabel('\DeltaT(z), K')
    set(gca,'FontSize',12)

    drawnow
end

fprintf('  ')
fprintf(['Converged! Max. T variation: ',num2str(max(double(deltaT(:,1)))),' K \n'])
% keyboard

% mode.T = mesh.T+mode.deltaT';

