
function mode=f_EvalHeatingTerms(geom,mesh,mode)

global TFasano EFasano xmolFasano alphaFasano

% funzione che calcola mediante post-processing di geom, mesh, mode, alcuni
% termini collegati al heating, e non solo
if mode.Elementi==1
    [Jn_x,Jn_y,Jp_x,Jp_y] = f_EvalCurrentDensityElementi(geom,mesh,mode);
else
    [Jn_x,Jn_y,Jp_x,Jp_y] = f_EvalCurrentDensity(geom,mesh,mode);
end

% Costanti
s_LoadConstants

% Calcolo resistività
mobn_n=mesh.mobn0_n;
mobp_n=mesh.mobp0_n;
sigma=qel.*(mode.elec.*mobn_n+mode.hole.*mobp_n);
if isfield(mesh,'IBTJ')==1 % expand for MTJ! (each one has a different sigma)
    if mode.flgBTJ_lithographic>0
        [~,iRagTJ]=min(abs(mesh.xgrid*1e4-mode.rAperture));
    else
        iRagTJ=mesh.nnxQW{1};
    end
    if mode.flgHeatTJ==1
        % Equivalent sigma for the TJ nodes is introduced in place of the
        % DD one!
        for iTJ=1:iRagTJ
            % sigmaTJ has been computed
            sigma(mesh.LBTJ(iTJ):mesh.RBTJ(iTJ))=mode.sigmaTJ(iTJ,1)/mode.CarrierNorm;
            %        sigma(mesh.LBTJ(iTJ):mesh.RBTJ(iTJ))=mode.sigmaTJ(iTJ)/mode.CarrierNorm;
        end
    end
end

% Calcolo termine Joule heating
HeatJoule_x=((Jn_x+Jp_x).^2)./sigma;
HeatJoule_y=((Jn_y+Jp_y).^2)./sigma;
if mode.quasi1D==1
    HeatJoule=HeatJoule_y;
    HeatJoule([mesh.nny-2 2*mesh.nny-2])=0;
else
    HeatJoule=HeatJoule_x+HeatJoule_y;
end

% Computes an equivalent HeatTJ with VTJ and ITJ (not used)
if isfield(mesh,'IBTJ')==1 && mode.flgHeatTJ==2
    % Equivalent HeatTJ=VTJ*Itot for the TJ nodes is introduced in place of
    % the DD one!
    for iTJ=1:iRagTJ
        HeatJoule(mesh.LBTJ(iTJ):mesh.RBTJ(iTJ))=squeeze(mode.HeatTJ(1,end,iTJ))/mode.CarrierNorm;
    end
end

% Per mia comodità, riporto sui nodi
mode.Jn_x=Jn_x;
mode.Jn_y=Jn_y;
mode.Jp_x=Jp_x;
mode.Jp_y=Jp_y;
mode.sigma=sigma;

% Calcolo Thomson heating
% Finding elements not in semiconductor regions
semiconductor=find(geom.semiconductor);
it=not(ismember(mesh.triangle(4,:),semiconductor));
%
% Initializing relevant quantities
sn=1; sp=1; % Sentaurus manual
kappan=1; kappap=1; % Sentaurus manual
Pn = -kappan.*kB./qel.*(5/2-sn+log(mesh.Nc./mode.elec)); % ThermoElectric Power, V/K
Pp = +kappap.*kB./qel.*(5/2-sp+log(mesh.Nv./mode.hole)); % ThermoElectric Power, V/K
%
dPn_t_x=Pn*mesh.gradx; dPn_t_x(it)=0; % V/(cm*K)
dPn_t_y=Pn*mesh.grady; dPn_t_y(it)=0; % V/(cm*K)
dPp_t_x=Pp*mesh.gradx; dPp_t_x(it)=0; % V/(cm*K)
dPp_t_y=Pp*mesh.grady; dPp_t_y(it)=0; % V/(cm*K)
%
dPn_x=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),dPn_t_x),1,mesh.nn);
dPn_y=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),dPn_t_y),1,mesh.nn);
dPp_x=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),dPp_t_x),1,mesh.nn);
dPp_y=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),dPp_t_y),1,mesh.nn);
%
HeatThomson_x=-Jn_x.*mesh.T.*dPn_x-Jp_x.*mesh.T.*dPp_x; % W/(cm^3)
HeatThomson_y=-Jn_y.*mesh.T.*dPn_y-Jp_y.*mesh.T.*dPp_y; % W/(cm^3)
%
% W/cm^3 -> W/um^3 : moltiplicazione per 1e-12;
mode.HeatJoule_x=HeatJoule_x*1e-12*mode.CarrierNorm;
mode.HeatJoule_y=HeatJoule_y*1e-12*mode.CarrierNorm;
mode.HeatJoule=HeatJoule*1e-12*mode.CarrierNorm;

mode.HeatRec_nr_bulk=qel.*mesh.Eg.*(mode.RSRH+mode.RAuger)*1e-12*mode.CarrierNorm;
mode.HeatRec_rad_bulk=qel.*mesh.Eg.*(mode.Rrad)*1e-12*mode.CarrierNorm;
mode.HeatThomson=(HeatThomson_x+HeatThomson_y)*1e-12*mode.CarrierNorm;
%'thom',keyboard
if(mode.oflg)
    mode.HeatRec_nr_QW=qel.*mesh.Eg.*(mode.RSRHQW+mode.RAugerQW)*1e-12*mode.CarrierNorm;
    mode.HeatRec_rad_QW=qel.*mesh.Eg.*(mode.RradQW)*1e-12*mode.CarrierNorm;
    mode.HeatRec_nrParziale=mode.HeatRec_nr_QW+mode.HeatRec_nr_bulk;
    mode.HeatRec_rad=mode.HeatRec_rad_QW+mode.HeatRec_rad_bulk;
    mode.HeatCapn=qel.*mesh.DeltaEcQW{1}.*(mode.Ccapn3D)*1e-12*mode.CarrierNorm;
    mode.HeatCapp=qel.*mesh.DeltaEvQW{1}.*(mode.Ccapp3D)*1e-12*mode.CarrierNorm;
    mode.HeatCap=mode.HeatCapn+mode.HeatCapp;
else
    mode.HeatRec_nr=mode.HeatRec_nr_bulk;
    mode.HeatRec_rad=mode.HeatRec_rad_bulk;
    mode.HeatRec_nrParziale=0;
    mode.HeatCap=0;
    mode.HeatRec_rad=0;
end
fat_RAD=mode.fat_RAD;
mode.HeatRec_nr=mode.HeatRec_nrParziale+mode.HeatCap+mode.HeatRec_rad*fat_RAD;
mode.HeatRec_Cap=mode.HeatCap;
mode.HeatRec_RAD=mode.HeatRec_rad*fat_RAD;
mode.HeatRec_13=mode.HeatRec_nrParziale;%*0;
%mode.HeatTh=mode.HeatThomson;
%'thom',keyboard

% R=reshape(mode.HeatRec_nrParziale,mesh.nny,mesh.nnx);
% figure, plot(mesh.xgrid,R)
% figure, plot(mesh.ygrid,R)
% Rr=reshape(mode.HeatRec_rad,mesh.nny,mesh.nnx);
% figure, plot(mesh.xgrid,Rr)
% figure, plot(mesh.ygrid,Rr)
% Rt=reshape(mode.HeatThomson,mesh.nny,mesh.nnx);
% figure, plot(mesh.xgrid,Rt)
% figure, plot(mesh.ygrid,Rt)


% Drude model

% questa parte è legata ai termini di assorbimento tramite il alpha, che
% costituisce un assorbimento ottico. per calcolare esplicitamente questo
% assorbimento ottico serve la standing wave. in aggiunta a questo, questi
% termini permettono di calcolare anche una variazione dell'indice di
% rifrazione del materiale dovuta ai portatori

%
% Define N2D, P2D charges from Poisson charge
% N2D=mode.N2D*1e6; P2D=mode.P2D*1e6;
% indQW=length(mesh.meffnQW); indQW=ceil(indQW/2);
% meffnQW=mesh.meffnQW(indQW); meffpQW=mesh.meffpQW(indQW);
%
% Converting electrons and holes in m^(-3)
elec=mode.elec*1e6*mode.CarrierNorm;
hole=mode.hole*1e6*mode.CarrierNorm;
%
% Converting mobilities from cm^2/(s V) to m^2/(s V)
mobn=mobn_n/1e4;
mobp=mobp_n/1e4;
%
% Saving refractive index from static: this is OK in GaAs, not general!
nrefr=sqrt(mesh.epsxx_n./eps0);

% Re-defining physical constants in meters
Clight=299792458; % Speed of light (m/s)
mu0=4*pi*1e-7; % Magnetic permeability (H/m)
eps0=1/(mu0*Clight^2); % Dielectric permittivity (F/m)
%
omega=2*pi*Clight/(mean(mode.vlambda)*1e-9); % average lambda


iprovaDrude=0;
if iprovaDrude==1
    %%%%%%% TUTTO IN cm
    % Converting electrons and holes in m^(-3)
    elec1=mode.elec;
    hole1=mode.hole;
    %
    % Converting mobilities from cm^2/(s V) to m^2/(s V)
    mobn1=mobn_n;
    mobp1=mobp_n;
    % Re-defining physical constants in meters
    Clight1=299792458*100; % Speed of light (cm/s)
    mu01=4*pi*1e-7/100; % Magnetic permeability (H/cm)
    eps01=1/(mu01*Clight1^2); % Dielectric permittivity (F/cm)
    %
    omega1=2*pi*Clight1/(mean(mode.vlambda)*1e-7); % average lambda
    alphaDrude=sqrt(mu0/eps0)./nrefr.*qel.^3./(omega.^2.*m0.^2).*(elec./(mesh.meffn.^2.*mobn)+hole./(mesh.meffp.^2.*mobp));
    alphaDe=sqrt(mu0/eps0)./nrefr.*qel.^3./(omega.^2.*m0.^2).*(1./(mesh.meffn.^2.*mobn));
    alphaDh=sqrt(mu0/eps0)./nrefr.*qel.^3./(omega.^2.*m0.^2).*(1./(mesh.meffp.^2.*mobp));
    
    alphaDe1=sqrt(mu01/eps01)./nrefr.*qel.^3./(omega1.^2.*m0.^2).*(1./(mesh.meffn.^2.*mobn1));
    alphaDh1=sqrt(mu01/eps01)./nrefr.*qel.^3./(omega1.^2.*m0.^2).*(1./(mesh.meffp.^2.*mobp1));
    
end
% mode.Deltan=-qel.^2./(2.*nrefr.*eps0.*omega.^2.*m0).*(elec./mesh.meffn+hole./mesh.meffp+N2D./meffnQW+P2D./meffpQW);

if isfield(mode,'FatDrudeRef')
    FatDrudeRef=mode.FatDrudeRef;
else
    FatDrudeRef=1;
end

% Drude model 
mode.Deltan=-FatDrudeRef*qel.^2./(2.*nrefr.*eps0.*omega.^2.*m0).*(elec./mesh.meffn+hole./mesh.meffp)*mode.CarrierNorm^(-2/3);
mode.Deltan(not(mesh.iq))=0;

% 'alpha TAROCCATO'
AlTarocco=mode.AlTarocco;
vEF=mode.EFn-mode.ecb;
vEF(isnan(vEF))=10; % to treat contacts or parts where NaN should appear
% alphainterp: cm^(-1)
% save alpha
%'calcolo di alpha',keyboard

%'qui Dru', keyboard
if mode.ABSh+mode.ABSe==0
    % 'vecchio Fasano', keyboard
    alphaFasanoInterp=interp3(TFasano,EFasano,xmolFasano,alphaFasano,mesh.T,vEF,mesh.xmol,'spline',0)*1e2;
    alphaDrude=sqrt(mu0/eps0)./nrefr.*qel.^3./(omega.^2.*m0.^2).*(elec./(mesh.meffn.^2.*mobn)+hole./(mesh.meffp.^2.*mobp));
    mode.alpha=AlTarocco.*(alphaFasanoInterp+alphaDrude);
    mode.alpha(not(mesh.iq))=0;
    %mode.alphaFasano=alphaFasanoInterp;
    %mode.alphaDrude=alphaDrude;
else
    ifFC=mode.ifFC;
    if ifFC==1
        elecABS=mode.elec;
        holeABS=mode.hole;
        
        if mode.flgBTJ==1
            for iTJ=1:iRagTJ
                elecABS(mesh.LBTJ(iTJ):mesh.RBTJ(iTJ))=elecABS(mesh.LBTJ(iTJ):mesh.RBTJ(iTJ))+mode.dop_dp(mesh.LBTJ(iTJ):mesh.RBTJ(iTJ));
                holeABS(mesh.LBTJ(iTJ):mesh.RBTJ(iTJ))=holeABS(mesh.LBTJ(iTJ):mesh.RBTJ(iTJ))+mode.dop_am(mesh.LBTJ(iTJ):mesh.RBTJ(iTJ));
            end
        end
    elseif ifFC==0
        elecABS=mode.dop_dp;
        holeABS=mode.dop_am;
    elseif ifFC==2
        elecABS=(mode.dop_dp+mode.elec)/2;
        holeABS=(mode.dop_am+mode.hole)/2;
    end
    elecABS=elecABS*1e-18*mode.CarrierNorm;
    holeABS=holeABS*1e-18*mode.CarrierNorm;
    
    mode.elecABSvelm=elecABS;   % used in VELM, T dependence still NOT included!
    mode.holeABSvelm=holeABS;
    
    T300=mode.T300;
    Temp=mesh.T-T300;
    
    Al=mode.ABS_Apor;
    
    % This is done TWICE: here and in f_CallVELM, with the same parameters
    % (there was a bug, such that (1+Temp/T300).^ABS_Texp was accounted
    % twice in VELM. This bug is now corrected, and ABS_Texp=0;
    % In summary: if you want to correct/change this formulas, do the same
    % in VELM!!!
    ABS_Texp=mode.ABS_Texp+mode.PerCoefExT*Temp;
    
    if Al>0
        N0=1; % per N=1 tutte coincidono
        InDe=1/log(Al*N0+1);
        if isfield(mode,'IOLDsw')
            if mode.IOLDsw==1
                N0=mode.ABS_Ader;
            end
            InDe=(Al*N0+1)/Al;
        end
%         InDe=1/log(Al+1);
        mode.elecABS=InDe*log(Al*elecABS+1).*(1+Temp/T300).^ABS_Texp;
        mode.holeABS=InDe*log(Al*holeABS+1).*(1+Temp/T300).^ABS_Texp;
    else
        mode.elecABS=elecABS.*(1+Temp/T300).^ABS_Texp;
        mode.holeABS=holeABS.*(1+Temp/T300).^ABS_Texp;
    end
    
    if isfield(mode,'Fat_PerCoefTemp')
        Fat_Perd=mode.Fat_Perd0+Temp*mode.Fat_PerCoefTemp;
        ABSh=mode.ABSh0.*Fat_Perd;
        ABSe=mode.ABSe0.*Fat_Perd;
    else
        ABSh=mode.ABSh;
        ABSe=mode.ABSe;
    end
    
    alN3=(mode.holeABS.*ABSh+mode.elecABS.*ABSe)*100;
%     alN=holecentro*7+eleccentro*3;
    mode.alpha=AlTarocco.*alN3;
    mode.alpha(not(mesh.iq))=0;
end

% FCA compuation
if(mode.oflg)
    Pst=reshape(mode.Pst,mode.nmodes,1)/1000; % mW
    
    % OptAbsorption (um)
    E2=mode.E2/1e8; % conversion cm^{-2} to um^{-2}
    
    for indMode=1:mode.nmodes
        if mode.quasi1D==1
            % 1D version excludes the contact nodes (not playing any role)
            OptPowerModes(indMode,:)=reshape([0 mode.Inviluppo_SW 0].'*E2(indMode,:),1,mesh.nn).';
        else
            OptPowerModes(indMode,:)=reshape(mode.Inviluppo_SW.'*E2(indMode,:),1,mesh.nn).';
        end
    end
    mode.OptPowerDensity=reshape(OptPowerModes.'*Pst,1,mesh.nn);
    mode.HeatOptAbs=mode.OptPowerDensity.*(mode.alpha*1e-6);    
end


return

