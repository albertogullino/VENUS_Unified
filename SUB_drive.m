mode.strName=strSave;
mode.structureName=structureName;

CurDyn=[];
flagCallVELM=1;

% if IPLOT==1
%     portfig=figure(kpar);
%     set(portfig,'pos',[1546         611         364         373])
% end

if IPAR==41
    if kpar==1
        save iplot IPLOT
    else
        load iplot
    end
end

%========================================================================76
% solve at equilibrium
%========================================================================76

if irest==0
    
    mesh=rectmesh(geom); % create tensorial mesh
    if mode.quasi1D==1
        %         mesh.puntatore=mesh.nny+(2:mesh.nny-1);
        mesh.puntatore=(2:mesh.nny-1);
    end
    
    T=mode.T0*ones(1,mesh.nn);
    
    mesh.DeltaT=T-mode.T0;
    % mesh.DeltaT=0;
    
    mesh.dz_Term=ParMore.dz_Term;
    
    mode.DeltaTmax=0;
    
    fil_str=[structureName,'.str'];
    mode.matgain=0;
%     mode.vlambda=850*ones(NUMERO_MODI,1); % inizializzazione
    mode.vlambda=ParOpt.lambdaNot*1e9*ones(NUMERO_MODI,1); % inizializzazione
    
    % Temporaneo: da sostituire con ciò che sarà preso dal codice Density
    % Matrix
    Dn=mode.matgain;
    mode.DeltaN=0;
    
    mode.v0_dd=dd;
    
    s_LoadConstants
    %'load', keyboard
    TKT=T;
    if isfield(mode,'KT')
        if mode.KT==0
            TKT=mode.T0*ones(size(T));
        end
    end
    
    Vt=kB.*TKT./qel;
    mode.Vt=Vt;
    
    % Material dependent quantities are computed and extracted
    mesh=loadmacro(geom,mesh,mode,T);
    
    mesh.Tgain=mesh.T;
    
    if IPLOT<=-1
        fprintf('Call "plotStructure_PLOT.m" to visualize mesh structure\n')
        
        plotStructure_IPLOT
        
        pausak
        if IPLOT==2
            % if IPLOT=2, only structure details are shown
            IPLOT=0; 
        end
    end
    
    fprintf('\n Thermodynamic Equilibrium (Poisson eq. solved)\n')
    mode=p_solve1(geom,mesh,mode); % first run: thermodynamic equilibrium
    
    % Store equlibrium quantities
    mode.n3Di=mode.elec;
    mode.p3Di=mode.hole;
    mode.phi_0=mode.phi;
    mesh.nint_equil=mesh.nint;
    % Carriers for FCA heating source (for the losses)
    mode.elecABS=mode.elec*1e-18*mode.CarrierNorm;
    mode.holeABS=mode.hole*1e-18*mode.CarrierNorm;
    
    mode.elecABSvelm=mode.elecABS;   % used in VELM
    mode.holeABSvelm=mode.holeABS;

    
    % nmodes=length(mode.Lm); % number of optical modes (VELM)
    mode.nmodes=NUMERO_MODI;
    mode.lambda=ones(mode.nmodes,1)*ParOpt.lambdaNot*1e9;
    
    % To initialize matgain, we use carriers in the central well
    %indCenter=floor((mesh.NMQW+1)/2);
    %nQW=median(mode.n2D{indCenter});
    %pQW=median(mode.p2D{indCenter});
    %TQW=300*ones(size(nQW));
    %[g_eq] = f_InterpGainEH(nQW,pQW,mode.Deltalam+mode.vlambda(1)*ones(size(nQW)),TQW);
    %mode.matgain=g_eq;
    g_eq=0;
    mode.matgain=g_eq;
    
    % initializations
    mode.iconv=0;
    semiconductor=find(geom.semiconductor);
    it=not(ismember(mesh.triangle(4,:),semiconductor));
    
    % compute electric field (in time domain), V/um
    efieldx_t=-mode.phi*mesh.gradx*1e-4; efieldx_t(it)=0;
    efieldy_t=-mode.phi*mesh.grady*1e-4; efieldy_t(it)=0;
    mode.efield_x=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),efieldx_t),1,mesh.nn);
    mode.efield_y=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),efieldy_t),1,mesh.nn);
    mode.Deltan=0; % Drude Delta refractive index (real part)
    mode.alpha=0;  % Drude Delta refractive index (imag part)
    
    
    % =============================================================================================100
    % solution vector uvet: {phi,elec,hole,elec_t,v_dd,i_dd}
    % =============================================================================================100
    % vector indexes:
    %
    %    in    -> pointer to Poisson equation
    %    in+nn -> pointer to electron continuity equation
    %    in+pp -> pointer to hole continuity equation
    %    in+tt -> pointer to trap equations
    %    in+qq -> pointer to circuit equations
    %    in+rr -> pointer to current equations
    %    iv+vv -> pointer to 2d electron continuity equations
    %    iv+ww -> pointer to 2d hole continuity equations
    %    in+ss -> pointer to photon rate equations
    %
    nm=geom.nm;               % number of active contacts
    if(isfield(mesh,'nnxQW'))
        nnQW=mesh.nnxQW{1};         % number of lateral points in the QW
    else
        nnQW=0;
    end
    NQW=mesh.NMQW;                % number of quantum wells
    %
    nl=mode.ntrap;                % number of trap levels
    
    if mode.EqPermutazione==0
        %
        nn=mesh.nn;                   % -> elec eqs.
        pp=nn+mode.nflg*nn;           % -> hole eqs.
        tt=pp+mode.pflg*nn;           % -> trap eqs.
        qq=tt+mode.tflg*nl*nn;        % -> circuit eqs.
        rr=qq+nm;                     % -> current eqs.
        vv=rr+nm;                     % -> n2D eqs.
        ww=vv+length(mesh.xgrid)*NQW;               % -> p2D eqs.
        ss=ww+length(mesh.xgrid)*NQW;               % -> optical rate eqs.
    else
        nn=mesh.nn;                             % -> elec eqs.
        pp=nn+mode.nflg*nn;                     % -> hole eqs.
        tt=pp+mode.pflg*nn;                     % -> trap eqs.
        vv=tt+mode.tflg*nl*nn;                  % -> n2D eqs
        ww=vv+length(mesh.xgrid)*NQW;           % -> p2D eqs
        ss=ww+length(mesh.xgrid)*NQW;           % -> optical rate eqs.
        qq=ss+mode.nmodes;
        rr=qq+nm;
    end
    %
    % total number of equations/unknowns
    neq=nn+mode.nflg*nn+mode.pflg*nn+mode.tflg*nl*nn+2*nm+mode.qflg*not(mode.firstrun)*2*NQW*length(mesh.xgrid)+mode.Oflg*mode.nmodes;
    
    % We first compute the solution vector with all the information we have at this point;
    % notice that here the solution is in frequency domain
    uvet=zeros(1,neq);
    uvet(1:nn)=mode.phi;
    if(mode.nflg), uvet(nn+(1:nn))=mode.elec; end
    if(mode.pflg), uvet(pp+(1:nn))=mode.hole; end
    if(isfield(mode,'v_dd')), uvet(qq+1)=mode.v_dd; end
    if(isfield(mode,'i_dd')), uvet(rr+1)=mode.i_dd; end
    if(isfield(mode,'n2D'))
        for indQW=1:NQW
            uvet(vv+(indQW-1)*length(mesh.xgrid)+(1:nnQW))=mode.n2D{indQW};
        end
    end
    if(isfield(mode,'p2D'))
        for indQW=1:NQW
            uvet(ww+(indQW-1)*length(mesh.xgrid)+(1:nnQW))=mode.p2D{indQW};
        end
    end
    
    if(mode.Oflg && isfield(mode,'Pst'))
        uvet(ss+(1:mode.nmodes))=mode.Pst*mode.Cpot/1e3;
    end
    %
    s_LoadConstants
    
    PTherm=0;
    
    mode_old=mode;
    indv=1;
    
    DeltaTempVELMLast=0;
    npMaxLast=0;
    Pst=0;
    indVELM=0;
else
    Glut4Dinc(mode)
end %irest

%'irest', keyboard
if irest>0
    if kpar==1
        save LSW Last_Workspace irest rad_setting IREST nomeSav iDyn kpar kpvet IPvet NP PMAT rad_setting rad_settingV
    else
        save LSW Last_Workspace irest rad_setting IREST nomeSav iDyn kpar kpvet IPvet NP PMAT rad_setting  rad_settingV MODEplot
    end
    %       !copy Last_Workspace.mat Last_Workspace_Old.mat
    %'qui', keyboard
    if irest==2
        eval(['!copy ', Last_Workspace,'.mat  ',  Last_Workspace,'_OLD.mat'])
    else
        eval(['!copy ', Last_Workspace,'_stop.mat  ',  Last_Workspace,'_OLDstop.mat'])
    end
    clear
    load LSW
    %        clear global
    close all
    %'ferma', keyboard
    
    %     if irest==2
    %         eval(['load ',Last_Workspace])
    %     else
    %         eval(['load ',Last_Workspace,'_stop'])
    %     end
    
    load LSW
    fis= strfind(Last_Workspace,'\');
    Last_Workspac=Last_Workspace(fis(end)+1:end);
    %'ver', keyboard
    %DirName=structure(1:fis(end));
    eval(['load ',nomeSW,'Contributi_',nomeSav])
    %        load Contributi
    ASSEGNO_Effetti
    eval(['settings_vari',rad_setting])
    if IPAR>0
        ppar=PMAT{IPAR}(kpar);
        stringa_parametro=[Svar{IPAR},'=',num2str(ppar)]
        eval([stringa_parametro,';'])
    end
    ASSEGNO_mode
    %        'Effetti cuora', keyboard
    le_VELMinfo=length(VELMInfo);
    %'ver', keyboard
    %        load Last_Workspace
    %clear MODEplot
    indv=indv+1;
    
    if irest==2
        mode.v0_dd=[mode.v0_dd,Vadd];
    end
end

if IOLDsw==0
    mesh.fCondTer=fCondTer;
    mesh.fCondTerZ=fCondTer*fCondTerZ;
else
    mesh.fCondTer=fCondTer;
    mesh.fCondTerZ=fCondTerZ;
end

if mode.flgBTJ==1
    mode.nBTJ=StrDD.nBTJ;
end

mesh.dndT=fatt_dndT*dndT;
mesh.dndT1D=dndT1D;

ContaBisez=0;   % initializes the bisection flag in case of no convergence
DTM0=mode.T300;
DeltaTmax=0;
I_mA=0;
if ~exist('Imassimo')
    Imassimo=1000;
end

CondPotBreak=0;
clear PtotV

mode.deltaT=0;
mode.IdriveON=0;            % becomes 1 when current driving is turned ON
fitStart=2; fitDegree=2;    % fitting of VELM parameters in quasi1D simulations
adiabaticCounter=1;
NVbias=length(mode.v0_dd);

iprimoRunna=1;
RunnaVELM=0;

while(indv<=NVbias && DeltaTmax<DTM0 && I_mA<abs(Imassimo) && CondPotBreak==0) % outer bias-loop %$$$$$$$$$$$$$$$$$$$$$$$$
    
    CustomDDstart % script to make the code start from any voltage (BE CAREFUL - Torrelli)
    
    elapsed_time=cputime;
    %    mode_old=mode;
    %    mesh_old=mesh;
    uvet_old=uvet;
    if mode.IdriveON==1
        i0_dd = mode.i0_dd(indv);
        fprintf(['Ibias = ',num2str(i0_dd*1e3),' mA\n'])
        v0=v_dd;
    else
        v0_dd = mode.v0_dd(indv);
        fprintf(['Vbias = ',num2str(v0_dd),' V\n'])
        v0=v0_dd;
    end
    
    if(mode.Tflg && abs(v0)>mode.minthermalvoltage && flagCallVELM) % OCCHIO!!! flagCallVELM
        Tprec=mesh.Tprec;
        % Recombination-only heating
        mode.iplotTerm=1;
        %        mode.iTfig=1;
        %           'DeltaT', keyboard
        if isfield(mode,'IsoThermal')
            IsoThermal=mode.IsoThermal;
        else
            IsoThermal=0;
        end
        if IsoThermal==0
            if isfield(mode,'quasi1D') && mode.quasi1D==1
                [DeltaTold,PTherm,T_Contributi,fattore_correttivo,condzTe]=f_ThermD1ANA(mesh,mode,StrTT,IPLOT);
                mode.fattore_correttivo(indv)=fattore_correttivo;
                mode.condzTe(indv,:)=condzTe;
                Tprec=0;
                
                StimaTempWU
                DeltaT=DeltaTold*Fat_STIMA_Temp;
                
                mode.FatMob=polyval(cot,max(max(DeltaT)));
                mode.cot=cot; % coefficients used in 1D simulations to fit the mob change with T, see fitmob.m
                mode.FatMob=1;
            else
% 			save TERM Tprec mesh mode StrTT IPLOT % save intermediate input of thermal solver for debugging
%                 [DeltaTold,Tprec,PTherm,T_Contributi]=f_ThermicFun(Tprec,mesh,mode,StrTT,IPLOT);
                if mode.ThermalFake==0
                   [DeltaTold,Tprec,PTherm,T_Contributi,fattore_correttivo]=f_ThermicFunMOD(Tprec,mesh,mode,StrTT,IPLOT);
                else
                   disp('Call the fake thermal solver')
                   [DeltaTold,Tprec,PTherm,T_Contributi,fattore_correttivo]=FakeThermalSolver(Tprec,mesh,mode,StrTT,IPLOT) ; 
               
                end                
%                  mode.FatMob=1;
                
                StimaTempWU
                DeltaT=DeltaTold*Fat_STIMA_Temp;
                mode.fattore_correttivo(indv)=fattore_correttivo; 
            end
        else
            if mode.quasi1D==1
                mesh.DeltaT=0;
            else
                DeltaT=zeros(1,mesh.nn);
            end
            PTherm=0;
            for tcj=1:5
                T_Contributi{tcj}=0;
            end
        end
        
        %[DeltaT,Tprec,PTherm,T_Contributi]=f_ThermicFun(Tprec,mesh,mode,StrTT,IPLOT);
        DeltaTJoule=T_Contributi{1};
        DeltaTRec_srhAu=T_Contributi{2};
        DeltaTRec_Ccap=T_Contributi{3};
        DeltaTRec_RAD=T_Contributi{4};
        DeltaTOptAbs=T_Contributi{5};
        DeltaT=C_Temp*reshape(DeltaT,1,mesh.nn);
    elseif(flagCallVELM==1)
        DeltaT=zeros(1,mesh.nn);
        DeltaTJoule=DeltaT; DeltaTOptAbs=DeltaT;
        DeltaTRec_srhAu=DeltaT; DeltaTRec_Ccap=DeltaT;
        DeltaTRec_RAD=DeltaT;
        Tprec=0; % initalization of Tprec=0 for thermic self-consistency
        
    end
    
%     if isfield(mode,'epsNLg')
%         nlG=mode.epsNLg;
%         mode.nlG=1;
%         if isfield(mode,'E2')
%             %       'qui per nlG', keyboard
%             UU = (mode.Pst.*mode.fPdif)*mode.E2;
%             uudu=1./(1+mode.epsNLg*UU);
%             mode.nlG=uudu(1:nnQW);
%         end
%     else
%         mode.nlG=1;
%     end
    
    
    mesh.DeltaT=DeltaT;
    
    
    if effetti(6)==0
        mesh.DeltaTvelm=DeltaT;
    else
        'TEMPERATURA Fissa in VELM'
    end
    
    DeltaTmax=max(DeltaT);
    if DeltaTmax-mode.DeltaTmax(end)>DTM0/4
        'incremento di temperatura troppo grande, esco '
        break
    end
    DeltaTmax_Joule=max(max(DeltaTJoule));
    DeltaTmax_srhAu=max(max(DeltaTRec_srhAu));
    DeltaTmax_Ccap=max(max(DeltaTRec_Ccap));
    DeltaTmax_RAD=max(max(DeltaTRec_RAD));
    DeltaTmax_OptAbs=max(max(DeltaTOptAbs));
    
    % giusto per non buttare via quel plot, faccio così.
    if(mode.qflg)
        n2Dtmp=0; p2Dtmp=0; n3Dtmp=0; p3Dtmp=0;
        for indQW=1:NQW
            n2Dtmp=n2Dtmp+mode.n2D{indQW};
            p2Dtmp=p2Dtmp+mode.p2D{indQW};
            n3Dtmp=n3Dtmp+mode.elec(mesh.inMQW{indQW});
            p3Dtmp=p3Dtmp+mode.hole(mesh.inMQW{indQW});
        end
        
        npMax=max(n2Dtmp+p2Dtmp)/2;
        nMax=max(n2Dtmp);
        pMax=max(p2Dtmp);
        n3Max=max(n3Dtmp);
        p3Max=max(p3Dtmp);
        %     nMax=max(mode.n2D);
        %     pMax=max(mode.p2D);
        %     n3Max=max(mode.elec(mesh.inQW));
        %     p3Max=max(mode.hole(mesh.inQW));
        mode.npMaxVet(indv)=npMax;
        mode.nMaxVet(indv)=nMax;
        mode.pMaxVet(indv)=pMax;
        mode.n3MaxVet(indv)=n3Max;
        mode.p3MaxVet(indv)=p3Max;
        
        if iold==2 && mode.Oflg==1
            % calcolo guadagno e derivate
            %    'qui derivate', keyboard
            NUOVA_Gain
            % carico dati tau di Vallone
            global CO coMh Tn Th
            if mode.taucarrier==1
                load saTaue
            end
        end
    end
    
    load stopDD
    if istop==1
        'pausa simulazione', keyboard
    end
    
    mode.rAperture=ParVet(1);   % needed to define the TJ radius
    
    if mode.qflg==1 && mode.Oflg==1
        %     flagPort=(flagCallVELM.*abs(DeltaTmax-DeltaTempVELMLast)>=mode.mintempVELM);
        DELTAT=abs(DeltaTmax-DeltaTempVELMLast)
        flagTemp=(DELTAT>=mode.mintempVELM);
        flagPort=(abs(npMax-npMaxLast)>=mode.minPorVELM);
        flagPot=(Pst>=.1);
        %        flagPot=(Pst*1000>=.1);
        %        if isfield(mode,'Scheck')
        %          flagTemp=(flagTemp &        max(full(mode.Scheck))<.8)
        %        end
        %flagTemp=0
        %        flagPot=0
        %       if isfield(mode,'Gmod')
        %        iGmo=length(find(mode.Gmod>0));
        %        if iGmo>2 & flagPot
        %         iVELM=1;
        %        end
        %       end
        %        if( (indv==1 | flagTemp | flagPort | flagPot | iVELM) & flagCallVELM==1)
        
        flagPot=0;
        FLAG_velm0=( indv==1 | (flagCallVELM & flagTemp) | iVELM | flagPot) & HARD;
        flagvv=0;
        if indv>1
            if ~isempty(VELMInfo)
                if mode.IdriveON==0
                    vlast=abs(mode.v0_dd(VELMInfo(indVELM).indVoltage));
                    flagvv=vlast>=abs(v0_dd);
                else
                    ilast=abs(mode.i0_dd(VELMInfo(indVELM).indVoltage));
                    flagvv=ilast>=abs(i0_dd);
                end
            end
        end
        
        FLAG_velm=(((FLAG_velm0 & ~flagvv & flagCallVELM) & HARD) | RunnaVELM)
         if RunnaVELM==1
%           'vede FLAG velm',               keyboard
          RunnaVELM=0; 
		 end 

%        keyboard
%        keyboard
        if FLAG_velm
            fprintf(' Calcolo VELM !!!!!! \n Calcolo VELM !!!!!! \n Calcolo VELM !!!!!! \n')
            
            if ivelMstop==1
                keyboard
            end
            
            %         iVELM=0;
            
            DeltaTempVELMLast=DeltaTmax;
            
            
            if effetti(7)==1
                VelmOptions.ianti_gui=0;  % tolgo anti_guiding
            end
            
            colordef black
            %                        'prima di Velm', keyboard
            iLoadVELM=0; % 1, check if VELM results are already computed; 0, always enters in VELM
            %
            % save SAVELM
            %
            if (iLoadVELM==1 && exist(['dati/VELM/velm_',strSave,'_quasi1D=',num2str(mode.quasi1D),'.mat'])==2 && v0_dd==0)
                % File exists.
                load(['velm_',strSave,'_quasi1D=',num2str(mode.quasi1D),'.mat'])
            else
                verbVELM=mode.verbVELM;
                % File does not exist.
                if isfield(mode,'quasi1D') && mode.quasi1D==1 
                    if mode.Oflg==1
                        [velm] = f_CallVELM_1D(mesh,mode,mode1,ParVet,VelmOptions,fil_str,mode.quasi1D);
                    end
                else
                    if v0_dd==0
                        save(['dati/VELM/VELM_',strSave],'mesh','mode','mode1','ParVet','VelmOptions','fil_str')
                    end
                    %if flgGEOM==0 | indv>1
                    [velm] = f_CallVELM(mesh,mode,mode1,ParVet,VelmOptions,fil_str);
                    %end
                    if isfield(velm,'UNOD')
                        VelmOptions.UNOD=velm.UNOD;
                    end
                end
                %
                if verbVELM<0
                    mode.verbVELM=0;
                end
            end
            
            colordef white
            if isfield(mode,'quasi1D') && mode.quasi1D==1
                velm.E2=ones(mode.nmodes,mesh.nnx)/(mode.AreaOx);
                velm.E2=ones(mode.nmodes,mesh.nnx)/(mode.AreaOpt);
                velm.SW=0;
            end
            %            mode.Lm=velm.Lm(1:mode.nmodes); % prendo solo il primo modo (per ora)
            mode.Lm=velm.Lm;%*1.31; % prendo solo il primo modo (per ora)
            mode.vlambda=velm.vlambda;
            mode.NQW=velm.NQW; % prendo solo il primo modo (per ora)
            mode.Gamma_z=velm.Gamma_z;
            mode.nindexQW=velm.nindexQW;
            mode.fPES=velm.fPES*mode.Cpot;
            vph=Clight./mode.nindexQW;
            mode.fPdif=TAROCCO*velm.fPdif/vph/mode.Cpot*mode.fPdifScaling;
            % mode.fPdif=velm.fPdif/2.3; % ATTENZIONE CHE HO MESSO UN /2 PER FAR TORNARE LE COSEEEEEEEEEEEEEEEEEEEEE
            mode.E2=velm.E2; % original
            
            % % setting non-zero field to avoid problems
            % indZero=find(mode.E2(:,1)<1e-5);
            % mode.E2(indZero,1)=10;
            % mode.E2=mode.E2+10;
            
            mode.Inviluppo_SW=velm.Inviluppo_SW;
            
            indVELM=indVELM+1;
			if indVELM==1
                clear mode1
            end
            
            % Torrelli Valerio: modification 12/12/2022
            % Saving the far field results
            if mode.quasi1D==0
                VELMInfo(indVELM).FF_Intensity = velm.If(:,1); % far field intensity (normalized to 1), this should change in current!
                VELMInfo(indVELM).FF_ApertureAngle = velm.theta; % semi-aperture angle
                VELMInfo(indVELM).FF_angle_at_e2 = velm.ffL; % semi-aperture angle corresponding to 1/e^2 of the peak intensity
            end           
            
            %            VELMInfo(indVELM).Lm=velm.Lm(1:mode.nmodes); % prendo solo il primo modo (per ora)
            VELMInfo(indVELM).Lm=velm.Lm; % prendo solo il primo modo (per ora)
            VELMInfo(indVELM).vlambda=velm.vlambda;
            VELMInfo(indVELM).NQW=velm.NQW; % prendo solo il primo modo (per ora)
            VELMInfo(indVELM).Gamma_z=velm.Gamma_z;
            VELMInfo(indVELM).nindexQW=velm.nindexQW;
            VELMInfo(indVELM).nindexRef=velm.nindexRef;
            VELMInfo(indVELM).fPES=velm.fPES;
            VELMInfo(indVELM).fPdif=velm.fPdif;
            VELMInfo(indVELM).E2=velm.E2;
            VELMInfo(indVELM).Inviluppo_SW=velm.Inviluppo_SW;
            VELMInfo(indVELM).SW=velm.SW;
            
            VELMInfo(indVELM).v0_dd=v0_dd;
            VELMInfo(indVELM).DeltaTmax=DeltaTmax;
            VELMInfo(indVELM).indVoltage=indv;
            
            mode1.LamVelm(indVELM)=velm.vlambda(1);
            %'salve lambda', keyboard
            mode1.TmVelm(indVELM)=DeltaTmax;
            
            mode1.NmVelm(indVELM)=nMax;
            mode1.PmVelm(indVELM)=pMax;
            mode1.IndVelm(indVELM)=indv;
            mode1.LmVelm(indVELM,:)=velm.Lm(1);
            mode1.E2Velm(indVELM,:)=velm.E2(1,:)';
            npMaxLast=npMax;
            %	                        'salva per Velm', keyboard
            salva_per_VELM
            
            flagVELMWasCalled=1;
            
            colordef white
            mode.indv=indv;
        elseif indVELM>=3 && mode.quasi1D==5
            if mode.IdriveON==0
                [mode]=velmFitting(indVELM,fitStart,fitDegree,VELMInfo,mode,v0_dd);
            elseif indVELM>=4
                [mode]=velmFitting(indVELM,fitStart,fitDegree,VELMInfo,mode,i0_dd);
            end
            %             'Interp',keyboard
        end
    end
    
    T=mode.T0+DeltaT;
    mesh.Tgain=T;
    %    'qui ver T', keyboard
    T=mode.T0+DeltaT*mode.C_Temp_DD;
    
    TKT=T;
    if isfield(mode,'KT')
        if mode.KT==0
            TKT=mode.T0*ones(size(T));
        end
    end
    
    Vt=kB.*TKT./qel;
    %    mode.Vt=Vt;
    
    mesh=loadmacro(geom,mesh,mode,T);
    
    mode.Vt=Vt;
    Vt_tr=pdeintrp(mesh.node,mesh.triangle(1:4,:),Vt.'); % T on triangles
    mode.Vt_tr=Vt_tr;
    
    mode_old=mode;
    mesh_old=mesh;
    
    mode=f_ComputeNeutrality(geom,mesh,mode);
    
    if(mode.Tflg==1)
        %        mode.DeltaTmax(indv)=DeltaTmax;
        mesh.Tprec=Tprec;
        
    end
    
    % calcolo carrier-dependent tau Vallone
    if(mode.taucarrier==1)
        indQW=1; % primo QW
        n2D=mode.n2D{indQW}/mesh.WMQW{indQW};
        p2D=mode.p2D{indQW}/mesh.WMQW{indQW};
        n3D=mode.elec(mesh.inMQW{indQW});
        p3D=mode.hole(mesh.inMQW{indQW});
        if(max(n2D) > 1e8)
            mode.Fat_cap_e=fit_tauEnl(n2D,n3D);
            mode.Fat_cap_h=fit_tauH(p2D,p3D);
            %            'VALLONE', keyboard
        end
    end
    
    
    iterMag1=0;
    t0_Newton=1;
    
    fprintf('\n Drift-diffusion analysis\n')
    mode.iconv(indv) = 0; % flag, 1 convergence, 2 convergence failure
    iter = 0; % number of outer iterations
    iter_t = 0; % number of inner iterations
    mode.res(indv) = Inf; % residual error
    % t = 1; % damping parameter % QUAAAAAA
    %
    
    if effetti(1)==1
        mode.fPES=VELMInfo(le_VELMinfo).fPES*1e3;
    end
    if effetti(2)==1
        mode.Gamma_z=VELMInfo(le_VELMinfo).Gamma_z;
    end
    if effetti(3)==1
        %            mode.fPdif=TAROCCO*velm.fPdif/vph/1000;
        mode.fPdif=TAROCCO*VELMInfo(le_VELMinfo).fPdif/vph/mode.Cpot;
    end
    if effetti(4)==1
        mode.E2=VELMInfo(le_VELMinfo).E2;
    end
    if effetti(5)==1
        mode.vlambda=VELMInfo(le_VELMinfo).vlambda;
    end
    fprintf('fPES Gamma_z  fPdif E2 lam\n')
    

    clear Res
    
    % Adiabatic section
    flgAdiab=0;
    if indv>1
        if v_dd<mode.Vadiab
            flgAdiab=0;
            mode.FatAdiab=mode.AdiabVec(1);
        elseif adiabaticCounter>=length(mode.AdiabVec)
            flgAdiab=0;
            mode.FatAdiab=1;
        else
            flgAdiab=1;
            mode.FatAdiab=mode.AdiabVec(adiabaticCounter);
            disp(['Adiabatic mode: ',num2str(adiabaticCounter)',' of ',num2str(length(mode.AdiabVec)),...
                ' (',num2str(mode.FatAdiab,'%.4e'),')'])
            keyboard
        end
    end
    mode.FatAdiab=1;

PrecTol0=10;
clear Res


    while(iter<mode.maxiter) % && iter_t<mode.maxiter_t) % inner while-loop
        
        uvet(pp+1:3*nn)=abs(uvet(pp+1:3*nn));
        uvet(nn+1:pp)=abs(uvet(nn+1:pp));

        if mode.Oflg==1
            uvet(ss+1)=abs(uvet(ss+1));
            uvet(vv:ss)=abs(uvet(vv:ss));
        end
        %
        indCountMultiplicationFactor=1; %%% temporaneo
        
        
        % assem jacobian matrix in time domain
        mode.ind_v0=indv;
        
        if mode.vv_dd>mode.VmaxStable
            mode.Shunt=0;
        end

        % Numerical Jacobian 
        if (mode.NumJac==1 && round(v0_dd,3)==round(mode.Vnum,3))
            mode_num=mode;
            uvet_num=uvet;
            
            save(['WorkSpace_Jacobian_',strSave,'_',nomeSav,'.mat'])
            disp('Numerical - run "Main_99_Jnum.m"!')
            keyboard
        end
        
        % Assem: the Newton scheme matrices are assembled - there is an OLD
        % version "assem_mp" working with multiprecision toolbox 
        mode.indv=indv;
        if mode.mp==0
            if mode.IdriveON==0
%                 Zmat=100;%*pi*(StrTT.Rox*1e-4)^2;
                if Zmat>0
                    fprintf(['Input impedance Rs = ',num2str(Zmat),' Ohm ! \n'])
                end
                mode.Zmat=Zmat;
                [Kmat0,Jmat0,Jmat1,Jmat2,uvet,rvet,mode,tvet]=assem_GBT(geom,mesh,mode,uvet,v0_dd);
            else
                mode.Zmat=mode.Ymat;
                [Kmat0,Jmat0,Jmat1,Jmat2,uvet,rvet,mode,tvet]=assem_GBT(geom,mesh,mode,uvet,i0_dd);
            end
        else
            if mode.IdriveON==0
                [Kmat0,Jmat0,Jmat2,uvet,rvet,mode,tvet]=assem_mp(geom,mesh,mode,uvet,v0_dd);
            else
                [Kmat0,Jmat0,Jmat2,uvet,rvet,mode,tvet]=assem_mp(geom,mesh,mode,uvet,i0_dd);
            end
        end
                         
        if(mode.Oflg)
            if full(mode.Scheck)>1
                iterMag1=iterMag1+1;
                if iterMag1>1
                    %         t0_Newton=.5
                    % break
                    
                end
            end
        end
        %
        % Numerical Jacobian section
        if mode.NumJac==1 && indv==10
            keyboard
            [JmatNum]=NumericalJacobian(uvet,rvet,geom,mesh,v0_dd,mode,tvet3);
        end
        %
        Jmat=Jmat0;
        
        [R,C]=dgsequ(Jmat);
        res = norm(R*rvet);
        
        MAXiter=mode.MAXiterMoveOn;
        
        if mode.MoveON && iter+iter_t>MAXiter && mode.ii_dd(end)*1e3 < mode.ImoveON
            res=0;
        end
%		if ContaBisez==9  && mode.ii_dd(end)*1e3<16
%          
%            RRvet(indv)=res ;
%            res=0 ; 
%            
%            
%        end
        
        %==============================================================================================100
        if(isnan(res))==100
            disp('================================================')
            disp('=         Residual error is NaN! :-(           =')
            disp('================================================')
            
            figure(112)
            subplot(3,1,1)
            grid on,box on
            plot(uvet,'o')
            ylabel('uvet')
            subplot(3,1,2)
            grid on,box on
            semilogy(abs(rvet),'o')
            ylabel('rvet')
            subplot(3,1,3)
            grid on,box on
            semilogy(abs(R*rvet),'o')
            ylabel('R\cdot rvet')
            drawnow
            break
        end
        if(mode.BankRose==0 || res<mode.res(indv)) % uncomment for Bank-Rose
        if(iter==0 && mode.report)
            fprintf('Iteration     Residual\n')
            fprintf(sprintf('%4i%20.4e\n',iter,res))
        end
        if(iter>=1 && mode.report)
            fprintf('Iteration     Residual            t\n')
            fprintf(sprintf('%4i%20.4e%15.2e\n',iter,res,t))
        end
        mode.res(indv)=res; % save residual

        Res(iter+1)=res;
        if mode.Oflg==1
            if max(mode.Scheck)>.2 & iprimoRunna==1
%                 RunnaVELM=1;
                RunnaVELM=0
                iprimoRunna=0;
            end
            if length(Res)>3 & max(mode.Scheck)>1
                PrecTol=abs(1-Res(end-1)/Res(end))*100
                %		 pausak
                if PrecTol<PrecTol0
                    %		'RES', keyboard
                    break
                end
            end
        end


        
        if(res<mode.tolconv)
            mode.iconv(indv)=1; 
            mode.res(indv) = res; 

%             figure(111)
%             set(gcf,'Position',[1.3578e3   9.8600e1   6.5920e2  8.6320e2])
%             subplot(311)
%             grid on,box on
%             plot(uvet,'o')
%             ylabel('uvet')
%             subplot(312)
%             grid on,box on
%             semilogy(abs(rvet),'o')
%             ylabel('rvet')
%             subplot(313)
%             grid on,box on
%             semilogy(abs(R*rvet),'o')
%             ylabel('R\cdot rvet')
%             drawnow

            break
        end
        % % % %      if(res<mode.tolconv*1e4)
        % % % %       if res>Res(iter)
        % % % %       mode.iconv(indv)=1; mode.res(indv) = res; break,
        % % % %      end
        % % % %      end
        
        uvet2 = uvet; % save solution, uncomment for Bank-Rose
        iter_t=0; t=t0_Newton; % restore damping parameter % QUAAAAAA
        %=========================================
        % solve matrix equation
%          R=(SM*C); 
%          C=R*SMinv; 
        
%         Jmat = Rd*R*Jmat*C*Cd;
        Jmat = R*Jmat*C;
%         [Lmat,Umat,pm,qm] = lu(Jmat);
        
%         if mode.C1>1
%             deltau = - (SM*C*Cd)*(qm*(Umat\(Lmat\(pm*(Rd*R*rvet)))));
% %             deltau = - (SM*C*Cd)*(Jmat\(Rd*R*rvet));
%         else
%             deltau = - C*Cd*(qm*(Umat\(Lmat\(pm*(Rd*R*rvet)))));
%         end

%         deltau = - C*(qm*(Umat\(Lmat\(pm*(R*rvet)))));
        
        if mode.mp==1
            deltau = - C*(Jmat\(R*rvet));
%             deltau = double(deltau);
        else
            [Lmat,Umat,pm,qm] = lu(Jmat);
            deltau = - C*(qm*(Umat\(Lmat\(pm*(R*rvet)))));
        end

        iter = iter + 1;
        uvet = uvet + t*deltau.'; % update
        %============================================
        elseif mode.BankRose==1 % uncomment for Bank-Rose
            t=t/2; % uncomment for Bank-Rose
            iter_t = iter_t + 1; % uncomment for Bank-Rose
            uvet = uvet2 + t*deltau.'; % update, uncomment for Bank-Rose
        end % uncomment for Bank-Rose
        tmpVar
    end % outer while-loop
    
    if(not(mode.iconv(indv)) && not(isnan(res)))
        disp('================================================')
        disp('=        Convergence failure!!!! :-(           =')
        disp('================================================')
        
        %        figure(112)
        %        subplot(211)
        %        grid on,box on
        %        plot(uvet,'o')
        %       subplot(212)
        %        grid on,box on
        %        plot(rvet,'o')
        %        drawnow
        
        if mode.Oflg==1
            fprintf('parte cancellata\n')
            Scheck=mode.gmod./mode.lmod';
            indWrongScheck=find(Scheck>=1);
            if ~isempty(indWrongScheck)
                disp('================================================')
                disp('Too much code breaks due to Scheck > 1 condition')
                disp(['Restarting code with Pst=',num2str(mode.ScheckMultiplicationFactor),'*Pst'])
                disp('================================================')
                
                uvet(ss+[indWrongScheck])=mode.ScheckMultiplicationFactor*uvet(ss+[indWrongScheck]);
                Pst=uvet(ss+[indWrongScheck]);
            end
        end
        %keyboard
    end
    if ((not(mode.iconv(indv)) || isnan(res)) && ContaBisez<ContaBisezMax*length(mode.ScheckMultiplicationFactor))
        ContaBisez=ContaBisez+1
        uvet=uvet_old;
        mode=mode_old;
        mesh=mesh_old;
        flagConv=0;
        flagCallVELM=0 ;
        if(mode.Oflg)
            if(flagVELMWasCalled==1)
                
                VELMInfo(indVELM)=[];
                flagVELMWasCalled=0;
                indVELM=indVELM-1 ;
            end
            flagCallVELM=0;
            
        end
        
        indv=indv-1;
        if mode.IdriveON==1
            if(indCountMultiplicationFactor==1)
                Ib1=mode.i0_dd(1:indv);
                IbNew=(mode.i0_dd(indv)+mode.i0_dd(indv+1))/2;
                Ib2=mode.i0_dd(indv+1:end);
                % aggiunge tensione
                mode.i0_dd=[Ib1,IbNew,Ib2];
                NVbias=NVbias+1;
            elseif(indCountMultiplicationFactor>=length(mode.ScheckMultiplicationFactor))
                indCountMultiplicationFactor=1;
                Ib1=mode.i0_dd(1:indv);
                IbNew=(mode.i0_dd(indv)+mode.i0_dd(indv+1))/2;
                Ib2=mode.v0_dd(indv+1:end);
                % aggiunge tensione
                mode.i0_dd=[Ib1,IbNew,Ib2];
                NVbias=NVbias+1;
            else
                disp('che gli facciamo fare qui ???')
            end
        else
            if(indCountMultiplicationFactor==1)
                Vb1=mode.v0_dd(1:indv);
                VbNew=(mode.v0_dd(indv)+mode.v0_dd(indv+1))/2;
                Vb2=mode.v0_dd(indv+1:end);
                % aggiunge tensione
                mode.v0_dd=[Vb1,VbNew,Vb2];
                NVbias=NVbias+1;
            elseif(indCountMultiplicationFactor>=length(mode.ScheckMultiplicationFactor))
                indCountMultiplicationFactor=1;
                Vb1=mode.v0_dd(1:indv);
                VbNew=(mode.v0_dd(indv)+mode.v0_dd(indv+1))/2;
                Vb2=mode.v0_dd(indv+1:end);
                % aggiunge tensione
                mode.v0_dd=[Vb1,VbNew,Vb2];
                NVbias=NVbias+1;
            else
                disp('che gli facciamo fare qui ???')
            end
        end
    else
        flagConv=1;
        flagCallVELM=1;
        ContaBisez=0;
    end
    if (ContaBisez>=ContaBisezMax)
        disp('ContaBisez')
        break
    end
    %
    % =====================================================================
    % save DD solution
    % =====================================================================
    if(flagConv==1) % if residual is not NaN
        %         mode.tolconv=mode.tolconv*2;
        if mode.mp==1
            uvet = double(uvet);
        end
        mode.phi=uvet(1:nn);
        iq = mesh.iq;
        mode.elec=abs(uvet(nn+(1:nn)));
        if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi')))
            tmp = mode.elec(iq)./mesh.Nc(iq);
            nferinv = invferdr(tmp,mode.tolinvferdr);
            mode.EFn(iq) = mode.ecb(iq) + Vt(iq).*nferinv; 
        end
        %
        mode.hole=abs(uvet(pp+(1:nn)));
        if((isfield(mode,'stats'))&&(strcmp(mode.stats,'Fermi')))
            tmp = mode.hole(iq)./mesh.Nv(iq);
            nferinv = invferdr(tmp,mode.tolinvferdr);
            mode.EFp(iq) = mode.evb(iq) - Vt(iq).*nferinv; 
        end
        %
        fprintf('Elapsed time: %g, s\n',(cputime-elapsed_time))
        mode.elapsed_time=(cputime-elapsed_time)/60;
        
        v_dd=uvet(qq+1); % saving voltage
        i_dd=uvet(rr+1)*mode.CarrierNorm; % saving current
        
        if(mode.Oflg==0)
            mode.ii_dd(indv)=i_dd;
            mode.vv_dd(indv)=v_dd;
            mode.vv0_dd(indv)=v0_dd;
            % compute electric field (in time domain), V/um
            efieldx_t=-mode.phi*mesh.gradx*1e-4; efieldx_t(it)=0;
            efieldy_t=-mode.phi*mesh.grady*1e-4; efieldy_t(it)=0;
            mode.efield_x=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),efieldx_t),1,mesh.nn);
            mode.efield_y=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),efieldy_t),1,mesh.nn);
            if(mode.qflg==1)
                mode.nQW{indv}=mode.n2D;
                mode.pQW{indv}=mode.p2D;
                % da togliere
                mode.IntCcapN(indv,:)=mode.IntCcapn;
                mode.IntCcapP(indv,:)=mode.IntCcapp;
                mode.IntRecN(indv,:)=mode.IntRecN;
                mode.IntRecP(indv,:)=mode.IntRecP;
                %               mode.IntCcapN(indv,:)=mode.IntCcapn;
                %               mode.IntCcapP(indv,:)=mode.IntCcapp;
                mode.IntREc(indv,:)=mode.IntRec;

                % calcolo leakage
                leakage
            end
            if(mode.Tflg==1)
                mode=f_EvalHeatingTerms(geom,mesh,mode);
            end
        else
            Pst=uvet(ss+[1:mode.nmodes]);
            if(max(mode.Scheck)>=1)
                Scheck=mode.Scheck;
                mode=mode_old;
                mesh=mesh_old;
                uvet=uvet_old;
                if(flagVELMWasCalled==1)
                   VELMInfo(indVELM).indVoltage=indv-1;
                    flagVELMWasCalled=0;
                end
                indv=indv-1;
                flagConv=0;
                % Vb1=mode.v0_dd(1:indv);
                % VbNew=(mode.v0_dd(indv)+mode.v0_dd(indv+1))/2;
                % VbNew=v0_dd;
                % Vb2=mode.v0_dd(indv+1:end);
                % aggiunge tensione
                % mode.v0_dd=[Vb1,VbNew,Vb2];
                % NVbias=NVbias+1;
                flagCallVELM=0;
                
                %%%%%%%%
                
                % indv=indv-1;
                % Vb1=mode.v0_dd(1:indv);
                % VbNew=(mode.v0_dd(indv)+mode.v0_dd(indv+1))/2;
                % Vb2=mode.v0_dd(indv+1:end);
                % mode.v0_dd=[Vb1,VbNew*[1 1],Vb2];
                % NVbias=NVbias+1;
                % flagCallVELM=0;
                % countScheck=countScheck+1;
                %%%%%%%
                
                countScheck=countScheck+1;
                if(countScheck>mode.maxScheckRepeat)
                    indWrongScheck=find(Scheck>=1);
                    disp('================================================')
                    disp('Too much code breaks due to Scheck > 1 condition')
                    disp(['Restarting code with Pst=',num2str(mode.ScheckMultiplicationFactor(indCountMultiplicationFactor)),'*Pst'])
                    disp('================================================')
                    uvet(ss+[indWrongScheck])=mode.ScheckMultiplicationFactor(indCountMultiplicationFactor)*uvet(ss+[indWrongScheck]);
                    Pst=uvet(ss+[indWrongScheck]);
                    indCountMultiplicationFactor=indCountMultiplicationFactor+1;
                    %                     keyboard
                end
            else
                countScheck=0;
                flagCallVELM=1;
                mode.vv0_dd(indv)=v0_dd;
                mode.vv_dd(indv)=v_dd;
                mode.ii_dd(indv)=i_dd;
                I_mA=i_dd*1000;
                Ptotal=sum(Pst(:));
                PtotV(indv)=Ptotal;
                if indv>10
                    diPot=diff(PtotV);
                    if diPot(end)<0
                        if Ptotal<PotMin & max(PtotV)>PotMin
                            CondPotBreak=1;
%                             keyboard
                        end
                    end
                end
                if mode.Oflg==1
                    mode.Pst_dd(:,indv)=Pst(:)/mode.Cpot*1e3;
                    mode.Psp_dd(:,indv)=mode.Psp(end)/mode.Cpot*1e3;
                    mode.PspBulk_dd(:,indv)=mode.PspBulk(end);
                    mode.lambda(:,indv)=mode.vlambda;
                else
                    mode.Pst_dd(:,indv)=0;
                    mode.Psp_dd(:,indv)=0;
                    mode.PspBulk_dd(:,indv)=0;
                    mode.lambda(:,indv)=0;
                end
                
                mode.DeltaTmax(indv)=DeltaTmax;
                mode.DeltaTmax_Joule(indv)=DeltaTmax_Joule;
                mode.DeltaTmax_srhAu(indv)=DeltaTmax_srhAu;
                mode.DeltaTmax_Ccap(indv)=DeltaTmax_Ccap;
                mode.DeltaTmax_RAD(indv)=DeltaTmax_RAD;
                mode.DeltaTmax_OptAbs(indv)=DeltaTmax_OptAbs;
                
                mode.PTherm(:,indv)=PTherm;

                mode.Gmod(:,indv)=full(mode.gmod);
                mode.Lmod(:,indv)=mode.lmod;
                mode.nQW{indv}=mode.n2D;
                mode.pQW{indv}=mode.p2D;
                                % da togliere
                mode.IntCcapN(indv,:)=mode.IntCcapn;
                mode.IntCcapP(indv,:)=mode.IntCcapp;
                mode.IntRecN(indv,:)=mode.IntRecN;
                mode.IntRecP(indv,:)=mode.IntRecP;
                %               mode.IntCcapN(indv,:)=mode.IntCcapn;
                %               mode.IntCcapP(indv,:)=mode.IntCcapp;
                mode.IntREc(indv,:)=mode.IntRec;
                
                % calcolo leakage
                leakage
                
                % compute electric field (in time domain), V/um
                efieldx_t=-mode.phi*mesh.gradx*1e-4; efieldx_t(it)=0;
                efieldy_t=-mode.phi*mesh.grady*1e-4; efieldy_t(it)=0;
                mode.efield_x=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),efieldx_t),1,mesh.nn);
                mode.efield_y=reshape(pdeprtni(mesh.node,mesh.triangle(1:4,:),efieldy_t),1,mesh.nn);
                mode=f_EvalHeatingTerms(geom,mesh,mode);
            end
        end
    end
    
    if IPLOT==100 && mode.qflg==1 && flagConv==1
        figure(portfig)
        len=1:length(mode.ii_dd);
        plot(1000*mode.ii_dd,1e-12*mode.nMaxVet(len)/mesh.NMQW,'o-',...
            1000*mode.ii_dd,1e-12*mode.pMaxVet(len)/mesh.NMQW,'+-')
        hold on
        plot(1000*mode.ii_dd,1e-18*mode.n3MaxVet(len),'--',...
            1000*mode.ii_dd,1e-18*mode.p3MaxVet(len),'--')
        xlabel('Current (mA)')
        ylabel('N-P density (PER well) 1e12/cm^2)')
    end
    %            figure(portfig)
    %            plot(1000*mode.ii_dd,1e-12*mode.nMaxVet/3,'o-',...
    %            1000*mode.ii_dd,1e-12*mode.pMaxVet/3,'+-')
    %            xlabel('Current (mA)')
    %            ylabel('N-P density 1e12/cm^2)')
    %       figure
    %               plot(1000*mode.ii_dd,1e-12*mode.nMaxVet(1:end-1),'o-',...
    %           1000*mode.ii_dd,1e-12*mode.pMaxVet(1:end-1),'+-')
    
    if(mode.Oflg)
        vind=[];
        for indplot=1:length(VELMInfo)
            vind=[vind,VELMInfo(indplot).indVoltage];
        end
        if IPLOT==1 && strcmp(mode.ContactPosition,'right')
            T0=mode.T0-273;
                
            figure(1234+kfig),clf
			% This figure contains IV, Rdd, LI, lambda(I)
            set(gcf,'Position',[ 50   150   800   800])
            
            subplot(2,2,1)
            hold on
            grid on
            box on
%             plot(Vmeas,Imeas,'r','LineWidth',2)
            
            % Torrelli: modification 16/11/2022
            % Adding the experimental data
            if strfind(mode.strName,'Stephan')
                load('Stephan_experimental_data')
                IndexT = find(data.LIV.T == T0);
                if isempty(IndexT)
                    error('Insert a valid T for experimental results!')
                end
                Vmeas=data.LIV.V{IndexT};
                Imeas=data.LIV.I{IndexT};
                Rmeas=data.LIV.Rdiff{IndexT};
                Lmeas=data.LIV.P{IndexT};
                
                Cur=ESFit.I;
                P_dB=10*log10(squeeze(ESFit.ModalP(:,IndexT,:)));
                LAM=squeeze(ESFit.lambda_highest(:,IndexT,:));
			else
				% Markus data, for example!
				load(['MarkusN_',num2str(Isize),'_T',num2str(double(T0)),'.mat'])
                pa=1:10:length(Imeas);
                
                Rmeas=diff(Vmeas)./diff(Imeas)*1e3;
                Rmeas=[NaN; Rmeas];

            end
            plot(Vmeas,Imeas,'-or','LineWidth',1.5,'MarkerSize',2)


            plot(mode.vv_dd,mode.ii_dd*1000,'b.','LineWidth',2)
            plot(mode.vv_dd(vind),mode.ii_dd(vind)*1000,'bo','LineWidth',2)
            axis([0 mode.vv_dd(end)+.2 0 mode.ii_dd(end)*1e3+1])
            xlabel('Voltage, V')
            ylabel('Current, mA')
            title(['Vcurrent = ',num2str(v_dd),' V'])
            legend('Measurements','Simulation','Location','Best')
            
            
            subplot(2,2,2)
            hold on
            grid on
            box on
            
            if length(mode.vv_dd)>1
                Rdiff = [NaN, diff(mode.vv_dd)./diff(mode.ii_dd)];
                % Rdiff(2) = NaN;
            else
                Rdiff = NaN;
            end
            
            % Torrelli: modification 16/11/2022
            % Adding the experimental data
            plot(Vmeas,Rmeas,'-or','LineWidth',1.5,'MarkerSize',2)
            plot(mode.vv_dd,Rdiff,'b.','LineWidth',2)
            plot(mode.vv_dd(vind),Rdiff(vind),'bo','LineWidth',2)
            % axis([0 mode.vv_dd(end)+.2 0 mode.ii_dd(end)*1e3+1])
            if mode.vv_dd(end)+.2>1.5
                xlim([1.5, mode.vv_dd(end)+.2])
                ylim([20 350])
            else
                xlim([0, mode.vv_dd(end)+.2])
            end
            xlabel('Voltage, V')
            ylabel('R_{diff}, \Omega')
            title(['Voltage ',num2str(indv),' of ',num2str(NVbias)])
            %        legend('Measurements','Simulation','Location','Best')
            drawnow
            
            
            subplot(2,2,3)
            hold on
            grid on
            box on
            plot(Imeas,Lmeas,'r','LineWidth',2)
            chold
            plot(Cur,10.^(P_dB/10),'s--','LineWidth',2)
            PPst=sum(mode.Pst_dd,1)+mode.Psp_dd;
            
            chold
            plot(mode.ii_dd*1000,mode.Pst_dd,'.','LineWidth',2)
            axis([0 mode.ii_dd(end)*1e3+1 0 max(sum(mode.Pst_dd,1))+.5])
%             legend('Measurements','Simulation','Location','Best')
            plot(mode.ii_dd(vind)*1000,PPst(1,vind),'bo','LineWidth',2)
            
            
            xlabel('Current, mA')
            ylabel('Optical power, mW')
            title(['Heat sink temperature: T=',num2str(T0),'°C'])
            drawnow
            
            subplot(2,2,4)
            grid on
            hold on
            box on
%                 figure(555),clf
%                 set(gcf,'position',[667 64 554 392])
            plot(mode.ii_dd*1e3,mode.lambda,'b.-','linewidth',1.5)
            plot(Cur,LAM,'ro','linewidth',1,'markersize',4)
            xlabel('Current, mA')
            ylabel('Emission wavelength, nm')
            xlim([0 mode.ii_dd(end)*1e3+0.01])
            drawnow
            title('VCSEL thermometer')

            % Carrier density and energy band diagram (central column)
            figure(439),clf
            set(gcf,'position',[667 64 554 392])
            grid on
            hold on
            box on
            plot(mesh.node(2,1:mesh.nny)*1e4,mode.elec(1:mesh.nny),'b','linewidth',1.5)
            plot(mesh.node(2,1:mesh.nny)*1e4,mode.hole(1:mesh.nny),'r','LineWidth',1.5)
            xlabel('z, \mum')
            ylabel('Elec. and holes, cm^{-3}')
            set(gca,'yscale','log')
            set(gca,'FontSize',14)
            xlim([StrTT.Tbuf+StrTT.Tdbr_inf-1 StrTT.Tbuf+StrTT.Tdbr_inf+StrTT.Tcav+1]),hold off
            drawnow
            title(['Red/blue: holes/electrons'])
            
            % Torrelli Valerio: modification 12/12/2022
            % Adding FF results
            
            % Saving the aperture in mode as well
            if mode.quasi1D==0
                mode.FF_angle_at_e2(indv) = velm.ffL;
                
                FFIntensity = velm.If(:,1);
                %             FFIntensity = VELMInfo(indVELM).FF_Intensity;
                % FFIndex = find(FFIntensity>=exp(-2)*max(FFIntensity));
                FFIndex = find(FFIntensity>=0.25*max(FFIntensity));
                FFIndex = max(FFIndex);
                mode.FF_angle_manual(indv) = velm.theta(FFIndex);
                %             mode.FF_angle_manual(indv) = VELMInfo(indVELM).FF_ApertureAngle(FFIndex);
                
                figure(6000),clf
                set(gcf,'position',[667 64 554 392])
                hold on
                % plot(mode.ii_dd*1e3,mode.FF_angle_at_e2,'-ob','linewidth',1.5,'MarkerSize',2)
                % plot(mode.ii_dd(vind)*1e3,mode.FF_angle_at_e2(vind),'ob','linewidth',1.5)
                plot(mode.ii_dd*1e3,mode.FF_angle_manual,'-ob','linewidth',1.5,'MarkerSize',2)
                plot(mode.ii_dd(vind)*1e3,mode.FF_angle_manual(vind),'ob','linewidth',1.5)
                hold off
                grid on
                box on
                xlabel('Current, mA')
                ylabel('\theta/2, °')
                set(gca,'FontSize',12)
                title('FF width @ 25%')
                xlim([0 mode.ii_dd(end)*1e3+0.01])
                drawnow
            end
            

            figure(440),clf
            set(gcf,'position',[1249  572  604  407])
            grid on
            hold on
            box on
            plot(mesh.node(2,1:mesh.nny)*1e4,mode.ecb(1:mesh.nny),'b.-','linewidth',2)
            plot(mesh.node(2,1:mesh.nny)*1e4,mode.evb(1:mesh.nny),'r.-','linewidth',2)
            plot(mesh.node(2,1:mesh.nny)*1e4,mode.EFn(1:mesh.nny),'k-.','linewidth',1)
            plot(mesh.node(2,1:mesh.nny)*1e4,mode.EFp(1:mesh.nny),'k--','linewidth',1)
            xlim([StrTT.Tbuf+StrTT.Tdbr_inf-1 StrTT.Tbuf+StrTT.Tdbr_inf+StrTT.Tcav+1]),hold off
            title(['I = ',num2str(mode.ii_dd(end)*1e3),' mA'])
            xlabel('z, \mum')
            ylabel('Band diagram, eV')
            set(gca,'FontSize',14)
            drawnow
            
            DensityCurrentPlot
            
            if IPLOT==1 && mode.Tflg==1
                % DeltaT (each contribution)
                figure(1235+kfig),clf
                set(gcf,'Position',[  634    73   579   365])
                hold on
                grid on
                box on
                plot(mode.ii_dd*1e3,mode.DeltaTmax,'k','LineWidth',2)
                plot(mode.ii_dd*1e3,mode.DeltaTmax_Joule,'c--','LineWidth',2)
                plot(mode.ii_dd*1e3,mode.DeltaTmax_OptAbs,'g--','LineWidth',2)
                plot(mode.ii_dd*1e3,mode.DeltaTmax_Ccap,'r--','LineWidth',2)
                plot(mode.ii_dd*1e3,mode.DeltaTmax_RAD,'b--','LineWidth',2)
                plot(mode.ii_dd*1e3,mode.DeltaTmax_srhAu,'y--','LineWidth',2)
                plot(mode.ii_dd(vind)*1e3,mode.DeltaTmax(vind),'ko','LineWidth',2)
                legend('\DeltaT','Joule','Opt. abs.','Ccap','Rad','SRH/Aug')
                xlabel('Current, mA')
                ylabel('Temperature rise, K')
            end
            
            % Gain vs Losses
            figure(120)
            hold on
            grid on
            box on
            plot(mode.ii_dd*1e3,mode.Gmod,'k','LineWidth',2)
            plot(mode.ii_dd*1e3,mode.Lmod,'r--','LineWidth',2)
            if(mode.nmodes>1)
                plot(mode.ii_dd(vind)*1e3,mode.Gmod(:,vind)','ko','LineWidth',2)
                plot(mode.ii_dd(vind)*1e3,mode.Lmod(:,vind)','ro','LineWidth',2)
            else
                plot(mode.ii_dd(vind)*1e3,mode.Gmod(vind),'ko','LineWidth',2)
                plot(mode.ii_dd(vind)*1e3,mode.Lmod(vind),'ro','LineWidth',2)
            end
            xlabel('Current, mA')
            ylabel('Modal gain vs losses, cm^{-1}')
            
            gipos=find(mode.Gmod(1,:)>0);
            if length(gipos)>1
                xli=mode.ii_dd(gipos([1 end]))*1e3 +[-.1,.1];
                axis([xli 0 max(max(mode.Lmod))] )
            end
            %        legend('Gain','Losses','Location','Best')
            %keyboard
            drawnow
                        
            if(indv>1 & IPLOT==100)
                figure(1236+kfig),clf
                set(gcf,'Position',[631   547   585   365])
                subplot(1,2,1)
                hold on
                grid on
                box on
                Resistence=diff(mode.vv_dd)./diff(mode.ii_dd);
                Resistence=[Resistence(1),Resistence];
                % plot(mode.ii_dd*1000,Resistence,Imeasinterp,Resmeasinterp)
                plot(Imeasinterp,Resmeasinterp,'r-','LineWidth',2)
                plot(Imeas_res,Rmeas,'.')
                lres=1:length(Resistence);
                plot(mode.ii_dd(lres)*1000,Resistence(lres),'go','LineWidth',2)
                ylim([40,120])
                %            xlim([0.5,1000*mode.ii_dd(end)+.1])
                xlabel('Current, mA')
                ylabel('Differential resistance, \Omega')
                %            legend('Measurements','Simulation','Location','Best')
                subplot(1,2,2)
                hold on
                grid on
                box on
                plot(Cur,LAM,'ro-','LineWidth',2)
                plot(mode.ii_dd*1000,mode.lambda,'b.','LineWidth',2)
                %            legend('Measurements','Simulation','Location','Best')
                if(mode.nmodes>1)
                    plot(mode.ii_dd(vind)*1000,mode.lambda(:,vind),'bo','LineWidth',2)
                else
                    plot(mode.ii_dd(vind)*1000,mode.lambda(vind),'bo','LineWidth',2)
                end
                %            xlim([0.,1000*mode.ii_dd(end)+.1])
                xlabel('Current, mA')
                ylabel('Modal wavelength, nm')
                drawnow
            end
        end
    else
        if IPLOT==1
            figure(1234+kfig),clf
            set(gcf,'Position',[268 533 1096 420])
            %             subplot(1,2,1)
            hold on
            grid on
            box on
            %             plot(Vmeas,Imeas,'r','LineWidth',2)
            plot(mode.vv_dd,mode.ii_dd*1000,'b.','LineWidth',2)
            xlabel('Voltage, V')
            ylabel('Current, mA')
            title(['Vcurrent = ',num2str(v0_dd),' V'])
            drawnow
%             keyboard
            %        legend('Measurements','Simulation','Location','Best')
            %             subplot(1,2,2)
            %             if(indv>1)
            %                 hold on
            %                 grid on
            %                 box on
            %                 Resistence=diff(mode.vv_dd)./diff(mode.ii_dd);
            %                 Resistence=[Resistence(1),Resistence];
            %                 plot(Imeasinterp,Resmeasinterp,'r-','LineWidth',2)
            %                 plot(mode.ii_dd*1000,Resistence,'b.','LineWidth',2)
            %                 ylim([75,120])
            %                 xlabel('Current, mA')
            %                 ylabel('Differential resistance, \Omega')
            %                 %            legend('Measurements','Simulation','Location','Best')
            %                 drawnow
            %             end
            
            % Carrier density and energy band diagram (central column)
            figure(439),clf
            set(gcf,'position',[667 64 554 392])
            grid on
            hold on
            box on
            plot(mesh.node(2,1:mesh.nny)*1e4,mode.elec(1:mesh.nny),'b','linewidth',1.5)
            plot(mesh.node(2,1:mesh.nny)*1e4,mode.hole(1:mesh.nny),'r','LineWidth',1.5)
            xlabel('z, \mum')
            ylabel('Elec. and holes, cm^{-3}')
            set(gca,'yscale','log')
            set(gca,'FontSize',14)
            xlim([StrTT.Tbuf+StrTT.Tdbr_inf-1 StrTT.Tbuf+StrTT.Tdbr_inf+StrTT.Tcav+1]),hold off
            drawnow
            title(['Red/blue: holes/electrons'])
            

            figure(440),clf
            set(gcf,'position',[1249  572  604  407])
            grid on
            hold on
            box on
            plot(mesh.node(2,1:mesh.nny)*1e4,mode.ecb(1:mesh.nny),'b.-','linewidth',2)
            plot(mesh.node(2,1:mesh.nny)*1e4,mode.evb(1:mesh.nny),'r.-','linewidth',2)
            plot(mesh.node(2,1:mesh.nny)*1e4,mode.EFn(1:mesh.nny),'k-.','linewidth',1)
            plot(mesh.node(2,1:mesh.nny)*1e4,mode.EFp(1:mesh.nny),'k--','linewidth',1)
            xlim([StrTT.Tbuf+StrTT.Tdbr_inf-1 StrTT.Tbuf+StrTT.Tdbr_inf+StrTT.Tcav+1]),hold off
            title(['I = ',num2str(mode.ii_dd(end)*1e3),' mA'])
            xlabel('z, \mum')
            ylabel('Band diagram, eV')
            set(gca,'FontSize',14)
            drawnow
            
            if mode.v0_dd>0
                DensityCurrentPlot
            end
        end
    end
    
    if(flagConv)
        if mode.mp==0
            savePlotD
        else
            mp_savePlotD
        end
        %
        if mode.Oflg==1 %&& mode.quasi1D==1
            PPst=sum(modePlot.Pst_dd,1)+(modePlot.Psp_dd/mode.frsp + modePlot.PspBulk_dd)*(1-mode.fat_RAD);
            PElec=modePlot.ii_dd*1e3.*modePlot.vv_dd;
            if mode.flgBTJ==1 && mode.flgHeatTJ==2
                % See further comments in the "settings_vari_"
                % mode.VTJ(indexBTJ,mode.ind_v0,:)=Vint; -> in assem_GBT.m
                if mode.quasi1D==0
                    xPerc=40;
                    % Here we want to compute the integral rho drho and divide it by the Area
                    [maxHeat,imax]=max(squeeze(mode.HeatTJ(1,end,:)));
                    avgHeat=maxHeat*xPerc/100;
                    iRagHeatTJ=find(squeeze(mode.HeatTJ(1,end,imax:end))<avgHeat,1)+imax;
%                     iRagHeatTJ=18;
                    
                    rho=mesh.xgrid;
                    drho=diff([0 rho]);
                    xdxN=rho.*drho;
                    
                    Area=pi*rho(iRagHeatTJ)^2;
                    
                    VTJ=squeeze(mode.VTJ(1,end,1:iRagHeatTJ))';
                    Vrho=VTJ*xdxN(1:iRagHeatTJ)'*2*pi/Area;
                else
                    % 1D case
                    VTJ=squeeze(mode.VTJ(1,end,1:2))';
                    Vrho=mean(VTJ);
                end
                PElec=modePlot.ii_dd*1e3.*(modePlot.vv_dd-mode.Fat_VTJ*Vrho);
            end
            PDiss=PElec(end)-PPst(end);
            mode.PDiss(indv)=PDiss;
            
            if (PPst(end)>mode.Pmin_Pfit)
                polGrad=2;
                xstart=-2-polGrad;
                
                if indv<=length(mode.v0_dd) && mode.IdriveON==0
                    xVolt=mode.v0_dd(indv+[xstart:0]);
                    yPdiss=mode.PDiss(indv+[xstart:0]);
                    
%                     coeffPdiss=polyfit(xVolt,log10(yPdiss),1);
                    coeffPdiss=polyfit(xVolt,yPdiss,polGrad);
                    if indv<length(mode.v0_dd)
%                         Pdissfit=10.^polyval(coeffPdiss,mode.v0_dd(indv+1));
                        Pdissfit=polyval(coeffPdiss,mode.v0_dd(indv+1));
                    else
                        Vbias=mode.v0_dd(end)+diff(mode.v0_dd([end-1,end]));
%                         Pdissfit=10.^polyval(coeffPdiss,Vbias);
                        Pdissfit=polyval(coeffPdiss,Vbias);
                    end
                    
                    mode.PDissPred(indv+1)=Pdissfit;
                    
                elseif (indv<=length(mode.i0_dd) && mode.IdriveON==1)
                    xCurr=mode.i0_dd(indv+[xstart:0]);
                    yPdiss=mode.PDiss(indv+[xstart:0]);
                    
%                     coeffPdiss=polyfit(xCurr,log10(yPdiss),1);
                    coeffPdiss=polyfit(xCurr,yPdiss,polGrad);
                    
                    if indv<length(mode.i0_dd)
%                         Pdissfit=10.^polyval(coeffPdiss,mode.i0_dd(indv+1));
                        Pdissfit=polyval(coeffPdiss,mode.i0_dd(indv+1));
                    else
                        Ibias=mode.i0_dd(end)+diff(mode.i0_dd([end-1,end]));
%                         Pdissfit=10.^polyval(coeffPdiss,Ibias);
                        Pdissfit=polyval(coeffPdiss,Ibias);
                    end
                    
                    mode.PDissPred(indv+1)=Pdissfit;
                end
            else
                mode.PDissPred(indv+1)=PDiss;
            end
            
            Tsoglia=1e3;
            if mode.DeltaTmax(indv)>Tsoglia
                
                FatTEMPfit=mode.PDiss(indv-1)./mode.PDiss(indv-2);
                if FatTEMPfit>1.2
                    FatTEMPfit=1.2;
                elseif FatTEMPfit<1.05
                    FatTEMPfit=1.05;
                end
%                 FatTEMPfit=1.10;
                mode.FatTEMP(indv)=FatTEMPfit;
                mode.PDissPred(indv)=mode.PDiss(indv)*FatTEMPfit;
            end
        end
        %
        if mode.Oflg==1 && mode.Pst_dd(end)>0.90 && mode.quasi1D==1
            if isfield(mode,'Pst_max')==1 && mode.Pst_dd(end)<mode.Pst_dd(end-1)
                mode.Pst_max=mode.Pst_dd(end-1);
            end
            if isfield(mode,'Pst_max')==1
                if (mode.Pst_dd(end)/mode.Pst_max)<0.90
                    keyboard, break
                end
            end
        end
        %
        saveModeVET=0;
        if saveModeVET==1 && mode.quasi1D==1
            meshvet{indv}=mesh;
            modevet{indv}=mode;
        else
            meshvet=[];
            modevet=[];
        end
%        'verifica save', keyboard
%         eval(['save ',Last_Workspace,'_Plot',rad_setting,' modePlot'])
        eval(['save ',Last_Workspace,'_Plot_',strSave,' modePlot'])
        MODEplot{kpar}=modePlot;
        eval(['save ',nomeSW,strSave,'_',nomeSav,'_',num2str(IPAR),'.mat MODEplot'])
        % save('Last_plot.mat','modePlot')
        if exist('V0_save')
            if v0_dd==V0_save
                save Riparto
            end
        end
    end
    
    if flgAdiab==0
        indv=indv+1;
    else
        adiabaticCounter=adiabaticCounter+1;
    end
    
    % CURRENT DRIVING SECTION
    if mode.qflg==1 && mode.Oflg==1
        if (mode.Idrive==1 && mode.IdriveON==0 && (sum(Pst)>mode.Pmin || mode.vv_dd(end)>mode.Vmin))
            mode.IdriveON=1;
            if mode.Ilog==1
                Iadd=[(abs(i_dd)+logspace(floor(log10(abs(i_dd))),-4,15))/mode.CarrierNorm mode.Istep:mode.Istep:mode.Imin-mode.Istep mode.Imin:mode.IstepLarge:Imassimo*1e-3/mode.CarrierNorm];
            else
                i_dd=i_dd/mode.CarrierNorm;
                Imax=Imassimo*1e-3/mode.CarrierNorm;
                
                Iadd=[(i_dd+mode.Istep):mode.Istep:(mode.Imin-mode.Istep) mode.Imin:mode.IstepLarge:Imax];
               % Iadd=[(mode.ii_dd(end)/mode.CarrierNorm+mode.Istep): mode.Istep : mode.Imax/mode.CarrierNorm/1e3] ; 
                if mode.Pmin_Pfit>20
                    Iadd=repelem(Iadd,2);
                end
            end
            mode.VIdrive=mode.vv_dd(end);   % last bias step at which Vdrive is used
            disp('Idrive ON')%,keyboard
            
            mode.i0_dd=[mode.ii_dd/mode.CarrierNorm Iadd+mode.vv_dd(end)*mode.Ymat];
            
            NVbias=length(mode.i0_dd);
            
            fitStart=2; % for VELM fitting
        end
    else
        if (mode.Idrive==1 && mode.IdriveON==0 && mode.vv_dd(end)>mode.Vmin)
            mode.IdriveON=1;
            mode.VIdrive=mode.vv_dd(end);   % last bias step at which Vdrive is used
            disp('Idrive ON'),keyboard
            
            i_dd=i_dd/mode.CarrierNorm;
            Imax=Imassimo*1e-3/mode.CarrierNorm;
            
            Iadd=[(i_dd+mode.Istep):mode.Istep:(mode.Imin-mode.Istep) mode.Imin:mode.IstepLarge:Imax];
    
            if mode.Pmin_Pfit>20
                mode.i0_dd=[mode.ii_dd/mode.CarrierNorm repelem(Iadd,2)];
            else
                % Repeated bias vector to avoid any fitting
                mode.i0_dd=[mode.ii_dd/mode.CarrierNorm Iadd];
                pausak
            end
            NVbias=length(mode.i0_dd);
            
            fitStart=2; % for VELM fitting
        end
    end
    
    % Comes back to sp once a certain voltage is reached
    if v_dd>mode.Vmp && mode.mp==1
        disp('Return to sp from mp!'),keyboard
        mode.mp=0;
    end
    
    %  'prima save dentro save', keyboard
    %eval(['save ',nomeSW,'Contributi_',nomeSav,' ',' MODEplot -append '])
    
    %  save Contributi MODEplot -append
    
end


indv=indv-1;

fprintf('Press "Continue" to save the simulation output\n')%,keyboard

if irest==0
%         save sadu MODEplot kpvet kpar IPvet PMAT NP rad_setting  rad_settingV
%         clear MODEplot kpvet kpar IPvet PMAT NP rad_setting  rad_settingV
%         eval(['save ',Last_Workspace])
%         load sadu
    % SAVE MESH AND MODE AS IN D1ANA        
    if ~exist('VELMInput')
        VELMInfo=[];
        VELMInput=[];
    end
    eval(['save ',Last_Workspace,'_',strSave,'_',nomeSav,'.mat geom mesh mode VELMInfo VELMInput MODEplot'])
else
    save ultimo2
end
Resistence=diff(mode.vv_dd)./diff(mode.ii_dd);
ind=5:length(Resistence);
