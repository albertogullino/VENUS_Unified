
function [velm] = f_CallVELM_1D(mesh,mode,mode1,ParVet,VelmOptions,fil_str,i1D)

ioldLoss=0;
if ioldLoss==1
    load saclearglo
    eval(stringaclear)
else
    load saclearglo1
    eval(stringaclear1)
end
%clear global
         
warning off
%==========================================================================
% VELM solver parameters
%==========================================================================
%'in VELM 1D', keyboard

verbVELM=abs(mode.verbVELM);
%verbVELM=2;

if(verbVELM==2)
    ifp_Opt=-10;  % ifp_Opt=-10: modalità verbose, da usare con strutture nuove
else
    ifp_Opt=-4;   % ifp_Opt=-4: modalità veloce
end
if(verbVELM==0)
    Ps.ifpstop=0;
end

fclose('all');
num_azim=VelmOptions.num_azim;
ivett=VelmOptions.ivett;
ipolar=VelmOptions.ipolar;
itutmir=VelmOptions.itutmir;
if isfield(VelmOptions,'imod_acc')==1
    imod_acc=VelmOptions.imod_acc;
else
    imod_acc=0;
end

global E2 ztot nz

if isfield(mode,'Ez')==1
    E2=mode.Ez;
    ztot=mode.ztot;
    nz=mode.nz;
else
    E2=0;
end


if isfield(VelmOptions,'ianti_gui')==1
    ianti_gui=VelmOptions.ianti_gui;
else
    ianti_gui=1;
end

if isfield(VelmOptions,'Pf')==1
    Ps.Pf=VelmOptions.Pf;
end

if isfield(VelmOptions,'dissfun')==1
    Ps.dissfun=VelmOptions.dissfun;
end
if isfield(VelmOptions,'igraef_new')==1
    Ps.igraef_new=VelmOptions.igraef_new;
end
if isfield(VelmOptions,'iBWnew')==1
    Ps.iBWnew=VelmOptions.iBWnew;
end
if isfield(VelmOptions,'iany')==1
    Ps.iany=VelmOptions.iany;
end

Ps.asspar=0;  %come prima
if mode.flgBTJ_lithographic>=1
    Ps.asspar=2;  %TJ vera
end


%'inati', keyboard
Ps.NP_sopra=15;
Ps.NP_sotto=2;

NPx=1e4*mesh.xgrid;

global EO ABS i1D

if i1D==0
    
    dox=ParVet(1);
    
    %'NPx', keyboard
    ired_int=1;    % riduce la zona dei fenomeni trasversali
    Fat_ox=2.5;
    if ired_int==1
        fianti=find(NPx<=Fat_ox*dox);
    else
        fianti=1:length(NPx);
    end
    
    % x=linspace(min(xdd),max(xdd),101);
    % calop.ro_in=x';
    
    if length(mesh.DeltaTvelm)>1
        TddFl=reshape(mesh.DeltaTvelm,mesh.nny,mesh.nnx);
        TddFl=TddFl(:,fianti);
    else
        TddFl=mesh.DeltaTvelm;
    end
    Tdd=flipud(TddFl);
    Tdd=Tdd(:,fianti);  % riduco la zona termica trasversale a dove il campo la sente
    
    
    T0=max(max(Tdd));
    
    % 'stop T0', keyboard
    if T0<mode.DT0/4
        % 'stop T0', keyboard
        %   if ParVet(1)<1
        %   'Ossido troppo piccolo: lo allargo a T=0 !!!!'
        %  ParVet(1)=1;
        % end
    end
    
    
    ro_in=NPx(fianti);
    
    ro_Fie=NPx;
    
    nnxVELM=351;    % original VELM for 2019JSTQE articles
    % nnxVELM=401;
    ro_campo=linspace(0,NPx(fianti(end)),nnxVELM);
    ro_campoF=linspace(0,NPx(end),nnxVELM);
    
    
    
    if length(mode.efield_y)>1
        Es=reshape(mode.efield_y,mesh.nny,mesh.nnx);
    else
        Es=mode.efield_y;
    end
    Es=Es(:,fianti);
    
    
    
    lambda=mean(mode.vlambda)*1e-9;
    
    if length(mode.Deltan)>1
        Deltan=reshape(mode.Deltan,mesh.nny,mesh.nnx);
        alphaDu=reshape(mode.alpha,mesh.nny,mesh.nnx);
        alphaDu=alphaDu(:,fianti);
        alpha_ax=alphaDu(:,1);
        fia=find(alpha_ax>1);
        alphaDue=zeros(size(alpha_ax));
        alphaDue(fia)=alphaDu(fia,end);
        alphaEnd=repmat(alphaDue,1,length(fianti));
        alpha=alphaDu-alphaEnd;
        k0=2*pi/lambda;
        Deltan=Deltan(:,fianti);
        ImagDeltan=alpha./(2.*k0).*(1+Tdd/300).^mode.ABS_Texp;
    else
        Deltan=0;
        ImagDeltan=0;
    end
    %'Deltan', keyboard
    
    
    Ex=reshape(mode.efield_x,mesh.nny,mesh.nnx); % rho field component
    Ey=reshape(mode.efield_y,mesh.nny,mesh.nnx); % z field component
    Ex=Ex(:,fianti);
    Ey=Ey(:,fianti);
    
    % for kz=1:length(calop.zeta)
    %  Tvelm(:,kz)=spline(xdd,Tdd(kz,:),x);
    % end
    %calop.N=epsTf(:,1)/epsTf(1,1);
    
    global inuT inuG
    inuT=1;
    inuG=1;
    
    NPx1=length(NPx);
    iLP=1-ivett;
    % fil_str='jim_ox1.str';
    
    for kk=1:length(ParVet)
        par_in{kk,1}=ParVet(kk);
    end
    % parametro con tag |P1 nel file str
    
    %'callvelm', keyboard
    
    %Dlam=[-0.2 1.3 7 0.3];   %scansione in frequenza per trovare il modo
    %alim=.15;               %kmax; 0.1 va bene per aperture normali (3-4 um)
    %NP_k=20;
    Dlam=VelmOptions.Dlam;   %scansione in frequenza per trovare il modo
    alim=VelmOptions.krel_max;               %kmax; 0.1 va bene per aperture normali (3-4 um)
    NP_k=VelmOptions.NP_k;
    
    % mod_max=[2];   % numero modi
    mod_max=[mode.nmodes];   % numero modi - Tibaldi
    numodiacc=0;        %numero armoniche azimutali accoppiate
    if isfield(VelmOptions,'numodiacc')==1
        numodiacc=VelmOptions.numodiacc;        %numero armoniche azimutali accoppiate
    end
    pasnu=2;             % passo armoniche, sempre 2 per strutture circolari
    
    dk1=alim/NP_k;
    %'qui modi', keyboard
    if ivett==1
        
        if num_azim==2
            Ev_or_Od='Both';   % 'Even' se solo modo fondamentale, 'Odd' se solo modo superiore
        elseif num_azim==1
            Ev_or_Od='Even';   % 'Even' se solo modo fondamentale, 'Odd' se solo modo superiore
        else
            Ev_or_Od=num2str(num_azim-2);
        end
        
    else
        if imod_acc==1
            pasnu=1;
            numodiacc=num_azim-1;
        end
        Ev_or_Od=num2str(num_azim-1);
        
        if num_azim==2
            Ev_or_Od='Both';   % 'Even' se solo modo fondamentale, 'Odd' se solo modo superiore
        elseif num_azim==1
            Ev_or_Od='Even';   % 'Even' se solo modo fondamentale, 'Odd' se solo modo superiore
        else
            Ev_or_Od=num2str(num_azim-2);
        end
    end
    
    %'qui acc', keyboard
    % sistemare
    lambda=lambda*1e6;  % in um
    rrat=3.6;    
    
    
    %'g0', keyboard
    %matDn=mode.matDn; % Tibaldi
    %g0=abs(matgain(end));  % in 1/cm % Tibaldi
    %g0=0;
    
    setVELMparameters %-- chiamare i parametri che non devo settare
    
    if isfield(mode,'DeltaN')==0
        Dn=0;
    else
        Dn=mode.DeltaN;
        %   if length(Dn)>1
        %    Dn=Dn(:,fianti);
        %   end
    end
    
    del_n_ag=min(Dn);
    if del_n_ag==0
        Nref=0;
    else
        del_n_ag0=del_n_ag;
        del_n_ag=min(Dn)-max(Dn);
        Nref=(Dn-Dn(end))/del_n_ag;
        del_n_ag=del_n_ag*mode.fat_ag;
        %del_n_ag=del_n_ag*(1-Nref_du(end));
        %   'cont dfat[', keyboard
        
    end
    
    
    
    matgain=mode.matgain; % Tibaldi
    
    if isfield(VelmOptions,'gain_gui')
        igain=VelmOptions.gain_gui;
        if igain==0
            matgain=0;
        end
    end
    if length(matgain)==1
        N=0;
        %    vph=c_light*100/rrat;
        %    g0=abs(matgain)/vph;
        g0=abs(matgain);
        iga=0;
        if iga==1
            roadd=8;
            ro_in=sort([ro_in roadd+.01]);
            fiN=find(ro_in<=roadd);
            N=zeros(size(ro_in));
            N(fiN)=1;
        else
            g0=0;
        end
        %g0=0;
        %    'quio', keyboard
    else
        matgain=matgain(:,fianti);
        if max(matgain)>-100
            g0=abs(min(matgain));  % in 1/cm % Tibaldi
            N=matgain-min(matgain); N=N./max(N);
            if N(1)<0
                N=-N;
            end
            g0=abs((matgain(end)));  % in 1/cm % Tibaldi
            N=matgain-matgain(end); N=N./max(N);
            if N(1)<0
                N=-N;
            end
        else
            roadd=dox+1;
            N=exp(-(ro_in/roadd).^4);
            g0=0;
            g0=abs(min(matgain));
        end
        %    if(max(matgain)<0)
        %%     g0=0;
        %     velm.Lm=mode.Lm; % prendo solo il primo modo (per ora)
        %     velm.vlambda=mode.vlambda;
        %     velm.NQW=mode.NQW; % prendo solo il primo modo (per ora)
        %     velm.Gamma_z=mode.Gamma_z;
        %     velm.nindexQW=mode.nindexQW;
        %     velm.fPES=mode.fPES;
        %     velm.fPdif=mode.fPdif;
        %     velm.E2=mode.E2; % original
        %    end
    end
    pimra=-g0/cconv_velm;
    calop.pimra=pimra;
    %'g0', keyboard
    calop.zeta=mesh.ygrid*1e4;
    
    
    dndT0=2e-4*mesh.fatt_dndT;
    dndT0=mesh.dndT;

    %     nmed=3.3;
    nmed=3.248; %media di n con |E|^2
    %Ps.dndT0=dndT0*(mode.T0-300);   % valore 1D
    
    T_por=flipud(Deltan-1i.*ImagDeltan);
    T_ter=dndT0*Tdd;
    Tg=2*nmed*(T_ter+T_por);
    %Tg=2*nmed*(T_ter);
    
    zTf=calop.zeta;
    
    dz=mesh.dz_Term/1000;
    Ps.itutmir=itutmir;
    
    Ps.zTemp_start=dz;
    
    zT=(max(zTf)-fliplr(zTf))+dz;
    fiT=find(diff(zT)<.9*max(diff(zT)));
    Ps.dndT0=dndT0*mean(Tdd(fiT,end));   % valore 1D
    zAbs=(max(zTf)-fliplr(zTf));
    calop.zeta=zT;
    
    EO.Es=flipud(Es);   % cambio direzione riferimento.
    EO.zT=zT;
    % 'ver EO', keyboard
    
    ABS.zAbs=zAbs;
    delec=reshape(mesh.dop_d,mesh.nny,mesh.nnx);
    Dope=delec(:,fianti(end))/1e18;
    %Dope=delec(:,1)/1e18;
    dhole=reshape(mesh.dop_a,mesh.nny,mesh.nnx);
    Doph=dhole(:,fianti(end))/1e18;
    %Doph=dhole(:,1)/1e18;
    
    elec=reshape(mode.elecABS,mesh.nny,mesh.nnx);
    %eleccentro=elec(:,end)/1e18;
    eleccentro=elec(:,fianti(end))/1e18*mode.CarrierNorm;
    hole=reshape(mode.holeABS,mesh.nny,mesh.nnx);
    %holecentro=hole(:,end)/1e18;
    holecentro=hole(:,fianti(end))/1e18*mode.CarrierNorm;
    
    xmo=reshape(mesh.xmol,mesh.nny,mesh.nnx);
    xmol=xmo(:,1);
    
    ABS.xmol=xmol;
    ABS.eleccentro=eleccentro;
    ABS.holecentro=holecentro;
    ABS.zT=zT;
    ABS.e=mode.ABSe;
    ABS.h=mode.ABSh;
    ABS.Dope=Dope;
    ABS.Doph=Doph;
    
    %'qui ABS f_callVelm', keyboard
    %'controllo T', keyboard
    %  zitr=zT;
    %  fitr=1;
    %Tz=dndT.*ones(size(zT));
    
    
    xT=NPx;
    
    %  'Tg', keyboard
    ax=dox*3;
    Ps.ax=ax;
    %ax=dox*1.3;
    
    
    
    
    if T0>0
        Ps.T=Tg;
        Ps.Tx=xT;
        % Ps.Tz=zT+dz;
        %   Ps.Tz=Tz;
        Ps.zT=calop.zeta;
        %   ndz=3;
        %   Ps.ndz=ndz;
        Ps.zitr=calop.zeta;
        %   Ps.fitr=fitr;
        
        % figure,
        % subplot(121), plot(xT,Tdd)
        % xlabel('ro')
        % ylabel('Delta T')
        %  %'qui Temp', keyboard
        %  subplot(122), plot(zT,Tdd)
        %   xlabel('z')
        %   pausak
        
        if(ifp_Opt==-10)
            
            figure,
            subplot(121), plot(xT,Tdd)
            xlabel('ro')
            ylabel('Delta Eps (Real)')
            %'qui Temp', keyboard
            subplot(122), plot(zT,Tdd)
            xlabel('z')
            title('real part')
            pausak
            figure,
            subplot(121), plot(xT,imag(Tg)*(2*k0/100))
            xlabel('ro')
            ylabel('Delta Eps (Imag)')
            %'qui Temp', keyboard
            subplot(122), plot(zT,imag(Tg)*(2*k0/100))
            title(' LOSSES')
            xlabel('z')
            
            % 'qui Temp', keyboard
        end
        
    else
        Tg=0;
        %    'qui Temp', keyboard
        
        
    end
    
    %Tg=0
    
    calop.epsT=Tg.';
    
    calop.ro_in=ro_in';
    %'ro_in', keyboard
    calop.N=N;
    calop.Nref=Nref;
    
    %calop.Nref=mode.matDn;
    calop.ndz=25;
    
    
    
    calop.del_n_ag=del_n_ag;
    
    if del_n_ag==0
        ianti_gui=0;
    else
        Dlai=Dlam;
        %Dlai(1:2)=Dlai(1:2)-1.3*calop.del_n_ag;
        calop.Dlam     =Dlai;
    end
    calop.ianti_gui=ianti_gui;
    
    if isfield(VelmOptions,'isoga')
        Ps.isoga=VelmOptions.isoga;
    else
        Ps.isoga=1;
    end
    
    %'cont matDn', keyboard
    % da riempire per SPP
    Par=0;
    Ps.Par=Par;
    calop.Ps=Ps;
    
    
    %calop.N=0;
    %calop.Nref=0;
    %calop.epsT=0;
    %calop.ianti_gui=0;
    %calop.pimra=0;
    %calop.del_n_ag=0;
    %' ho messo tutto a ZERO in f_CallVELM',
    % calop.epsT=calop.epsT*10;
    
else
    %% 1D !!!!, metto parametri tipici per 1 modo 3D, tanto non lo calcolo
    MODC0=0;
    ro_campo=0;
    ro_in=0;
    iLP=1;
    mod_max=0;
    numodiacc=0;
    dk1=.01;
    alim=.12;
    Ev_or_Od='Even';
    pasnu=2;
    Dlam=[0 1 3 0 .1];
    Nref=0;
    del_n_ag=0;
    pimra=0;
    
    lambda=mean(mode.vlambda)*1e-9; %  from nm to m
    
    if length(mode.Deltan)>1
        alpha=mode.alpha*1e-6;  % from cm-3 to m-3
        k0=2*pi/lambda;
        Deltan=mode.Deltan;
        ImagDeltan=alpha./(2.*k0);
    else
        Deltan=0;
        ImagDeltan=0;
    end
    
    
    for k=1:length(ParVet)
        par_in{k}=ParVet(k);
    end
    NPx1=10;
    lambda=mean(mode.vlambda)*1e-3;
    Ps.i1D=1;
    Ps.itutmir=itutmir;    % defined in "VELMparameters.m", 0 for reduced thermal, 1 for full
    
    setVELMparameters_1D                                 %-- chiamare i parametri che non devo settare
    
    if isfield(mode,'DeltaN')==0   % From 4D LUT, "NUOVA_Gain.m"
        Dn=0;
    else
        Dn=mode.DeltaN;
    end
    
    del_n_agInf=Dn(end);
    del_n_ag0=Dn(1);
    del_n_ag=del_n_ag0-del_n_agInf;
    if del_n_ag==0
        Nref=0;
    else
        Nref=(Dn-del_n_agInf)/del_n_ag;
        
        mode.fat_ag=1;
        del_n_ag=del_n_ag*mode.fat_ag;
        
        %del_n_ag=del_n_ag*(1-Nref_du(end));
        %   'cont dfat[', keyboard
        
    end
    
    % Material gain
    matgain=mode.matgain; % Tibaldi
    matgain=matgain(end);    % VELM double column
    
    if isfield(VelmOptions,'gain_gui')
        igain=VelmOptions.gain_gui;
        if igain==0
            matgain=0;
        end
    end
    if length(matgain)==1
        N=0;
        %    vph=c_light*100/rrat;
        %    g0=abs(matgain)/vph;
        g0=abs(matgain);
        iga=0;
        if iga==1
            roadd=8;
            ro_in=sort([ro_in roadd+.01]);
            fiN=find(ro_in<=roadd);
            N=zeros(size(ro_in));
            N(fiN)=1;
        else
            g0=0;
        end
        %g0=0;
        %    'quio', keyboard
    else
        matgain=matgain(:,fianti);
        if max(matgain)>-100
            g0=abs(min(matgain));  % in 1/cm % Tibaldi
            N=matgain-min(matgain); N=N./max(N);
            if N(1)<0
                N=-N;
            end
            g0=abs((matgain(end)));  % in 1/cm % Tibaldi
            N=matgain-matgain(end); N=N./max(N);
            if N(1)<0
                N=-N;
            end
        else
            roadd=dox+1;
            N=exp(-(ro_in/roadd).^4);
            g0=0;
            g0=abs(min(matgain));
        end
        %    if(max(matgain)<0)
        %%     g0=0;
        %     velm.Lm=mode.Lm; % prendo solo il primo modo (per ora)
        %     velm.vlambda=mode.vlambda;
        %     velm.NQW=mode.NQW; % prendo solo il primo modo (per ora)
        %     velm.Gamma_z=mode.Gamma_z;
        %     velm.nindexQW=mode.nindexQW;
        %     velm.fPES=mode.fPES;
        %     velm.fPdif=mode.fPdif;
        %     velm.E2=mode.E2; % original
        %    end
    end
    pimra=-g0/cconv_velm;
    
    calop.pimra=pimra;
    calop.ianti_gui=ianti_gui;
    calop.Nref=Nref;
    calop.del_n_ag=del_n_ag;
    calop.ndz=25;
    
%     calop.zeta=mesh.ygrid*10;
%     calop.zeta=mesh.node*1e4;
    calop.zeta=mesh.ygrid*1e4;
    
    dndT0=mesh.dndT;
    dndT=mesh.dndT1D;
    
        nmed=3.3;
%     nmed=3.248; %media di n con |E|^2

    %Ps.dndT0=dndT0*(mode.T0-300);   % valore 1D
    
    T300=mode.T300;
    
%     puntatore=mesh.puntatore;
    puntatore=1:mesh.nny;
    
    deltaT=fliplr(mesh.DeltaTvelm(puntatore));
    
    DT_Peltier=mode.T0-T300;
    
    Tdd=deltaT+DT_Peltier;  % riduco la zona termica trasversale a dove il campo la sente
    
    T_por=0;
    if isempty(mode.ii_dd)==0
        T_por=flipud(Deltan-1i.*ImagDeltan); % carriers, look LUTs
        T_por=T_por(puntatore); % carriers, look LUTs
    end
%     T_por=0;
    
    T_ter=dndT0*DT_Peltier+dndT*deltaT;
    if length(T_ter)==1    % at first VELM iteration, Tdd is a single value, bring it on a constant vector!
        T_ter=T_ter*ones(size(calop.zeta));
    end
%     Tg=2*nmed*(T_ter+T_por);    % VENUS 3
    Tg=(T_ter+T_por);
    calop.epsT=Tg;
    
    zTf=calop.zeta(puntatore);   % excluding the contacts...
    zAbs=max(zTf)-fliplr(zTf);
    %  zAbs=zTf;
    ABS.zAbs=zAbs;
    
    %  global EO ABS
    
    dz=0;
    zT=(max(zTf)-fliplr(zTf))+dz;
    calop.zeta=zT;
	
if(ifp_Opt==-10)


 figure, 

 ylabel('Incremento Temperatura')
  %'qui Temp', keyboard
  plot(zT,Tdd)
   xlabel('z')
%   title('real part')
   pausak

  % 'qui Temp', keyboard
end	
		
	
    EO.zT=zT;
    
    if isfield(mode,'Fat_PerCoefTemp')
        Fat_Perd=mode.Fat_Perd0+Tdd*mode.Fat_PerCoefTemp;
        ABSh=mode.ABSh.*Fat_Perd;
        ABSe=mode.ABSe.*Fat_Perd;
    else
        ABSh=mode.ABSh;
        ABSe=mode.ABSe;
    end
    
    ABS.e=ABSe;
    ABS.h=ABSh;
    
    elecABS=mode.elecABSvelm(puntatore);%/1e18*mode.CarrierNorm;
    holeABS=mode.holeABSvelm(puntatore);%/1e18*mode.CarrierNorm;
    
    ABS.xmol=mesh.xmol(puntatore)';
    
    if isfield(mode,'PerCoefExT')
        ABS_Texp=mode.ABS_Texp+mode.PerCoefExT*Tdd'; % PerCoefExt=0
    else
        ABS_Texp=mode.ABS_Texp;
    end
    
    ABS.eleccentro=elecABS';
    ABS.holecentro=holeABS';
    
    ABS.eleccentro=ABS.eleccentro.*(1+Tdd'./T300).^ABS_Texp;
    ABS.holecentro=ABS.holecentro.*(1+Tdd'./T300).^ABS_Texp;
    
    ABS.Tlosses=mode.ABSTlosses;
    
end


%  ' prima di calop in set_opt ', keyboard
if ifp_Opt==-10
    ' prima di calop in set_opt ', keyboard
end
[St_wave,camzT,cam2_0,gain_0,lamod_0,modc,xro,fian,fPdif_0,fPES_0,fPGA,...
    tyPmod_0,omP0_0,eta_eff_0,Tper,Pa,Ppol,Plot,Eo_x,Eo_y,ord,velm1D]=...
    caloptdd_1D(MODC0,ro_campo,ro_in,fil_str,iLP,ifp_Opt,calop);

fiInv=find(zAbs<=camzT(end,1));
zInviluppo=zAbs(fiInv);
cazfit=zeros(size(zAbs));

cazfit(fiInv)=spline(camzT(:,1),(camzT(:,2)),zInviluppo);
cazfit=cazfit/cazfit(1);
% finv=find(ABS.zAbs<min(camzT(:,1)) | ABS.zAbs>max(camzT(:,1)));
% if length(finv)>0
%  cazfit(finv)=0;
% end
Inviluppo_SW=fliplr(cazfit);
velm1D.Inviluppo_SW=Inviluppo_SW;


velm=velm1D;


%'dopo SW  fiine', keyboard

