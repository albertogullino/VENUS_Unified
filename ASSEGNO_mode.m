
%ParVet=[Oxide,Contact]; % metterli diversi tra loro, anche con differenza 1e-8;
if ~exist('iAndrea')
    iAndrea=0;
end



if iAndrea==1
    if irel==1
        ParVet=[Oxide,Relief]; % metterli diversi tra loro, anche con differenza 1e-8;
    else
        ParVet=[Oxide];         % metterli diversi tra loro, anche con differenza 1e-8;
    end
else
    
    Oxide=OX(Isize)/2;          % oxide radius, um (already defined in "settings")
    
%     if exist('DelOx')
%         OX1=OX+.4;              % (already defined in "settings")
%         OX2=OX1+DelOx;
%         OX3=OX2+DelOx;
%         %         if mode.flgBTJ==1
%         %             clear OX2 OX3
%         %         end
%     end
    
    %     COnt=MEsa-Width_Contact*2;
    %     Contact=COnt(Isize)/2;      % Contact radius, um
    %     Contact_e=Contact+Width_Contact;
%     VelmOptions.krel_max=KMax(Isize);
    ParVet=[Oxide,Contact,Contact_e]; % metterli diversi tra loro, anche con differenza 1e-8;
    
    if exist('OX2')
        ParVet=[Oxide,Contact,Contact_e,OX1(Isize)/2,OX2(Isize)/2];
    end
    if exist('OX3')
        ParVet=[Oxide,Contact,Contact_e,OX1(Isize)/2,OX2(Isize)/2,OX3(Isize)/2];
    end
    if exist('OX4')
        ParVet=[Oxide,Contact,Contact_e,OX1(Isize)/2,OX2(Isize)/2,OX3(Isize)/2,OX4(Isize)/2];
    end
    
    %     Oxide=3;
    
    
    if isfield(mode,'quasi1D')
        
        if mode.quasi1D==1  % quasi 1D VENUS (VCSEL and resistor 1D)
            ParVet=[Oxide Oxide Oxide Oxide Oxide Oxide];
            
            if mode.flgBTJ==1
                                ParVet=[Oxide Oxide Oxide];
%                 ParVet=[Oxide Oxide Oxide Oxide Oxide Oxide];
            end
        end
        
    end
    
    if irel==1
        ParVet=[ParVet RelSize/2];
        % keyboard
        
    end
    
    if isfield(mode,'flgBTJ_lithographic')
        if mode.flgBTJ_lithographic>0
            ParVet(4)=PulLayTJ;
            %  ParVet(5)=ThickEtchTJ;
            %  ParVet(6)=P6;
            % 	'qui parvet', keyboard
        end
    end
    
    
    mode.AreaOx=pi*(Oxide*1e-4)^2;
    mode.AreaOpt=pi*((Oxide+Oxide*mode.OptScaling/100)*1e-4)^2;
    
    
    %'NIM', keyboard
    %     VelmOptions.num_azim=NUM_Azim; % example, 3 is 0, 1, 2
    %     if NUMERO_MODI<=2
    %         VelmOptions.Dlam=[-1 2.5 5 0 .4]; % ORIGINAL
    %         VelmOptions.Dlam=[1 8 5 0 .4]; % ORIGINAL, 2 modes, dox=2 um
    %     end
    if isfield(mode,'quasi1D') && mode.quasi1D==1
        VelmOptions.dlam = .02;
        VelmOptions.NPlam = 10;
        VelmOptions.dndT0=2.3e-4;
    end
    %    'Parvet ',keyboard
end

if exist('iBard')
    if iBard==1
        ParVet=[ParVet nLC];
    end
end


% 'ass;,', keyboard
mode.IOLDsw=IOLDsw;
mode.FLos=FLos;
mode.Qc=Qc;
mode.Deltalam=Deltalam; % in nm
mode.ExpH=ExpH;
mode.ExpE=FattoreExpE*ExpH;
mode.STIMA_Temp=STIMA_Temp;

%'qui', keyboard
if IHILS==1
    mode.NxCoe=NxCoe;
    mode.CoeffHilsum0=N_X;
else
    if IOLDsw==0
        mode.CoeffHilsum=N_X+(T0-T300)*NxCoe*1e18; % Coefficient N_X for the Hilsum model
    else
        mode.CoeffHilsum=N_X+(T0-T300)*NxCoe; % Coefficient N_X for the Hilsum model
    end
    
    if isfield(mode,'NxCoe')
        mode=rmfield(mode,'NxCoe');
    end
    
    if isfield(mode,'CoeffHilsum0')
        mode=rmfield(mode,'CoeffHilsum0');
    end
end




%'N_X', keyboard
if exist('tHCoe')
    tauHt=tauH+(T0-T300)*tHCoe; %
else
    tauHt=tauH; %
end
mode.tausH=tauHt*1e-12;
mode.tausE=FattoreTauE*tauHt*1e-12;
mode.T_tauscat=T_tauscat;
mode.Tcap_EXP=Tcap_EXP;
mode.tauRat=tauRat;
mode.fat_gain=fat_gain;
if exist('epsNLg')
    mode.epsNLg=epsNLg;
end
if exist('fat_gainG')
    mode.fat_gainG=fat_gainG;
end


%mesh.fCondTer=fCondTer;
%mesh.fCondTerZ=fCondTer*Fat_CondZ;
mode.T0=T0;
mode.T300=T300;
mode.Exp_Temp0=Exp_Temp0;
mode.Tstart=Tstart;
mode.iTappo=iTappo;
VelmOptions.itutmir=itutmir; % 1: caso termico completo; 0: caso ridotto (+ veloce)
mode.fat_RAD=fat_RAD;
mode.AlTarocco=AlTarocco;



if exist('ABS_Texp')
    mode.ABS_Texp=ABS_Texp;
end

mode.Fat_Dop=Fat_Dop;
mode.Fat_N2=Fat_N2;
mode.Fat_P2=Fat_P2;
mode.fat_ag=fat_ag;

effetti=effetti0;
%                                              5                          7                 9
%          [fPES  Gam_z fPdif  E2 lambda TempVELM anti_gui  Diffus  Str  Lut]

if Effetti>0
    effetti(Effetti)=1
    if Effetti==8
        FAT_Diff_E=2;
        FAT_Diff_H=2;
    elseif Effetti==9
        iStruttura=6;
    elseif Effetti==10
        iLUT=6;
    end
end



mode.effetti=effetti;
%'effetti', keyboard


NOMELUT=LUT{iLUT};
mode.Fasano=[nomeSR,'Fasano'];
mode.GLUT=[NOMELUT,'_Der.mat'];
geom.GLUTm=[NOMELUT,'_more.mat'];
% geom.GLUTm='LUT4D_Feb24_Markus_nMark_xW=0_xB=0.2_more.mat';
% geom.GLUTm='LUT4D_Feb24_Markus_nMark_xW=0_xB=0.38_more.mat';
mode.FAT_idiffusioneQW_E=FAT_Diff_E;
mode.FAT_idiffusioneQW_H=FAT_Diff_H;

mode.CN_Auger=CN_Auger;
mode.FatNP_Auger=FatNP_Auger;
mode.CTemp_Auger=CTemp_Auger;
mode.CTemp_Ion=CTemp_Ion;
mode.Brad = Brad; % Brad 3D
mode.FatAuger23D=FatAuger23D;


mode.C_Temp=C_Temp;
mode.C_TempGain=C_TempGain;
mode.Fat_regeneration=Fat_regeneration;
mode.Auger_broad=Auger_broad;
mode.C_Temp_DD=C_Temp_DD;


if Fat_regeneration<0
    mode.AugerExpulsion=0; % se 1, attiva il modello di Auger expulsion
end

if Zmat>=0
    mode.Zmat=Zmat; % linear network embedding nonlinear device
    structureName=STR{iStruttura}
else
    mode.Zmat=0; % linear network embedding nonlinear device
    structureName=STR{abs(Zmat)}
end
rad_setting
%'qui', keyboard
if iloop==0
    eval(['save ',nomeSW,'Contributi_',nomeSav,' ',' effetti IPAR PMAT Svar kpar Effetti LUT STR -append'])
    %    save Contributi effetti IPAR PMAT Svar kpar Effetti LUT STR -append
end


if FLos>10
    mode.mintempVELM=500; % degrees, minimum temperature such that VELM is called
end

mode.ABS_Apor=ABS_Apor;

% Al=ABS_Apor;
% N0=ABS_Ader;
% N00=6;
% InDe=(Al*N0+1)/Al;
% alco=InDe*log(Al*N00+1);
% alco3=N00;
% if Al>0
% Fat_Perd=Fat_Perd0*alco3/alco
% else
%  Fat_Perd=Fat_Perd0;
% end



Fat_Perd_mod=Fat_Perd0+Fat_PerCoefTemp*(mode.T0-T300);
mode.Fat_Perd0=Fat_Perd0;

mode.PerCoefExT=PerCoefExT;


if exist('Fat_PerCoefTemp')
    %   mode.Fat_PerCoefTemp=Fat_PerCoefTemp;
end
%'cont mode', keyboard

if isfield(mode,'ABSe0')
    mode.ABSe=mode.ABSe0*Fat_Perd_mod;
    mode.ABSh=mode.ABSh0*Fat_Perd_mod;
end



%PE=Fat_Perd*log(Al*N00+1);

% x=linspace(0,10,101);
% x0=1;
% InDe=(Al*N0+1)/Al;
% yy0=InDe*log(Al*x+1);
% y(:,kpar)=yy0*mode.ABSh;

% figure, plot(x,y,x,x*mode.ABSh0*Fat_Perd0)
%'ver', keyboard



mode.ABS_Ader=ABS_Ader;
mode.TARde=TARde;

if exist('ifFC')
    mode.ifFC=ifFC;
end

%ABS_Apor=.5;
%ABS_Ader=1;