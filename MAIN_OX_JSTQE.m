%========================================================================76
% pin VCSEL (JSTQE article)
%========================================================================76
clc
clear
clear global
% close all


% try

colordef white
dbstop if error

addpathVENUS    % add the folders needed for a VENUS simulation

% Flag to avoid slow geom generation of lithographic structures (load geom file)
flgGEOM=1;     % 1, carica vecchia geom
flgGEOM=0;     % carica vecchia geom

% Dynamic Analysis
iDyn = 0;                   % flag, 1 to turn on the dynamic analysis
fvet = 5*logspace(8,10,81);   % frequencies for small-signal simulation, Hz
% fvet = logspace(7,11,121);   % frequencies for small-signal simulation, Hz
% CurDynRef = [1:3:13];       %     values of current where small-signal analysis is performed
CurDynRef = 1:1:3;       % values of current where small-signal analysis is performed

minPerc=10;

IOLDsw=0;

iSavNome=1;

% IPLOT=-1;  % Structure details + live plots of simulation results
IPLOT=1;  % live plots of simulation results
% IPLOT=-2;   % Structure details plot
% IPLOT=0;    % No intermediate plots

% Imassimo=1;  % massima corrente analizzata
PotMin=.1;        % potenza finale

%% 
prompt = 'Insert the prefix to append at the begin pmat file: ';
nomeSav = input(prompt,'s'); % Suffix appended at the end of the save file, to distinguish the various tries
%
nomeSave=[nomeSW,nomeSav];
%eval(['save ',nomeSave,' h'])

iStruttura=8;  % 6,7,8; look into "List_strNames.m" for a complete list of str files 
 
%  Last_Workspac='dud';  
Last_Workspac='LW';  

Last_Workspace=[nomeSW,Last_Workspac];

% radi='_JSTQE_dotLU3'; % JSTQE results
% radi='_OX'; % OC-VCSEL: use for comparison with TJ!
radi='_OXelementi'; % OC-VCSEL: use for comparison with TJ!
% radi='_OXelementi_1D'; % OC-VCSEL: use for comparison with TJ!

rad_settingV{100}=radi;    %vale anche per IPAR = 0
for k=1:70
    rad_settingV{k}=radi;
end


%% Investigated parameters

% VELMinput='VelmSe';
VELMinput='VelmSa';
% VELMinput='Velmdu1';


iPrimo=1;

% IPvet=[0];  33
% IPvet=[2 3 23 55 25 27 10 22 5 44 42 33];   % Auger
% IPvet=[10 22 5 44 42 33];  % Set 2
% IPvet=[23 55 2 3 25 27];  % Set 1
IPvet=11;

IREST=0;  % HA EFFETTO solo con IPvet=0 e 30
%  IREST=2;  % HA EFFETTO solo con IPvet=0 e 30  USARE IREST=2
Isize=4;
% Isize=2;

% Room temperature
Temperature = 20;       % Celsius
T0 = Temperature+273;   % Kelvin
FLos=1;

HARD=1;       % = 0 toglie VELM        
% effetti=[0 0 0 1 0];  %Effetti_E
% effetti=[0 0 0 0 0];  %Effetti_0
% effetti0=[0 0 1 0 0];  %Effetti_fPdif
% effetti0=[0 0 0 0 0 0];  %Effetti_fPdif
effetti0=[0 0 0 0 0 0 0 0 0 0];  %Effetti_Temp
%                     5
%                                              5                          7                 9
%          [fPES  Gam_z fPdif  E2 lambda TempVELM anti_gui  Diffus  Str  Lut]
            
effetti=effetti0;

%'dopo primo save', keyboard

%iloop=1          % =1 provo i loops, to launch a simulation set it to 0
if IPvet(1)~=0
    iba=0;
    if iba==0
        iloop=input(' Controllo loops ? [1, SI; Enter, NO] ');
    else
        iloop=1;
    end
    if length(iloop)==0
        iloop=0;
    end
else
    iloop=0;
end

if iloop==0 && iSavNome==1
    Fun_Nome
    save 	 nomeSave nomeSW nomeSR nomeLUT
%    save ultimonome nomeSave nomeSW nomeSR
    eval(['save ',nomeSW,'Contributi_',nomeSav,' ',' HARD nomeSave nomeSR nomeSW effetti0 VELMinput'])
end

%% Definition of the bias condition (cell of tension: DDstr{itemp})
clear DDstr

% Bias values BELOW threshold condition
DDin=[0:0.4:1.2 1.3:0.1:1.5];

% Different bias conditions ABOVE threshold are investigated for different temperatures
TTve=[110:-30:20];
TTve=20

Tpelt=TTve;
for itemp=1:length(TTve)
    Temperaturei=TTve(itemp);
    if Temperaturei == 110
        Vadd=[1.55:0.02:2.4];
        Imassimo=8;  % massima corrente analizzata
    elseif Temperaturei == 80
        Vadd=[1.55:0.03:2.6];
        Imassimo=11;  % massima corrente analizzata
    elseif Temperaturei == 50
        Vadd=[1.55:0.04:2.8];
        Imassimo=13;  % massima corrente analizzata
    elseif Temperaturei == 20
%         Imassimo=16;  % massima corrente analizzata
        Vadd=[1.55:0.05:3.1];
%         Vadd=[1.62:0.02:1.8 1.85:0.05:2.3];
%         Imassimo=7;
        Imassimo=14;
    end
    DD0=[DDin Vadd];
%     DDstr{itemp}=DD0;
end

%keyboard

iold=2;
ivelMstop=0;
vv_double=0;

flagCallVELM=0;

iVELM=0;  % se 1, calcolo sempre VELM
ILOOP=1;         % flag per settings_vari: DEVE ESSERE SEMPRE A 1 nel programma parametrico

Gs.QWorientation='X'; % orient. quantum well (parallel to x (rho): VCSEL)
% geom.QWorientation='Y'; % orient. quantum well (parallel to y (z): microROD)
mode.iderGLUT=1;      % = 0 fa vecchio fit
mode.iderGLUTnew=2;   % = 0 fa vecchio fit

%% Variation of default parameters
SETTA_loopsTemp       % returns PMAT (cell containing the information for each parameter)
                      % IPvet extract the needed parameter

%keyboard
Ikcont=0;
Iclear=0;

for IPAR=IPvet
    
    if IPAR==0
        kpvet=1;
        rad_setting = rad_settingV{100};
    else
        kpvet=1:NP(IPAR);
        rad_setting = rad_settingV{IPAR};
        % VALORI di Riferimento in settings vari
        
    end
    eval(['settings_vari',rad_setting])
   
    % Save parameters in PMAT and save IPvet (i.e., to be analysed parameters)
    if iPrimo==1
        eval(['save ',nomeSave,'_pmat  PMAT IPvet tit nomeSave strSave irel nomeSW nomeSR nomeLUT'])
    end

    clear MODE MESH VInf VInp VelmOptions MODEplot
    
    dd=DD0;
    %'ver ee', keyboard
    
    for kpar = kpvet
        
        eval(['save facendo-',nomeSav,' IPAR kpar kpvet IPvet'])
        
        % In case they are different, it means that there some
        % inconsistence between SETTA_loopsTemp and the input given by the
        % user
        if exist('DDstr')
            if length(DDstr)~=length(kpvet)
                'Aggiustare lunghezza DDstr', keyboard 
            end
            dd=DDstr{kpar};
        end
        Ikcont = Ikcont+1;
        indv = 1;
        kfig = kpar-1;
        
%         mode1=0;
%         clear mode VELMInfo VELMInput VO
        mode1.a=0;
        clear mode VELMInfo VELMInput VO
        
        clear global Ps
        
        if iloop==0
            save sapar IPAR kpar
        end
        
        % Temperature dependence of lifetime
        if IPAR==6
            Tcap_EXP=1;  % Fat_cap=(1+(T-T0)/T_tauscat0).^Tcap_EXP; Inf: exponential
        end
        
        % assegno variabile corrente
        if IPAR>0
            % kpar take parameter value from PMAT{IPAR}
            ppar=PMAT{IPAR}(kpar);
            
            (['IPAR kpar =  ',num2str([IPAR, kpar])]),
            stringa_parametro=[Svar{IPAR},'=',num2str(ppar)]
            
            eval([stringa_parametro,';'])
        end
        eval(['settings_vari',rad_setting])
        %'controllo parametri', keyboard
        
        ASSEGNO_mode
        
        % In case of testing loop is activated
        if iloop==1
            pausak
        end
        
%        'controllo Par', keyboard
        
        % In case of testing loop is de-activated
        if iloop==0
            
            if (Ikcont==1 || Iclear==1)
                if IREST==0
                    'RIGENERO STRUTTURA'
                    structureName;
					if flgGEOM==0
						Sub_Str_from_VELM
					else
						fprintf('flgGEOM=1! Loading structure from "geom" file\n'),keyboard
						fis= strfind(structureName,'\');
						strName=structureName(fis(end)+1:end);
						DirName=structureName(1:fis(end));
						%%%% save geom
						load([DirName,'geom_' strName])
						ParMore=StrDD.ParMore;
						mode.nBTJ=StrDD.nBTJ;
						fprintf('Structure loaded\n'),keyboard
                    end
					
                    if mode.oflg==1
                        Glut4Dinc(mode)
                    end
                end
                
            else
                if IPAR==12
                    'Aggiorna LUT'
                    
                    Glut4Dinc(mode)
                elseif (IPAR==1 || IPAR==13 || IPAR==21 || IPAR==36 || IPAR==37 || IPAR==41 || IPAR==45 || IPAR==46)
                    'RIGENERO STRUTTURA'
                    structureName;
                    Sub_Str_from_VELM
                end
            end
            
        end %iloop
        if mode.oflg==1
            geom.QWorientation=Gs.QWorientation;
        end
        
        % 'ver par', keyboard
        if iloop==0
            
            ticSUB=tic;
            SUB_drive % for both current and voltage driving
            tocSUB=toc(ticSUB)/60        

            MODE{kpar}=mode;
            MESH{kpar}=mesh;
            
            if mode.oflg
                VInf{kpar}=VELMInfo;
                VO{kpar}=VelmOptions;
            end
            
            if exist('VELMInput')
                VInp{kpar}=VELMInput;
            end
            
            %  eval(['save ',nomeSave,num2str(IPAR),'.mat MODEplot VInf VInp VO'])
            %  eval(['save ',nomeSave,num2str(IPAR),'.mat MODEplot VO'])
            
            fprintf('fine ciclo\n')
            %keyboard
%             clear modePlot VELMInput VelmOptions VELMInfo mode mesh
            if (IPAR==1 || IPAR==13 || IPAR==21 || IPAR==29)
                Iclear=1;
                clear geom
            end
            % keyboard
        else
            %   mode.ExpE
            %   mode.Deltalam
            %   mode.CoeffHilsum
            %  mode.T_tauscat
            %  mode.tauE
            
            if IPAR==12
                NOMELUT
            elseif IPAR==13
                structureName
            end
            pausak
            
        end
        %pausak
        % close all
    end
    if IPAR<IPvet(end)
        kpvet=1:NP(IPAR+1);
        %'qui', keyboard
    end
    
end
%  catch
