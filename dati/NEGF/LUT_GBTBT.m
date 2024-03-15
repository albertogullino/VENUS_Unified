% load LUT_GBTBTonT % only needed for the Tvet values in the LUT
% coefficienti del fit della corrente ricavatao in Main_3
% load Interp_Fit_AlGaAs_Temp     % Main_3, T dependence

iLtunn=mode.iLtunn; % 1: fitting on the tunneling length (see "TJcharacteristic_test.m")

if strfind(mode.strName,'Julian')
    if iLtunn==1
%         load Interp_btj_GaAs_Julian_Nd=3e19_ltunn
        indTJ_R=mesh.inBTJ{1}(1);
        indTJ_L=mesh.inBTJ{2}(1);

        dA=strrep(num2str(mesh.dop_a(indTJ_L)),'+','');
        dD=strrep(num2str(mesh.dop_d(indTJ_R)),'+','');

        if strfind('Full',mode.ionization)==1
            strIoniz='';
        else
            strIoniz='_Incom';
        end
        prompt = 'Any suffix for the NEGF LUT? ';
        sffx_append=['_',input(prompt,'s')];
        if mesh.xmol(indTJ_R)==0
% %             load(['Interp_btj_GaAs_Julian_Nd=',dD,strIoniz,'_BGNdouble_ltunn.mat'])
%             load(['Interp_btj_GaAs_Julian_Nd=',dD,strIoniz,'_DeltaE=-200meV_ltunn.mat'])
            load(['Interp_btj_GaAs_Julian_Nd=',dD,strIoniz,sffx_append,'ltunn.mat'])
        else
            load(['Interp_btj_AlGaAs_Julian_Nd=',dD,'_x=',num2str(mesh.xmol(indTJ_R)),strIoniz,sffx_append,'_ltunn.mat'])
        end
        clog=ctunn;
        I0=0;

        mode.J0=J0;
        mode.L0=L0;
    else
%         load Interp_btj_GaAs_Julian_log
        fprintf('Warning! Julian str with an OLD NEGF curve: check it!\n'), pausak

        load Interp_Fit_AlGaAs_Temp_log
    end
elseif strfind(mode.strName,'Chalmers')
    % 2020Tibaldi_PRAPP
    prompt = 'Any suffix for the NEGF LUT? ';
    sffx_append=['_',input(prompt,'s')];
    load(['Interp_btj_testChalmers2020',sffx_append,'_ltunn.mat'])
    
    clog=ctunn;
    I0=0;

    mode.J0=J0;
    mode.L0=L0;
else
    % 2024Gullino_PJ
    load Interp_Fit_AlGaAs_Temp_log
end

if exist('T_vet')
    mode.T_NEGF=T_vet;
else
    mode.T_NEGF=293;
end
% mode.cBTBT=c;
mode.cBTBT=clog;
mode.I0_NEGF=I0;