%'ididn', keyboard

if exist('iDyn')
    if iDyn==1
        mode.fvet=fvet;
        Cf=(1000*mode.ii_dd(end));
        CfPerc = (1-CurDynRef/Cf)*100;
        
        [DifCur,fi]= min(abs(CfPerc));
        if DifCur < minPerc
            CurDynRef(fi)=1000;
            CurDyn=[CurDyn Cf];
            
            %             save(['NoTherm_',num2str(CurDynRef(fi)),'_mA'])
            
            dinamica
            %             dinamicaOLD
            
            leCuD=length(CurDyn);
            modePlot.AM(leCuD,:)=AM;
            modePlot.Pss(leCuD,:)=sum(Pst_ss);
            modePlot.Yss(leCuD,:)=Y_ss;
            modePlot.Css(leCuD,:)=C_ss;
            modePlot.Gss(leCuD,:)=G_ss;
            modePlot.fvet=fvet;
            modePlot.CurDyn=CurDyn;
        end
        
    end
end


%'entro savePlot', keyboard
modePlot.ii_dd=mode.ii_dd;
%                modePlot.len=len;
if(mode.oflg)
    modePlot.nMaxVet=mode.nMaxVet;
    modePlot.pMaxVet=mode.pMaxVet;
    modePlot.n3MaxVet=mode.n3MaxVet;
    modePlot.p3MaxVet=mode.p3MaxVet;
    modePlot.nQW=mode.nQW;
    modePlot.NMQW=mesh.NMQW;
    
    modePlot.Pst_dd=mode.Pst_dd;
    modePlot.Psp_dd=mode.Psp_dd;
    modePlot.PspBulk_dd=mode.PspBulk_dd;
    
    modePlot.vind=vind;
    
end
%                modePlot.Vmeas=Vmeas;
%                modePlot.Imeas=Imeas;
%                modePlot.Imeas_res=Imeas_res;
%                modePlot.Rmeas=Rmeas;
modePlot.indv=indv;
modePlot.vv0_dd=mode.vv0_dd;
modePlot.DeltaTmax=mode.DeltaTmax;
if mode.Tflg==1
    modePlot.DeltaTmax_Joule=mode.DeltaTmax_Joule;
    modePlot.DeltaTmax_OptAbs=mode.DeltaTmax_OptAbs;
    modePlot.DeltaTmax_Ccap=mode.DeltaTmax_Ccap;
    modePlot.DeltaTmax_RAD=mode.DeltaTmax_RAD;
    modePlot.DeltaTmax_srhAu=mode.DeltaTmax_srhAu;
end

mode.DeltaTmax_Joule(indv)=DeltaTmax_Joule;
mode.DeltaTmax_srhAu(indv)=DeltaTmax_srhAu;
mode.DeltaTmax_Ccap(indv)=DeltaTmax_Ccap;
mode.DeltaTmax_RAD(indv)=DeltaTmax_RAD;
mode.DeltaTmax_OptAbs(indv)=DeltaTmax_OptAbs;

modePlot.DeltaTmax=mode.DeltaTmax;

modePlot.vv_dd=mode.vv_dd;
modePlot.NVbias=NVbias;

if mode.oflg==0
    fprintf(['\nCurrent: ',num2str(abs(mode.ii_dd(end))*1e3),' mA | Voltage: ',num2str(abs(mode.vv_dd(end))), ' V\n\n'])
else
    fprintf(['\nCurrent: ',num2str(abs(mode.ii_dd(end))*1e3),' mA | Voltage: ',num2str(abs(mode.vv_dd(end))), ' V | Opt. power: ',num2str(sum(mode.Pst_dd(:,end),1),'%.2f'), ' mW | dT(max): ',num2str(abs(mode.DeltaTmax(end)),'%.1f'), ' K\n\n'])
end

if mode.oflg==1
    modePlot.Gmod=mode.Gmod;
    modePlot.Lmod=mode.Lmod;
    modePlot.oflg=mode.oflg;
    modePlot.FLos=mode.FLos;
    modePlot.nmodes=mode.nmodes;
    modePlot.lambda=mode.lambda;
    
    modePlot.PTherm=mode.PTherm;
    if mode.v0_dd==0
        modePlot.PDiss=mode.PDissPred;
    end
    modePlot.matgain(indv,:)=mode.matgain;
    modePlot.fPdif(indv,:)=mode.fPdif;
    modePlot.E2(indv,:,:)=mode.E2;
    %                'in savePlot', keyboard
    modePlot.dn(indv,1:length(mode.DeltaN))=mode.DeltaN';
    modePlot.Tqw(indv,:)=mesh.T(mesh.inMQW{1})';
    modePlot.Elqw(indv,:)=mode.nQW{end}{2}';
    modePlot.Hoqw(indv,:)=mode.pQW{end}{2}';
    if isfield(mode,'Fat_cap_e')
        modePlot.Fate(indv,:)=mode.Fat_cap_e';
        modePlot.Fath(indv,:)=mode.Fat_cap_h';
    else
        modePlot.Fate(indv,:)=ones(size(mode.nQW{end}{2}'));
        modePlot.Fath(indv,:)=ones(size(mode.nQW{end}{2}'));
    end
    
    if isfield(mode,'Fat_cap')
        modePlot.FatCap(indv,:)=mode.Fat_cap';
    else
        modePlot.FatCap(indv,:)=ones(size(mode.nQW{end}{2}'));
    end
    modePlot.yQW=mesh.yMQW{2}*1e4;
end
modePlot.y=mesh.ygrid*1e4;
modePlot.x=mesh.xgrid*1e4;


Rcdu=reshape(mode.ecb,mesh.nny,mesh.nnx);
Rvdu=reshape(mode.evb,mesh.nny,mesh.nnx);
REfcdu=reshape(mode.EFn,mesh.nny,mesh.nnx);
REfvdu=reshape(mode.EFp,mesh.nny,mesh.nnx);
REldu=reshape(mode.elec,mesh.nny,mesh.nnx);
RHodu=reshape(mode.hole,mesh.nny,mesh.nnx);
REldo0=reshape(mesh.dop_d,mesh.nny,mesh.nnx);
RHodo0=reshape(mesh.dop_a,mesh.nny,mesh.nnx);
REldo=reshape(mode.dop_dp,mesh.nny,mesh.nnx);
RHodo=reshape(mode.dop_am,mesh.nny,mesh.nnx);

if mode.oflg==1
    REl2du=reshape(mode.N2D,mesh.nny,mesh.nnx);
    RHo2du=reshape(mode.P2D,mesh.nny,mesh.nnx);
    
    RRaug=sum(reshape(mode.RAugerQW,mesh.nny,mesh.nnx));
    RRsrh=sum(reshape(mode.RSRHQW,mesh.nny,mesh.nnx));
    RRlea=sum(reshape(mode.RLeakageQW,mesh.nny,mesh.nnx));
    Cn=sum(reshape(mode.Ccapn3D,mesh.nny,mesh.nnx))/NQW*1e-9*mesh.vWMQW{1};
    Cp=sum(reshape(mode.Ccapp3D,mesh.nny,mesh.nnx))/NQW*1e-9*mesh.vWMQW{1};
    
    Cnd=(reshape(mode.Ccapn3D,mesh.nny,mesh.nnx))/NQW*1e-9*mesh.vWMQW{1};
    Cpd=(reshape(mode.Ccapp3D,mesh.nny,mesh.nnx))/NQW*1e-9*mesh.vWMQW{1};
else
    Cnd=0;
    Cpd=0;
end

if mode.Tflg==1
    JN_X=reshape(mode.Jn_x,mesh.nny,mesh.nnx);
    JN_Y=reshape(mode.Jn_y,mesh.nny,mesh.nnx);
    JP_X=reshape(mode.Jp_x,mesh.nny,mesh.nnx);
    JP_Y=reshape(mode.Jp_y,mesh.nny,mesh.nnx);
    
    J_X=JN_X+JP_X;
    J_Y=JN_Y+JP_Y;
    J_Xn=JN_X;
    J_Yn=JN_Y;
    
    modePlot.JX(indv,:,:)=J_X;
    modePlot.JY(indv,:,:)=J_Y;
    modePlot.JXn(indv,:,:)=JN_X;
    modePlot.JYn(indv,:,:)=JN_Y;
    modePlot.JXp(indv,:,:)=JP_X;
    modePlot.JYp(indv,:,:)=JP_Y;
    
else
    [Jn_x,Jn_y,Jp_x,Jp_y] = f_EvalCurrentDensity(geom,mesh,mode);
    
    JN_X=reshape(Jn_x,mesh.nny,mesh.nnx);
    JN_Y=reshape(Jn_y,mesh.nny,mesh.nnx);
    JP_X=reshape(Jp_x,mesh.nny,mesh.nnx);
    JP_Y=reshape(Jp_y,mesh.nny,mesh.nnx);
    
    mode.Jn_x=Jn_x;
    mode.Jn_y=Jn_y;
    mode.Jp_x=Jp_x;
    mode.Jp_y=Jp_y;
    
    mode.JN_X=JN_X;
    mode.JN_Y=JN_Y;
    mode.JP_X=JP_X;
    mode.JP_Y=JP_Y;
    
end

%		'in savePlot', keyboard
%		RRaug=(reshape(mode.RAugerQW,mesh.nny,mesh.nnx));
%		RRaug=RRaug+(reshape(mode.RAuger,mesh.nny,mesh.nnx));
%		RRsrh=(reshape(mode.RSRHQW,mesh.nny,mesh.nnx));
%		RRsrh=RRsrh+(reshape(mode.RSRH,mesh.nny,mesh.nnx));
%		RRrad=(reshape(mode.RradQW,mesh.nny,mesh.nnx));
%		RRrad=RRrad+(reshape(mode.Rrad,mesh.nny,mesh.nnx));
%		Cnd=(reshape(mode.Ccapn3D,mesh.nny,mesh.nnx))/NQW*1e-9*mesh.vWMQW{1};
%		Cpd=(reshape(mode.Ccapp3D,mesh.nny,mesh.nnx))/NQW*1e-9*mesh.vWMQW{1};

RRaug=(reshape(mode.RAuger,mesh.nny,mesh.nnx));
RRsrh=(reshape(mode.RSRH,mesh.nny,mesh.nnx));
RRrad=(reshape(mode.Rrad,mesh.nny,mesh.nnx));

Rtot=RRaug+RRsrh+RRrad;
Rtn=Rtot-Cnd;
Rtp=Rtot-Cpd;

% figure, plot(x,RRaug,x,RRsrh,'--',x,RRlea,'.')
%xm=x*1e-2;
% WQW=7.9e-7;
% Iaug=qel*2*pi*trapz(x,x.*RRaug)*1000*WQW;
% Ilea=qel*2*pi*trapz(x,x.*RRlea)*1000*WQW;
% Isrh=qel*2*pi*trapz(x,x.*RRsrh)*1000*WQW;
% Icn=qel*2*pi*trapz(x,x.*Cn)*1000*WQW;
% Icp=qel*2*pi*trapz(x,x.*Cp)*1000*WQW;


ifigura=0;
if ifigura==1
    figure, plot(y,REfcdu(:,[1 end]),'linewidth',2)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(y,Rcdu(:,[1 end]))
    title('Elettroni')
    
    figure, plot(y,REfvdu(:,[1 end]),'linewidth',2)
    hold on
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(y,Rvdu(:,[1 end]))
    title('Lacune')
    
    figure,plot(mesh.ygrid*10000,mesh.xmol(1:mesh.nny))
    figure,plot(mesh.ygrid*10000,mode.elec(1:mesh.nny))
    figure,plot(mesh.ygrid*10000,mode.hole(1:mesh.nny))
end

modePlot.Rtn(indv,:,:)=Rtn;
modePlot.Rtp(indv,:,:)=Rtp;

if mode.oflg==1
    modePlot.Cn(indv,:)=Cn;
    modePlot.Cp(indv,:)=Cp;
    modePlot.El2D(indv,:)=REl2du(:,1)';
    modePlot.Ho2D(indv,:)=RHo2du(:,1)';
end
modePlot.Ec(indv,:)=Rcdu(:,1)';
modePlot.Ev(indv,:)=Rvdu(:,1)';
modePlot.EcH(indv,:)=Rcdu(:,fix(end/2))';
modePlot.EvH(indv,:)=Rvdu(:,fix(end/2))';
modePlot.EFc(indv,:)=REfcdu(:,1)';
modePlot.EFv(indv,:)=REfvdu(:,1)';
modePlot.El(indv,:)=REldu(:,1)';
modePlot.Ho(indv,:)=RHodu(:,1)';

%                'qui save', keyboard

modePlot.Eldr(indv,:)=REldo(:,1)';
modePlot.Hodr(indv,:)=RHodo(:,1)';
modePlot.Eldr0(indv,:)=REldo0(:,1)';
modePlot.Hodr0(indv,:)=RHodo0(:,1)';

%		REldo0=reshape(mesh.dop_d.elec,mesh.nny,mesh.nnx);
%		RHodo0=reshape(mesh.dop_a,mesh.nny,mesh.nnx);
%		REldo=reshape(mode.dop_dp.elec,mesh.nny,mesh.nnx);
%		RHodo=reshape(mode.dop_am,mesh.nny,mesh.nnx);


%                'stop batbn', keyboard

modePlot.xmol=mesh.xmol(1:mesh.nny);

if indv>1 && mode.oflg==1
    modePlot.Ic(indv,:)=curr_n;
    modePlot.Iv(indv,:)=curr_p;
    fP=abs((curr_n)+(curr_p));
    dz=[0 diff(z)];
    Cint=fP.*dz;
    curM=sum(Cint(3:end-3))/sum(dz(3:end-3));
    modePlot.CurM(indv)=curM;
    modePlot.Ind=Ind;
    modePlot.Zl=Zl;
    modePlot.zLim=zLim;
    modePlot.IntCcapN=mode.IntCcapN;
    modePlot.IntCcapP=mode.IntCcapP;
    modePlot.IntRec=mode.IntREc;
    modePlot.IntRecN=mode.IntRecN;
    modePlot.IntRecP=mode.IntRecP;
    modePlot.Fleak=mode.Fleak;
end
% da togliere

if mode.oflg==1
    modePlot.nQW=mode.nQW;
    modePlot.pQW=mode.pQW;
end
modePlot.irel=irel;

modePlot.fat_gain=mode.fat_gain;
modePlot.Last_Workspace=Last_Workspace;
modePlot.fCondTer=fCondTer;
modePlot.fCondTerZ=fCondTerZ;
modePlot.fatt_dndT=fatt_dndT; %

if isfield(mode,'CoeffHilsum')
    modePlot.CoeffHilsum= mode.CoeffHilsum; % parte bene
end

if isfield(mode,'NxCoe')
    modePlot.NxCoe= mode.NxCoe; % parte bene
    modePlot.CoeffHilsum0= mode.CoeffHilsum0; % parte bene
end

modePlot.Tcap_EXP= mode.Tcap_EXP;
modePlot.ExpH= mode.ExpH;
modePlot.ExpE= mode.ExpE;
modePlot.fat_RAD= mode.fat_RAD;
modePlot.fat_gain= mode.fat_gain;
modePlot.AlTarocco= mode.AlTarocco;
modePlot.T_tauscat= mode.T_tauscat;
modePlot.tauRat= mode.tauRat;
modePlot.Deltalam= mode.Deltalam; % in nm ; % forse 3.5 è l'ideale per lo Standard
modePlot.TAROCCO= TAROCCO;
%
modePlot.idiffusioneQW= mode.idiffusioneQW; % 0 nulla per port qw; 1: solo diffusione;  2 DD; 3 DD con gamma
modePlot.iambipolar= mode.iambipolarQW;
modePlot.FAT_idiffusioneQW_E= mode.FAT_idiffusioneQW_E;
modePlot.FAT_idiffusioneQW_H= mode.FAT_idiffusioneQW_H;
modePlot.tausE= mode.tausE;
modePlot.tausH= mode.tausH;
if mode.oflg==1
    modePlot.Gamma= sum(mode.Gamma_z);
end
modePlot.VelmOptions= VelmOptions;
% VelmOptions.ivett=1;
modePlot.Exp_Temp0= mode.Exp_Temp0;

modePlot.fat_ag= mode.fat_ag;     % fattore antiguiding
modePlot.ParVet= ParVet;     % fattore antiguiding


%    modePlot.Contact=Contact;
modePlot.Relief=Relief;
modePlot.Oxide=Oxide;
modePlot.ExpH=ExpH;
modePlot.N_X=N_X;
modePlot.ABSe=mode.ABSe;
modePlot.ABSh=mode.ABSh;
modePlot.ABS_Texp=mode.ABS_Texp;
modePlot.T_tauscat=T_tauscat;
modePlot.Tcap_EXP=Tcap_EXP;
modePlot.tauRat=tauRat;
modePlot.T0=T0;
modePlot.iLUT=iLUT;
modePlot.iStruttura=iStruttura;
modePlot.iTappo=iTappo;
modePlot.FattoreTauE=FattoreTauE;
modePlot.AlTarocco=AlTarocco;
modePlot.CN_Auger=CN_Auger;
modePlot.FatNP_Auger=FatNP_Auger;
modePlot.CTemp_Auger=CTemp_Auger;
modePlot.CTemp_Ion=CTemp_Ion;
modePlot.Fat_regeneration=Fat_regeneration;
modePlot.Auger_broad=Auger_broad;
modePlot.C_Temp=C_Temp;

if isfield(mode,'TempEst')
    modePlot.TempEst=mode.TempEst;
end


modePlot.GLUT=mode.GLUT;
modePlot.LUT=LUT;
modePlot.NOMELUT=NOMELUT;
modePlot.structureName=structureName;
modePlot.settings=rad_setting;

modePlot.node=mesh.node;
modePlot.triangle=mesh.triangle;
modePlot.nn=mesh.nn;
modePlot.nt=mesh.nt;
modePlot.zox=StrDD.zox;
iRa=imag(StrDD.raggi);
rox=iRa(iRa>0);
%'qui rox', keyboard
modePlot.rox=StrTT.Rox;
modePlot.Contact_i=StrTT.ro_met;
modePlot.Contact_e=StrTT.ro_mesa;
modePlot.StrTT=StrTT;

if mesh.NMQW>0
    modePlot.yMQW=mesh.yMQW;
    modePlot.vWMQW=mesh.vWMQW;
    efield_rho = -diff(mode.phi(mesh.inMQW{1}))./diff(mesh.node(1,mesh.inMQW{1}));
    modePlot.efield_rho(indv,:)=efield_rho;
end
modePlot.geom=geom;
modePlot.Fasano=mode.Fasano;
if exist('Isize')
    modePlot.Isize=Isize;
end
if exist('Isize')
    modePlot.Isize=Isize;
end


if isfield(mode,'JnQW')
    for kp=1:NQW
        modePlot.JnQW(indv,:,kp)=mode.JnQW{kp};
        modePlot.JpQW(indv,:,kp)=mode.JpQW{kp};
    end
end

%figure,plot(mesh.xgrid(1:mesh.nnxQW{1}-1),efield)

%Efie=reshape(mode.efieldy,mesh.nny,mesh.nnx);
%'fine saveplot', keyboard
Epot=reshape(mode.phi,mesh.nny,mesh.nnx);
E0=Epot(:,1);
Ez=-diff(E0)./diff(mesh.ygrid');
modePlot.efield_z(indv,:)=Ez;
modePlot.Phi(indv,:,:)=Epot;
modePlot.elec(indv,:,:)=REldu;
modePlot.hole(indv,:,:)=RHodu;
Tedu=reshape(mesh.T,mesh.nny,mesh.nnx);
modePlot.Temp(indv,:,:)=Tedu;
if mode.Tflg==1
    if max(size(mesh.Tprec))>1
        modePlot.Tprec(indv,:,:)=mesh.Tprec; % input for the thermal solver
    end
    modePlot.HeatJoule(indv,:,:)=reshape(mode.HeatJoule,mesh.nny,mesh.nnx);
    modePlot.HeatRec_RAD(indv,:,:)=reshape(mode.HeatRec_RAD,mesh.nny,mesh.nnx);
    modePlot.HeatRec_Cap(indv,:,:)=reshape(mode.HeatRec_Cap,mesh.nny,mesh.nnx);
    modePlot.HeatOptAbs(indv,:,:)=reshape(mode.HeatOptAbs,mesh.nny,mesh.nnx);
    modePlot.HeatRec_13(indv,:,:)=reshape(mode.HeatRec_13,mesh.nny,mesh.nnx);
    
    HeatRec_RAD=reshape(mode.HeatRec_RAD,mesh.nny,mesh.nnx);
    HeatRec_Cap=reshape(mode.HeatRec_Cap,mesh.nny,mesh.nnx);
    HeatOptAbs=reshape(mode.HeatOptAbs,mesh.nny,mesh.nnx);
    HeatRec_nr=reshape(mode.HeatRec_13,mesh.nny,mesh.nnx);
end
Mob=reshape(mesh.mobn0_n,mesh.nny,mesh.nnx);
modePlot.Mob(indv,:,:)=Mob;
Nc=reshape(mesh.Nc,mesh.nny,mesh.nnx);
modePlot.Nc(indv,:,:)=Nc;
EFn=reshape(mode.EFn,mesh.nny,mesh.nnx);
modePlot.EFn(indv,:,:)=EFn;
EFp=reshape(mode.EFp,mesh.nny,mesh.nnx);
modePlot.EFp(indv,:,:)=EFp;
ecb=reshape(mode.ecb,mesh.nny,mesh.nnx);
evb=reshape(mode.evb,mesh.nny,mesh.nnx);
modePlot.ecb(indv,:,:)=ecb;
modePlot.evb(indv,:,:)=evb;