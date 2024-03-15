% set(gca,'FontSize',14,'FontName','Arial','box','on')
%==========================================================================
% Code for the caluclation of the band diagram of
% InGaAs/GaAsP QW emitting at 1.06 um here using model solid theory
% the strain used comes from composition file
%==========================================================================
%% costants
[GaAs,InAs,AlAs]=Costant_GaP_InAs_GaAs(mesh);
%% Initialization


%xmol_barr=1; %980nm
Thick_well=0.6; %
%E_lam=linspace(1,1.2,100);
% for i=1:length(E_lam)
%     Emis_lambd=E_lam(i);
%% retrieve molar composition and strain
[Thick_barr, strain,strain_b]= Composition_strain_mol_frac(xmol_well,Thick_well,xmol_barr); % strain well is used
% depending on the type of strain we have to consider hh or lh. tensile-> +
% compressive-> -
%  strain=0;
%  strain_b=0;
%% band calculation
%==========================================================================
% Well InGaAs 
%==========================================================================
% AlAs always from the same reference y=1
%well.Eg= 0.36 +0.629*xmol_well + 0.426*xmol_well.^2; 
well.Eg= xmol_well*GaAs.Eg + (1-xmol_well)*InAs.Eg; 


well.av=xmol_well*GaAs.av+(1-xmol_well)*InAs.av;
well.ac=xmol_well*GaAs.ac+(1-xmol_well)*InAs.ac-xmol_well*(1-xmol_well)*2.61; % here  bowing parameters non zero
well.b=xmol_well*GaAs.b+(1-xmol_well)*InAs.b;
well.C11=xmol_well*GaAs.C11+(1-xmol_well)*InAs.C11;
well.C12=xmol_well*GaAs.C12+(1-xmol_well)*InAs.C12;
well.Ev_av=xmol_well*GaAs.Ev_av+(1-xmol_well)*InAs.Ev_av;
well.Ev=xmol_well*GaAs.Ev+(1-xmol_well)*InAs.Ev;
well.Ec=xmol_well*GaAs.Ec+(1-xmol_well)*InAs.Ec;
well.so=xmol_well*GaAs.so+(1-xmol_well)*InAs.so-xmol_well*(1-xmol_well)*0.15;
well.mn=xmol_well*GaAs.mn+(1-xmol_well)*InAs.mn-xmol_well*(1-xmol_well)*0.0091;  
Ehv_w=2*well.av*(1-well.C12/well.C11)*strain;
Et_w=-2*well.b*(1+2*well.C12/well.C11)*strain;
Ehc_w=2*well.ac*(1-well.C12/well.C11)*strain;

%VB avar compressive case
Ehh_w=well.Ev_av+Ehv_w-0.5*Et_w+well.so/3; % error in zwang! - Et is correct !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 %Anche Chuang del 1999 mette il meno!
%CB
Ec_w=well.Ev_av+well.so/3+ well.Eg+Ehc_w;
% %Harrison model
% Ehh_w=well.Ev+Ehv_w-Et_w/2;
% 
% Ec_w=well.Ec+Ehc_w;
%==========================================================================
% Barrier GaAsP
%==========================================================================
% AlAs always from the same reference x=1
%xmol_barr=linspace(0,1,101);
%barr.Eg=2.776-1.469*xmol_barr+0.108*xmol_barr.^2;
barr.Eg=xmol_barr*GaAs.Eg+(1-xmol_barr)*AlAs.Eg;
%plot(xmol_barr,quadEg,xmol_barr,linEg,'.')
%keyboard
barr.av=xmol_barr*GaAs.av+(1-xmol_barr)*AlAs.av;
barr.ac=xmol_barr*GaAs.ac+(1-xmol_barr)*AlAs.ac;
barr.b=xmol_barr*GaAs.b+(1-xmol_barr)*AlAs.b;
barr.C11=xmol_barr*GaAs.C11+(1-xmol_barr)*AlAs.C11;
barr.C12=xmol_barr*GaAs.C12+(1-xmol_barr)*AlAs.C12;
barr.Ev_av=xmol_barr*GaAs.Ev_av+(1-xmol_barr)*AlAs.Ev_av;
barr.Ev=xmol_barr*GaAs.Ev+(1-xmol_barr)*AlAs.Ev;
barr.Ec=xmol_barr*GaAs.Ec+(1-xmol_barr)*AlAs.Ec;
barr.so=xmol_barr*GaAs.so+(1-xmol_barr)*AlAs.so;
barr.mn=xmol_barr*GaAs.mn+(1-xmol_barr)*AlAs.mn;
Ehv_b=2*barr.av*(1-barr.C12/barr.C11)*strain_b;
Et_b=-2*barr.b*(1+2*barr.C12/barr.C11)*strain_b;
Ehc_b=2*barr.ac*(1-barr.C12/barr.C11)*strain_b;

%VB avar  tensile case
Elh_b=barr.Ev_av+Ehv_b+Et_b/4-barr.so/6+ 0.5*sqrt(barr.so^2+barr.so*Et_b+(9/4)*Et_b^2);
%CB
Ec_b=barr.Ev_av+barr.so/3+ barr.Eg+Ehc_b;

% %Harrison model
% Elh_b=barr.Ev+Ehv_b+Et_b/2;
% 
% Ec_b=barr.Ec+Ehc_b;

gap=Ec_w-Ehh_w;
%  AlAs=well.Eg;

%% final quantity + plot of BD
Delta_Ev=-Ehh_w+Elh_b;
Delta_Ec=Ec_b-Ec_w;%band offset
Qc=Delta_Ec/(Delta_Ec+-Delta_Ev);

