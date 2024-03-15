function macro=mw_phmat(material,xmol,T,mesh)

%        T=295*ones(size(T));
%        Tg=295*ones(size(T));
    if isfield(mesh,'CarrierNorm')
        CarrierNorm=mesh.CarrierNorm;
    else
        CarrierNorm=1;
    end
        Tg=T;
   if isfield(mesh,'GapTemp')
    if mesh.GapTemp==0
     Tg=mesh.T0*ones(size(T));
    end
   end 

%'mw_phmat', keyboard
%%%% last update 12/06/2015
%Purpose: get physical properties of materials at microwave frequencies.
%
%Input variables:
%material            String which defines the material. Avaiable materials are:
%                    XY-LiNbO3','YX-LiNbO3', ...
%molar               Vector of molar fraction [initial final]
%X1                  x points of the layer
%Ny                  mesh points in the column
%
%Output variables:
%espr                Relative permittivity tensor (3-by-3 matrix). EPSR must be
%                    of the form:
%                    [epsrxx epsrxy 0     ]
%                    [epsryx epsryy 0     ]
%                    [0      0      epsrzz]
%mur                 Relative permeability tensor (3-by-3 matrix).
%                    MUR must be diagonal
%tgdelta             Loss tangent tensor of dielectric materials (3-by-3 matrix)
%sigma               Electric conductivity tensor of conductor materials (3-by-3 matrix), S/cm.
% loading constants
s_LoadConstants
Nxmol=length(xmol);
onesxmol=ones(1,Nxmol);
macro.xmol=xmol;
%
switch material
    
    case 'AlGaAs' % molar fraction dependent aluminum gallium arsenide
        %%
        %-- Quantum-corrected model parameters
        %
        macro.tauscatn = 1e-12*onesxmol; % quantum well capture time, s
        macro.tauscatp = 1e-12*onesxmol; % quantum well escape time, s
        %
        macro.tauscatn = mesh.tausE*onesxmol; % quantum well capture time, s
        macro.tauscatp = mesh.tausH*onesxmol; % quantum well escape time, s
        %
        %-- Incomplete ionization parameters
        T300=mesh.T300;
        [Ed,Ea]=ActivEner(xmol);
%         Ed=-150;
%         Ea=-150;
        macro.DeltaEa=1e-3*Ea.*(1+mesh.CTemp_Ion*(T-T300)./T300); % (eV)
        macro.DeltaEd=1e-3*Ed.*(1+mesh.CTemp_Ion*(T-T300)./T300); % (eV)

%        macro.DeltaEa=1e-3*25.*(1+mesh.CTemp_Ion*(T-300)./300); % (eV)
%        macro.DeltaEd=1e-3*6.*(1+mesh.CTemp_Ion*(T-300)./300); % (eV)        
        
        %-- Recombination models parameters
        macro.Etrap = 0*onesxmol; % (eV)
        % % SRH (Ioffe, xmol=0.1)
        macro.taun = mesh.taun*onesxmol.*(T300./T).^mesh.tauExp; % (s)
        macro.taup = mesh.taup*onesxmol.*(T300./T).^mesh.tauExp; % (s)
        macro.taunQW = mesh.taunQW*onesxmol.*(T300./T).^mesh.tauExp; % (s)
        macro.taupQW = mesh.taupQW*onesxmol.*(T300./T).^mesh.tauExp; % (s)
        % SRH (Michalzik book, chapter 3, Debernardi)
        % macro.taun = 5e-9*onesxmol.*(300./T); % (s)
        % macro.taup = 5e-9*onesxmol.*(300./T); % (s)
        %
        %- Radiative (Ioffe)
        macro.brad = mesh.Brad*onesxmol*CarrierNorm; % (cm^3/s)
        %
        %- Auger (Ioffe, xmol=0.1)
%         macro.Cnnp = 1.2e-31*onesxmol; % (cm^6/s)
%         macro.Cppn = 8.5e-31*onesxmol; % (cm^6/s)
%	'keyb Auger', keyboard
 %       [Cn,Cp]=fun_auger(xmol,T-mesh.T0);
 %       macro.Cnnp =Cn*10;
%	macro.Cppn =Cp*10;
	
        macro.Cnnp =mesh.CN_Auger*exp(mesh.CTemp_Auger*(T-T300)./100)*1e-30*CarrierNorm^2;
        %macro.Cnnp =mesh.CN_Auger*abs(1+mesh.CTemp_Auger*(T-mesh.T0)./300)*1e-30;
        macro.Cppn =mesh.FatNP_Auger.*macro.Cnnp;
        % Auger (Michalzik book, chapter 3, Debernardi)
        % macro.Cnnp = 3.5e-30*onesxmol; % (cm^6/s)
        % macro.Cppn = 3.5e-30*onesxmol; % (cm^6/s)
        %
        %%%%%%%%%% Nodes-dependent quantities %%%%%%%%%%
        %-- Bandgap (eV) (1985Adachi_JAP, eqs. (20)-(22))
        %-- Temperature dependent bandgap (eV) (1976Aspnes_PRB,
        %   http://www.ioffe.ru/SVA/NSM/Semicond/AlGaAs/bandstr.html)
        Eg0_G=1.519+1.155.*xmol+0.370.*xmol.^2;
        alpha_G=5.41e-4; beta_G=204;
        Eg_G=Eg0_G-(alpha_G.*Tg.^2)./(beta_G+Tg);
        %
        Eg0_X=1.981+0.124.*xmol+0.144.*xmol.^2;
        alpha_X=4.6e-4; beta_X=204;
        Eg_X=Eg0_X-(alpha_X.*Tg.^2)./(beta_X+Tg);
        %
        Eg0_L=1.815+0.69.*xmol; % mind the error in Ioffe!
        alpha_L=6.05e-4; beta_L=204;
        Eg_L=Eg0_L-(alpha_L.*Tg.^2)./(beta_L+Tg);
        
        %
%' qui Eg', keyboard        
        Eg = min(Eg_G,Eg_X);
        %
        % Calculating Eg in GaAs
        Eg0_G_GaAs=1.519;
        alpha_G=5.41e-4; beta_G=204;
        Eg_G_GaAs=Eg0_G_GaAs-(alpha_G.*Tg.^2)./(beta_G+Tg);
        %
        Eg0_X_GaAs=1.981;
        alpha_X=4.6e-4; beta_X=204;
        Eg_X_GaAs=Eg0_X_GaAs-(alpha_X.*Tg.^2)./(beta_X+Tg);
        %
        Eg0_L_GaAs=1.815; % mind the error in Ioffe!
        alpha_L=6.05e-4; beta_L=204;
        Eg_L_GaAs=Eg0_L_GaAs-(alpha_L.*Tg.^2)./(beta_L+Tg);
        %
        Eg_GaAs = min(Eg_G_GaAs,Eg_X_GaAs);

        %
        % saving bandgaps
        macro.Eg_G = Eg_G;
        macro.Eg_X = Eg_X;
        macro.Eg_L = Eg_L;
        macro.Eg = Eg;
        %
        %-- Effective DoS masses (1985Adachi_JAP, Table II)
        meffn_G = 0.067+0.083.*xmol;
        meffn_X = 0.850-0.140.*xmol;
        meffn_L = 0.560+0.100.*xmol;
        %
        % DoS
        Nc_G = 2*(2*pi*kB*T*m0).^(3/2)/(h.^3).*1e-6.*meffn_G.^(3/2)/CarrierNorm;
        Nc_X = 2*(2*pi*kB*T*m0).^(3/2)/(h.^3).*1e-6.*meffn_X.^(3/2)/CarrierNorm;
        Nc_L = 2*(2*pi*kB*T*m0).^(3/2)/(h.^3).*1e-6.*meffn_L.^(3/2)/CarrierNorm;
        %
        % Equivalent DoS
        Vt=kB.*T./qel;
        DeltaE_XG = Eg_X-Eg_G;
        DeltaE_LG = Eg_L-Eg_G;
        Nc_1 = Nc_G.*(1+Nc_X./Nc_G.*exp(-DeltaE_XG./Vt)+Nc_L./Nc_G.*exp(-DeltaE_LG./Vt));
        DeltaE_GX = Eg_G-Eg_X;
        DeltaE_LX = Eg_L-Eg_X;
        Nc_2 = Nc_X.*(1+Nc_G./Nc_X.*exp(-DeltaE_GX./Vt)+Nc_L./Nc_X.*exp(-DeltaE_LX./Vt));
        %
        indX=find(Eg_G>Eg_X); % index where band X is dominant
        Nc = Nc_1;
        Nc(indX) = Nc_2(indX);
        meffn = (Nc./(2.*(2.*pi.*kB.*T.*m0).^(3/2)/(h.^3).*1e-6)).^(2/3);
        %
        meffp=0.55+0.26*xmol; % effective hole mass
        Nv=2*(2.*pi.*kB.*T.*m0).^(3/2)/(h.^3)*1e-6.*meffp.^(3/2)/CarrierNorm; 
        meffp = (Nv./(2.*(2.*pi.*kB.*T.*m0).^(3/2)/(h.^3).*1e-6)).^(2/3);
        %
        % saving DoSs and effective masses
        macro.Nc=Nc;
        macro.Nc_G=Nc_G;
        macro.Nc_X=Nc_X;
        macro.Nc_L=Nc_L;
        %
        macro.meffn=meffn;
        macro.meffn_G=meffn_G;
        macro.meffn_L=meffn_X;
        macro.meffn_X=meffn_L;
        %
        macro.Nv=Nv;
        macro.meffp=meffp;
        %
        %
%        indReg2=find(xmol>0.45);
        indReg2=find(xmol>0.425);
        
        %%%%% Electron affinity (eV)
        
       if mesh.Qc==0


        affinity=4.07-1.1.*xmol;
        affinity(indReg2)=3.64-0.14.*xmol(indReg2);
        %
        % Temperature dependence from 2005 Adachi (?????)
         Eaff=204;
         Daff=2.75e-4;
         affinity=affinity+Daff.*(Tg.^2./(Tg+Eaff)-300^2/(300+Eaff)); 
        else 
         affinity = -mesh.Qc*(Eg-Eg_GaAs)+4.07;
        end
             
        macro.affinity=affinity;
        %
        %%%%%%%%%% Triangles-dependent quantities %%%%%%%%%%
        %-- Dielectric constant
        macro.epsrxx=12.90-2.84.*xmol;
        %
        %-- Low-field mobility (cm^2/s/V)
        % Electron low-field mobility
%        mobnint = 8000-22000.*xmol+10000.*xmol.^2;
%        mobnint(indReg2)= -255+1160*xmol(indReg2)-720*xmol(indReg2).^2;
% new Roland        
        mobnint = 8000-24000.*xmol+13000.*xmol.^2;
        mobnint(indReg2)= 1200*(xmol(indReg2)-.45).^2+148;
        
        %
        % Temperature dependence from 2005 Adachi, p. 325
        mobnint=mobnint.*(T300./T).^mesh.ExpE*mesh.FatMob;
        macro.mobnint=mobnint;
        %
        % Hole low-field mobility
         mobpint1=400-700*xmol+450*xmol.^2;  %Calciati?  Sentaurus
         mobpint2=370-970*xmol+740*xmol.^2;  % Joffe
         mobpint3=400-775*xmol+535*xmol.^2;   % Roland         
         mobpint=mobpint3;
        %
        % Temperature dependence from 2005 Adachi, p. 335
        mobpint=mobpint.*(T300./T).^mesh.ExpH*mesh.FatMob;
        macro.mobpint=mobpint;
        
        
    case 'GaAsP' 
        %%
        'GaAsP material parameters!'%, pausak
        %-- Quantum-corrected model parameters
        %
        macro.tauscatn = 1e-12*onesxmol; % quantum well capture time, s
        macro.tauscatp = 1e-12*onesxmol; % quantum well escape time, s
        %
        macro.tauscatn = mesh.tausE*onesxmol; % quantum well capture time, s
        macro.tauscatp = mesh.tausH*onesxmol; % quantum well escape time, s
        %
        %-- Incomplete ionization parameters
        T300=mesh.T300;
        [Ed,Ea]=ActivEner(xmol);
        macro.DeltaEa=1e-3*Ea.*(1+mesh.CTemp_Ion*(T-T300)./T300); % (eV)
        macro.DeltaEd=1e-3*Ed.*(1+mesh.CTemp_Ion*(T-T300)./T300); % (eV)

%        macro.DeltaEa=1e-3*25.*(1+mesh.CTemp_Ion*(T-300)./300); % (eV)
%        macro.DeltaEd=1e-3*6.*(1+mesh.CTemp_Ion*(T-300)./300); % (eV)        
        
        %-- Recombination models parameters
        macro.Etrap = 0*onesxmol; % (eV)
        % % SRH (Ioffe, xmol=0.1)
        macro.taun = mesh.taun*onesxmol.*(T300./T).^mesh.tauExp; % (s)
        macro.taup = mesh.taup*onesxmol.*(T300./T).^mesh.tauExp; % (s)
        macro.taunQW = mesh.taunQW*onesxmol.*(T300./T).^mesh.tauExp; % (s)
        macro.taupQW = mesh.taupQW*onesxmol.*(T300./T).^mesh.tauExp; % (s)
        % SRH (Michalzik book, chapter 3, Debernardi)
        % macro.taun = 5e-9*onesxmol.*(300./T); % (s)
        % macro.taup = 5e-9*onesxmol.*(300./T); % (s)
        %
        %- Radiative (Ioffe)
        macro.brad = mesh.Brad*onesxmol*CarrierNorm; % (cm^3/s)
        %
        %- Auger (Ioffe, xmol=0.1)
%         macro.Cnnp = 1.2e-31*onesxmol; % (cm^6/s)
%         macro.Cppn = 8.5e-31*onesxmol; % (cm^6/s)
%	'keyb Auger', keyboard
 %       [Cn,Cp]=fun_auger(xmol,T-mesh.T0);
 %       macro.Cnnp =Cn*10;
%	macro.Cppn =Cp*10;
	
        macro.Cnnp =mesh.CN_Auger*exp(mesh.CTemp_Auger*(T-T300)./100)*1e-30*CarrierNorm^2;
        %macro.Cnnp =mesh.CN_Auger*abs(1+mesh.CTemp_Auger*(T-mesh.T0)./300)*1e-30;
        macro.Cppn =mesh.FatNP_Auger.*macro.Cnnp;
        % Auger (Michalzik book, chapter 3, Debernardi)
        % macro.Cnnp = 3.5e-30*onesxmol; % (cm^6/s)
        % macro.Cppn = 3.5e-30*onesxmol; % (cm^6/s)
        %
        %%%%%%%%%% Nodes-dependent quantities %%%%%%%%%%
        %-- Bandgap (eV) (1985Adachi_JAP, eqs. (20)-(22))
        %-- Temperature dependent bandgap (eV) (1976Aspnes_PRB,
        %   http://www.ioffe.ru/SVA/NSM/Semicond/AlGaAs/bandstr.html)
        Eg0_G=1.519+1.155.*xmol+0.370.*xmol.^2;
        alpha_G=5.41e-4; beta_G=204;
        Eg_G=Eg0_G-(alpha_G.*Tg.^2)./(beta_G+Tg);
        %
        Eg0_X=1.981+0.124.*xmol+0.144.*xmol.^2;
        alpha_X=4.6e-4; beta_X=204;
        Eg_X=Eg0_X-(alpha_X.*Tg.^2)./(beta_X+Tg);
        %
        Eg0_L=1.815+0.69.*xmol; % mind the error in Ioffe!
        alpha_L=6.05e-4; beta_L=204;
        Eg_L=Eg0_L-(alpha_L.*Tg.^2)./(beta_L+Tg);
        
        %
%' qui Eg', keyboard        
        Eg = min(Eg_G,Eg_X);
        %
        % Calculating Eg in GaAs
        Eg0_G_GaAs=1.519;
        alpha_G=5.41e-4; beta_G=204;
        Eg_G_GaAs=Eg0_G_GaAs-(alpha_G.*Tg.^2)./(beta_G+Tg);
        %
        Eg0_X_GaAs=1.981;
        alpha_X=4.6e-4; beta_X=204;
        Eg_X_GaAs=Eg0_X_GaAs-(alpha_X.*Tg.^2)./(beta_X+Tg);
        %
        Eg0_L_GaAs=1.815; % mind the error in Ioffe!
        alpha_L=6.05e-4; beta_L=204;
        Eg_L_GaAs=Eg0_L_GaAs-(alpha_L.*Tg.^2)./(beta_L+Tg);
        %
        Eg_GaAs = min(Eg_G_GaAs,Eg_X_GaAs);

        %
        % saving bandgaps
        macro.Eg_G = Eg_G;
        macro.Eg_X = Eg_X;
        macro.Eg_L = Eg_L;
        macro.Eg = Eg;
        %
        %-- Effective DoS masses (1985Adachi_JAP, Table II)
        meffn_G = 0.067+0.083.*xmol;
        meffn_X = 0.850-0.140.*xmol;
        meffn_L = 0.560+0.100.*xmol;
        %
        % DoS
        Nc_G = 2*(2*pi*kB*T*m0).^(3/2)/(h.^3).*1e-6.*meffn_G.^(3/2)/CarrierNorm;
        Nc_X = 2*(2*pi*kB*T*m0).^(3/2)/(h.^3).*1e-6.*meffn_X.^(3/2)/CarrierNorm;
        Nc_L = 2*(2*pi*kB*T*m0).^(3/2)/(h.^3).*1e-6.*meffn_L.^(3/2)/CarrierNorm;
        %
        % Equivalent DoS
        Vt=kB.*T./qel;
        DeltaE_XG = Eg_X-Eg_G;
        DeltaE_LG = Eg_L-Eg_G;
        Nc_1 = Nc_G.*(1+Nc_X./Nc_G.*exp(-DeltaE_XG./Vt)+Nc_L./Nc_G.*exp(-DeltaE_LG./Vt));
        DeltaE_GX = Eg_G-Eg_X;
        DeltaE_LX = Eg_L-Eg_X;
        Nc_2 = Nc_X.*(1+Nc_G./Nc_X.*exp(-DeltaE_GX./Vt)+Nc_L./Nc_X.*exp(-DeltaE_LX./Vt));
        %
        indX=find(Eg_G>Eg_X); % index where band X is dominant
        Nc = Nc_1;
        Nc(indX) = Nc_2(indX);
        meffn = (Nc./(2.*(2.*pi.*kB.*T.*m0).^(3/2)/(h.^3).*1e-6)).^(2/3);
        %
        meffp=0.55+0.26*xmol; % effective hole mass
        Nv=2*(2.*pi.*kB.*T.*m0).^(3/2)/(h.^3)*1e-6.*meffp.^(3/2)/CarrierNorm; 
        meffp = (Nv./(2.*(2.*pi.*kB.*T.*m0).^(3/2)/(h.^3).*1e-6)).^(2/3);
        %
        % saving DoSs and effective masses
        macro.Nc=Nc;
        macro.Nc_G=Nc_G;
        macro.Nc_X=Nc_X;
        macro.Nc_L=Nc_L;
        %
        macro.meffn=meffn;
        macro.meffn_G=meffn_G;
        macro.meffn_L=meffn_X;
        macro.meffn_X=meffn_L;
        %
        macro.Nv=Nv;
        macro.meffp=meffp;
        %
        %
%        indReg2=find(xmol>0.45);
        indReg2=find(xmol>0.425);
        
        %%%%% Electron affinity (eV)
        
       if mesh.Qc==0


        affinity=4.07-1.1.*xmol;
        affinity(indReg2)=3.64-0.14.*xmol(indReg2);
        %
        % Temperature dependence from 2005 Adachi (?????)
         Eaff=204;
         Daff=2.75e-4;
         affinity=affinity+Daff.*(Tg.^2./(Tg+Eaff)-300^2/(300+Eaff)); 
        else 
         affinity = -mesh.Qc*(Eg-Eg_GaAs)+4.07;
        end
             
        macro.affinity=affinity;
        %
        %%%%%%%%%% Triangles-dependent quantities %%%%%%%%%%
        %-- Dielectric constant
        macro.epsrxx=12.90-2.84.*xmol;
        %
        %-- Low-field mobility (cm^2/s/V)
        % Electron low-field mobility
%        mobnint = 8000-22000.*xmol+10000.*xmol.^2;
%        mobnint(indReg2)= -255+1160*xmol(indReg2)-720*xmol(indReg2).^2;
% new Roland        
        mobnint = 8000-24000.*xmol+13000.*xmol.^2;
        mobnint(indReg2)= 1200*(xmol(indReg2)-.45).^2+148;
        
        %
        % Temperature dependence from 2005 Adachi, p. 325
        mobnint=mobnint.*(T300./T).^mesh.ExpE*mesh.FatMob;
        macro.mobnint=mobnint;
        %
        % Hole low-field mobility
         mobpint1=400-700*xmol+450*xmol.^2;  %Calciati?  Sentaurus
         mobpint2=370-970*xmol+740*xmol.^2;  % Joffe
         mobpint3=400-775*xmol+535*xmol.^2;   % Roland         
         mobpint=mobpint3;
        %
        % Temperature dependence from 2005 Adachi, p. 335
        mobpint=mobpint.*(T300./T).^mesh.ExpH*mesh.FatMob;
        macro.mobpint=mobpint;        
        
        %%
%         %------ Quantities indepdendent from the molar fraction ------%
%         %
%         %-- Scattering time for quantum correction
%         macro.tauscatn = mp('3e-12')*ones(size(xmol)); % (s)
%         macro.tauscatp = mp('1e-10')*ones(size(xmol)); % (s)
%         %
%         %-- Incomplete ionization parameters
%         macro.DeltaEa = mp('0.026')*ones(size(xmol)); % (eV)
%         macro.DeltaEd = mp('0.005')*ones(size(xmol)); % (eV)
%         %
%         %-- Recombination models parameters
%         %- Shockley-Read-Hall (Ioffe)
%         macro.Etrap = mp('0')*ones(size(xmol)); % (eV)
%         macro.taun = mp('5e-9')*ones(size(xmol)); % (s)
%         macro.taup = mp('20e-9')*ones(size(xmol)); % (s)
%         %
%         %- Radiative (Ioffe)
%         macro.brad = mp('1.8e-10')*ones(size(xmol)); % (cm^3/s)
%         %
%         %- Auger (Ioffe)
%         macro.Cnnp = mp('1.2e-30')*ones(size(xmol)); % (cm^6/s)
%         macro.Cppn = mp('8.5e-30')*ones(size(xmol)); % (cm^6/s)
%         %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%------ Quantities depending on the molar fraction ------%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %-- VEGARD'S LAW: p(x,y) = xy*p(GaAs)+x(1-y)*p(GaP)+...
%         %                          (1-x)y*p(InAs)+(1-x)(1-y)p(InP)
%         %
%         %-- Dielectric constant (Ioffe @300 K, static)
%         epsr_GaAs = mp('12.9');
%         epsr_GaP = mp('11.1');
%         epsr_InAs = mp('15.15');
%         epsr_InP = mp('12.5');
%         epsr = xmol.*ymol.*epsr_GaAs + xmol.*(1-ymol).*epsr_GaP + ...
%               (1-xmol).*ymol.*epsr_InAs + (1-xmol).*(1-ymol).*epsr_InP;
%         macro.epsr=epsr;
%         %
%         %-- Energy Gap (Ioffe, @300 K; [eV]) 
%         Eg_GaAs = mp('1.424'); 
%         Eg_GaP = mp('2.26');
%         Eg_InAs = mp('0.354'); 
%         Eg_InP = mp('1.344');
%         Eg = xmol.*ymol.*Eg_GaAs + xmol.*(1-ymol).*Eg_GaP + ...
%               (1-xmol).*ymol.*Eg_InAs + (1-xmol).*(1-ymol).*Eg_InP;
%         macro.Eg = Eg;
%         %
%         %-- Affinity (Affinity Rule, @300 K; [eV])
%         CBO = mp('0.65');
%         affinity_InAs = mp('4.9');
%         affinity = affinity_InAs - CBO*(Eg-Eg_InAs);
%         macro.affinity = affinity;
%         %
%         %-- Effective mass (Ioffe, @300 K; [*m0])
%         meffn_GaAs = mp('0.063');
%         meffn_GaP = mp('0.378'); 
%         meffn_InAs = mp('0.023'); 
%         meffn_InP = mp('0.08');
%         meffn = xmol.*ymol.*meffn_GaAs + xmol.*(1-ymol).*meffn_GaP + ...
%                 (1-xmol).*ymol.*meffn_InAs + (1-xmol).*(1-ymol).*meffn_InP;
%         macro.meffn = meffn;
%         %
%         meffp_GaAs = mp('0.51');
%         meffp_GaP = mp('0.79'); 
%         meffp_InAs = mp('0.41'); 
%         meffp_InP = mp('0.6');
%         meffp = xmol.*ymol.*meffp_GaAs + xmol.*(1-ymol).*meffp_GaP + ...
%                 (1-xmol).*ymol.*meffp_InAs + (1-xmol).*(1-ymol).*meffp_InP;
%         macro.meffp = meffp;
%         % Reduced effective mass
%         meff_red = (meffn.*meffp)/(meffn+meffp);
%         macro.meff_red = meff_red;
%         %
%         % Kane parameters for BTBT tunnelling
%         A = qel^2*sqrt(meff_red*m0)./(mp('18*pi')*(h/(mp('2*pi')))^2*sqrt(Eg*qel))*('1e-2');
%         macro.A = A;        % [V^-2*cm^-1*s^-1]
%         B = mp('pi')*sqrt((Eg.*qel).^3).*sqrt(meff_red.*m0)./(2*qel*h/(mp('2*pi')))*('1e-2');
%         macro.B = B;        % [V/cm]
%         %
%         %
%         % Density of States (DoS)
%         macro.Nc=2*(mp('2*pi')*kB*T*m0).^mp('3/2')/(h^3).*mp('1e-6').*macro.meffn.^mp('3/2'); 
%         macro.Nv=2*(mp('2*pi')*kB*T*m0).^mp('3/2')/(h^3).*mp('1e-6')*macro.meffp.^mp('3/2');
%         %
%         %-- Low-field mobility (Ioffe + Vegard, @300 K, cm^2/V*s)
%         mobnint_GaAs = mp('8500');
%         mobnint_GaP = mp('250');
%         mobnint_InAs = mp('4e4');
%         mobnint_InP = mp('5400');
%         mobnint = xmol.*ymol.*mobnint_GaAs + xmol.*(1-ymol).*mobnint_GaP + ...
%                 (1-xmol).*ymol.*mobnint_InAs + (1-xmol).*(1-ymol).*mobnint_InP;
%         macro.mobnint = mobnint;
%         %
%         mobpint_GaAs = mp('400');
%         mobpint_GaP = mp('150');
%         mobpint_InAs = mp('500');
%         mobpint_InP = mp('200');
%         mobpint = xmol.*ymol.*mobpint_GaAs + xmol.*(1-ymol).*mobpint_GaP + ...
%                 (1-xmol).*ymol.*mobpint_InAs + (1-xmol).*(1-ymol).*mobpint_InP;
%         macro.mobpint = mobpint;
    case 'InGaAs' 
        %%
        'InGaAs material parameters!'%, pausak
        %-- Quantum-corrected model parameters
        %
        macro.tauscatn = 1e-12*onesxmol; % quantum well capture time, s
        macro.tauscatp = 1e-12*onesxmol; % quantum well escape time, s
        %
        macro.tauscatn = mesh.tausE*onesxmol; % quantum well capture time, s
        macro.tauscatp = mesh.tausH*onesxmol; % quantum well escape time, s
        %
        %-- Incomplete ionization parameters
        T300=mesh.T300;
        [Ed,Ea]=ActivEner(xmol);
        macro.DeltaEa=1e-3*Ea.*(1+mesh.CTemp_Ion*(T-T300)./T300); % (eV)
        macro.DeltaEd=1e-3*Ed.*(1+mesh.CTemp_Ion*(T-T300)./T300); % (eV)

%        macro.DeltaEa=1e-3*25.*(1+mesh.CTemp_Ion*(T-300)./300); % (eV)
%        macro.DeltaEd=1e-3*6.*(1+mesh.CTemp_Ion*(T-300)./300); % (eV)        
        
        %-- Recombination models parameters
        macro.Etrap = 0*onesxmol; % (eV)
        % % SRH (Ioffe, xmol=0.1)
        macro.taun = mesh.taun*onesxmol.*(T300./T).^mesh.tauExp; % (s)
        macro.taup = mesh.taup*onesxmol.*(T300./T).^mesh.tauExp; % (s)
        macro.taunQW = mesh.taunQW*onesxmol.*(T300./T).^mesh.tauExp; % (s)
        macro.taupQW = mesh.taupQW*onesxmol.*(T300./T).^mesh.tauExp; % (s)
        % SRH (Michalzik book, chapter 3, Debernardi)
        % macro.taun = 5e-9*onesxmol.*(300./T); % (s)
        % macro.taup = 5e-9*onesxmol.*(300./T); % (s)
        %
        %- Radiative (Ioffe)
        macro.brad = mesh.Brad*onesxmol*CarrierNorm; % (cm^3/s)
        %
        %- Auger (Ioffe, xmol=0.1)
%         macro.Cnnp = 1.2e-31*onesxmol; % (cm^6/s)
%         macro.Cppn = 8.5e-31*onesxmol; % (cm^6/s)
%	'keyb Auger', keyboard
 %       [Cn,Cp]=fun_auger(xmol,T-mesh.T0);
 %       macro.Cnnp =Cn*10;
%	macro.Cppn =Cp*10;
	
        macro.Cnnp =mesh.CN_Auger*exp(mesh.CTemp_Auger*(T-T300)./100)*1e-30*CarrierNorm^2;
        %macro.Cnnp =mesh.CN_Auger*abs(1+mesh.CTemp_Auger*(T-mesh.T0)./300)*1e-30;
        macro.Cppn =mesh.FatNP_Auger.*macro.Cnnp;
        % Auger (Michalzik book, chapter 3, Debernardi)
        % macro.Cnnp = 3.5e-30*onesxmol; % (cm^6/s)
        % macro.Cppn = 3.5e-30*onesxmol; % (cm^6/s)
        %
        %%%%%%%%%% Nodes-dependent quantities %%%%%%%%%%
        %-- Bandgap (eV) (1985Adachi_JAP, eqs. (20)-(22))
        %-- Temperature dependent bandgap (eV) (1976Aspnes_PRB,
        %   http://www.ioffe.ru/SVA/NSM/Semicond/AlGaAs/bandstr.html)
        Eg0_G=1.519+1.155.*xmol+0.370.*xmol.^2;
        alpha_G=5.41e-4; beta_G=204;
        Eg_G=Eg0_G-(alpha_G.*Tg.^2)./(beta_G+Tg);
        %
        Eg0_X=1.981+0.124.*xmol+0.144.*xmol.^2;
        alpha_X=4.6e-4; beta_X=204;
        Eg_X=Eg0_X-(alpha_X.*Tg.^2)./(beta_X+Tg);
        %
        Eg0_L=1.815+0.69.*xmol; % mind the error in Ioffe!
        alpha_L=6.05e-4; beta_L=204;
        Eg_L=Eg0_L-(alpha_L.*Tg.^2)./(beta_L+Tg);
        
        %
%' qui Eg', keyboard        
        Eg = min(Eg_G,Eg_X);
        %
        % Calculating Eg in GaAs
        Eg0_G_GaAs=1.519;
        alpha_G=5.41e-4; beta_G=204;
        Eg_G_GaAs=Eg0_G_GaAs-(alpha_G.*Tg.^2)./(beta_G+Tg);
        %
        Eg0_X_GaAs=1.981;
        alpha_X=4.6e-4; beta_X=204;
        Eg_X_GaAs=Eg0_X_GaAs-(alpha_X.*Tg.^2)./(beta_X+Tg);
        %
        Eg0_L_GaAs=1.815; % mind the error in Ioffe!
        alpha_L=6.05e-4; beta_L=204;
        Eg_L_GaAs=Eg0_L_GaAs-(alpha_L.*Tg.^2)./(beta_L+Tg);
        %
        Eg_GaAs = min(Eg_G_GaAs,Eg_X_GaAs);

        %
        % saving bandgaps
        macro.Eg_G = Eg_G;
        macro.Eg_X = Eg_X;
        macro.Eg_L = Eg_L;
        macro.Eg = Eg;
        %
        %-- Effective DoS masses (1985Adachi_JAP, Table II)
        meffn_G = 0.067+0.083.*xmol;
        meffn_X = 0.850-0.140.*xmol;
        meffn_L = 0.560+0.100.*xmol;
        %
        % DoS
        Nc_G = 2*(2*pi*kB*T*m0).^(3/2)/(h.^3).*1e-6.*meffn_G.^(3/2)/CarrierNorm;
        Nc_X = 2*(2*pi*kB*T*m0).^(3/2)/(h.^3).*1e-6.*meffn_X.^(3/2)/CarrierNorm;
        Nc_L = 2*(2*pi*kB*T*m0).^(3/2)/(h.^3).*1e-6.*meffn_L.^(3/2)/CarrierNorm;
        %
        % Equivalent DoS
        Vt=kB.*T./qel;
        DeltaE_XG = Eg_X-Eg_G;
        DeltaE_LG = Eg_L-Eg_G;
        Nc_1 = Nc_G.*(1+Nc_X./Nc_G.*exp(-DeltaE_XG./Vt)+Nc_L./Nc_G.*exp(-DeltaE_LG./Vt));
        DeltaE_GX = Eg_G-Eg_X;
        DeltaE_LX = Eg_L-Eg_X;
        Nc_2 = Nc_X.*(1+Nc_G./Nc_X.*exp(-DeltaE_GX./Vt)+Nc_L./Nc_X.*exp(-DeltaE_LX./Vt));
        %
        indX=find(Eg_G>Eg_X); % index where band X is dominant
        Nc = Nc_1;
        Nc(indX) = Nc_2(indX);
        meffn = (Nc./(2.*(2.*pi.*kB.*T.*m0).^(3/2)/(h.^3).*1e-6)).^(2/3);
        %
        meffp=0.55+0.26*xmol; % effective hole mass
        Nv=2*(2.*pi.*kB.*T.*m0).^(3/2)/(h.^3)*1e-6.*meffp.^(3/2)/CarrierNorm; 
        meffp = (Nv./(2.*(2.*pi.*kB.*T.*m0).^(3/2)/(h.^3).*1e-6)).^(2/3);
        %
        % saving DoSs and effective masses
        macro.Nc=Nc;
        macro.Nc_G=Nc_G;
        macro.Nc_X=Nc_X;
        macro.Nc_L=Nc_L;
        %
        macro.meffn=meffn;
        macro.meffn_G=meffn_G;
        macro.meffn_L=meffn_X;
        macro.meffn_X=meffn_L;
        %
        macro.Nv=Nv;
        macro.meffp=meffp;
        %
        %
%        indReg2=find(xmol>0.45);
        indReg2=find(xmol>0.425);
        
        %%%%% Electron affinity (eV)
        
       if mesh.Qc==0


        affinity=4.07-1.1.*xmol;
        affinity(indReg2)=3.64-0.14.*xmol(indReg2);
        %
        % Temperature dependence from 2005 Adachi (?????)
         Eaff=204;
         Daff=2.75e-4;
         affinity=affinity+Daff.*(Tg.^2./(Tg+Eaff)-300^2/(300+Eaff)); 
        else 
         affinity = -mesh.Qc*(Eg-Eg_GaAs)+4.07;
        end
             
        macro.affinity=affinity;
        %
        %%%%%%%%%% Triangles-dependent quantities %%%%%%%%%%
        %-- Dielectric constant
        macro.epsrxx=12.90-2.84.*xmol;
        %
        %-- Low-field mobility (cm^2/s/V)
        % Electron low-field mobility
%        mobnint = 8000-22000.*xmol+10000.*xmol.^2;
%        mobnint(indReg2)= -255+1160*xmol(indReg2)-720*xmol(indReg2).^2;
% new Roland        
        mobnint = 8000-24000.*xmol+13000.*xmol.^2;
        mobnint(indReg2)= 1200*(xmol(indReg2)-.45).^2+148;
        
        %
        % Temperature dependence from 2005 Adachi, p. 325
        mobnint=mobnint.*(T300./T).^mesh.ExpE*mesh.FatMob;
        macro.mobnint=mobnint;
        %
        % Hole low-field mobility
         mobpint1=400-700*xmol+450*xmol.^2;  %Calciati?  Sentaurus
         mobpint2=370-970*xmol+740*xmol.^2;  % Joffe
         mobpint3=400-775*xmol+535*xmol.^2;   % Roland         
         mobpint=mobpint3;
        %
        % Temperature dependence from 2005 Adachi, p. 335
        mobpint=mobpint.*(T300./T).^mesh.ExpH*mesh.FatMob;
        macro.mobpint=mobpint;        
        
        %%
%         mode.tauncap=[mp('0.3e-12'); mp('1e-12'); mp('1e-12')];%.*mtau(I_mtau)];         % electron capture time for each state, (s) B->WL; WL->ES ES->GS
%         mode.taupcap=[mp('0.1e-12'); mp('0.1e-12'); mp('0.1e-12')];   % hole capture time for each state, (s) B->WL; WL->ES ES->GS
%         
%         %------ Quantities indepdendent from the molar fraction ------%
%         %
%         %-- Scattering time for quantum correction
%         macro.tauscatn=[mp('3') mp('1') mp('1')]*mp('1e-12'); % WL, ES, GS (s)
%         macro.tauscatp=[mp('100') mp('0.1') mp('0.1')]*mp('1e-12'); % WL, ES, GS
%         %
%         %-- Incomplete ionization parameters
%         macro.DeltaEa=mp('0.026')*ones(size(xmol)); % (eV)
%         macro.DeltaEd=mp('0.005')*ones(size(xmol)); % (eV)
%         %
%         %-- Recombination models parameters
%         %- Shockley-Read-Hall (Ioffe, xmol=0.1)
%         macro.Etrap = mp('0')*ones(size(xmol)); % (eV)
%         macro.taun = mp('5e-9')*ones(size(xmol)); % (s)
%         macro.taup = mp('20e-9')*ones(size(xmol)); % (s)
%         %
%         %- Radiative (Ioffe)
%         macro.brad = mp('6.22e-10')*ones(size(xmol)); % (cm^3/s)
%         %
%         %- Auger (Ioffe, xmol=0.1)
%         macro.Cnnp = mp('1.2e-31')*ones(size(xmol)); % (cm^6/s)
%         macro.Cppn = mp('8.5e-31')*ones(size(xmol)); % (cm^6/s)
%         %
%         %------ Quantities depending on the molar fraction ------%
%         %
%         %-- Dielectric constant (Ioffe @300 K, static)
%         macro.epsr=mp('15.10')-mp('2.87').*xmol+mp('0.67').*xmol.^2;
%         %
%         %-- Indexes for molar fraction in Gamma and X bands
%         %-- Low-field mobility (cm^2/s/V)
%         % Electron mobility
%         mobnint = (mp('40')-mp('80.7').*xmol+mp('49.2').*xmol.^2)*mp('1e3');
%         macro.mobnint=mobnint;
%         %
%         % Hole mobility
%         macro.mobpint=mp('350')*ones(size(xmol));
%         %
%         %-- Bandgap (eV)
%         Eg = mp('0.36')+mp('0.63').*xmol+mp('0.43').*xmol.^2;
%         macro.Eg=Eg;
%         %
%         %-- Effective mass
%         macro.meffn=mp('0.023')+mp('0.037').*xmol+mp('0.003').*xmol.^2; % effective electron mass
%         macro.meffp=mp('0.41')+mp('0.1')*xmol; % effective hole mass
%         %
%         % Density of States (DoS)
%         macro.Nc=2*(2*mp('pi')*kB*T*m0).^(mp('3/2'))/(h.^3).*mp('1e-6').*macro.meffn.^(mp('3/2'));
%         macro.Nv=2*(2*mp('pi')*kB*T*m0).^(mp('3/2'))/(h.^3).*mp('1e-6').*macro.meffp.^(mp('3/2'));
%         %
%         %%%%% Electron affinity (eV)
%         affinity=mp('4.9')-mp('0.83').*xmol; % affinity, eV
%         macro.affinity=affinity;
%         
    case 'AlOx' % oxide of AlAs (AlGaAs,xmol=0)
        %%
        %-- "AlOx": "x" stands for "unknown fraction of oxygen"
        %-- n in [1.6,1.7]; use n=1.65; epsr=n^2.
        epsr=1.65^2*onesxmol;
        macro.epsrxx=epsr;

    case 'Polyamide' % passivation material
        %%
        epsr=3.5*onesxmol;
        macro.epsrxx=epsr;

    case 'vacuum' % free-space
        %%
        epsr=1*onesxmol;
        macro.epsrxx=epsr;
        macro.DeltaEa=0; macro.DeltaEd=0;
        
    case 'Au' % gold
        %%
        epsr=1*onesxmol;
        macro.epsrxx=epsr;

    otherwise
        error('Unknown material'),
end
% *********************************************************************************************100
