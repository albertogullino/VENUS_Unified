function [G,Es,Dep,Rsp,eD1,hD1,EFcv,EFvv,Ren,EFsplit]=pf_functionFit(eDensityv,hDensityv,lambda,Tvet,indT,Ban,mesh,mode,Pargain)

% pf_functionFit : "parfor" (pf) function with "infittimento" (Fit).

global PPEc PPEv
global PPEcV PPVccV PPEvV PPVvvV
% global PPVcc PPVvv % non-screened potentials
global PPWcc PPWvv % screened potentials
global C0 V0

global qel T EFc EFv
global kBTev K1

if(mode.LUT==0)
    global EFsave ECVtotsave
end

T=Tvet(indT);

PPEcV=Ban.PPEcV;
PPEvV=Ban.PPEvV;
PPVccV=Ban.PPVccV;
PPVvvV=Ban.PPVvvV;

    FatRen=1;

if isfield(mode,'iren')
    iren=mode.iren;
    FatRen=mode.FatRen;
else
    iren=1;
end

OPTIONS = optimset('TolX',1e-9,'MaxIter',1000);
Nb=Ban.Nb;
gammak=mesh.gammak;

kdkf=Ban.kdkf;
wkrho=Ban.wkrho;
wktheta=Ban.wktheta;
VKRHO=Ban.VKRHO;
VKTHETA=Ban.VKTHETA;
nbc=mesh.ncb;
nbv=mesh.nvb;

kgrid=Ban.kgrid;
kgridf=Ban.kgridf;
nkf=length(kgridf);
SBC=Ban.SBC;
SBV=Ban.SBV;
SBCf0=Ban.SBCf;
SBVf0=Ban.SBVf;

C0=mesh.C0;
V0=mesh.V0;
iline0=mode.iline;
iline=iline0;
iland=-1;

if strcmp(iline0(1:3),'Lan')
 iline='lorentzian';
% iline='nMark';
 if strcmp(iline0(end-1:end),'ed')
  iland=0;
 else
  iland=1;
 end
 iLaNew=0; 
 if strcmp(iline0(end-2),'i')
  iLaNew=1;  
 end
 
end

tnm=mesh.tnm;  %s  for non markovian effects
ganm=1/tnm;  %in 1/s
ganm0=ganm;

s_LoadConstants

omegavet=2*pi*c_light./lambda;
const_HF = 1/hbar * 1/(2*pi)^2;

k0=Nb*2*pi./lambda;
vic=1:mesh.ncb;
viv=1:mesh.nvb;

%'Pargain', keyboard
kBTev=Pargain(indT).kBTev;
kBT=Pargain(indT).kBT;
M2d=Pargain(indT).M2d;
Cost_Rsp=Pargain(indT).Cost_Rsp;
M2esd=Pargain(indT).M2esd;
ECV_tot0=Pargain(indT).ECV_tot0;


if mode.ieh_equal==1
    eDensityv=1;
end

% Initialization of variables to save memory and computation time
EFvv=zeros(1,length(hDensityv));
hD1=zeros(1,length(hDensityv));
EFcv=zeros(1,length(hDensityv));
eD1=zeros(1,length(hDensityv));
Ren=zeros(length(hDensityv),nkf);
Rsp=zeros(length(hDensityv),length(eDensityv));
G=zeros(length(hDensityv),length(eDensityv),length(omegavet));
Dep=zeros(length(hDensityv),length(eDensityv),length(omegavet));
Es=zeros(length(hDensityv),length(eDensityv),length(omegavet));

EFnsave=zeros(length(eDensityv),length(omegavet));

for idh=1:length(hDensityv)
    hDensity=hDensityv(idh); % 1/cm^2
    
    % Notice that the quasi-Fermi level for holes, EFv, is computed
    % assuming a potential well equal to that of electrons. With this
    % convention, low hole densities correspond to negative EFv values,
    % whereas high densities to positive EFv. This is completely analog to
    % electrons  / conduction band
    [EFv] = fminbnd(@(EFv) f_rho(kgrid,SBV,EFv,kBTev,hDensity,-V0,mode),-mesh.Eg,-2*V0,OPTIONS);
    
    for ide=1:length(eDensityv)
        
        if mode.ieh_equal==0
            eDensity=eDensityv(ide); % 1/cm^2
        else
            eDensity=hDensity; % 1/cm^2
        end
        
        [EFc] = fminbnd(@(EFc) f_rho(kgrid,SBC,EFc,kBTev,eDensity,C0+mesh.Eg,mode), 0,mesh.Eg*2,OPTIONS);
        
        %==================================================================
        % Many-body effects: gap renormalization
        %==================================================================
        DelGap = zeros(nbc*nbv,nkf);
        if iren==1
            
%            'IREN', keyboard
            % Computing Lindhard epsilon for plasma screening
            % Veprek PhD, (5.26)
            [epsfun_PD] = f_EvalLindhardScreening(kgrid,wkrho,wktheta,VKRHO,VKTHETA,vic,viv);
            
            % Correction including measured gap, Veprek PhD, (5.30).
            eps_fatt=1./epsfun_PD-1;

%            eps_fatt=-1;
            
            
            
            % Evaluating and fitting screened potentials, Veprek PhD,
            % (5.25)
            for ic = 1:mesh.ncb
                Vcc=ppval(PPVccV{ic},kgrid);
                Wcc = log(-Vcc.*eps_fatt);
                PPWccV{ic} = spline(log(kgrid),Wcc);
            end
            for iv = 1:mesh.nvb
                Vvv=ppval(PPVvvV{iv},kgrid);
                Wvv = log(-Vvv.*eps_fatt);
                PPWvvV{iv} = spline(log(kgrid),Wvv);
            end
            
            % Gap renormalization can be obtained by integrating the
            % screened diagonal Hartree-Fock contribution
            
            % Renormalization for electrons
            RenGapE = @(k,theta) HF_diagE(k,theta);
            for ic = 1:mesh.ncb
                PPWcc=PPWccV{ic};
                PPEc=PPEcV{ic};
                for iq = 1:length(kgrid)
                    K1 = kgrid(iq);
                    RenormIntegrand=RenGapE(VKRHO,VKTHETA);
                    DelGapIntegrand=wkrho*RenormIntegrand*wktheta';
                    RR = DelGapIntegrand*const_HF;
                    DelEc0(ic,iq) = RR;
                end
                DelEc(ic,:)=spline(kgrid,DelEc0(ic,:),kgridf);
            end
            EFco=EFc;
            
      %      'fine ren E', keyboard
            % Many-body gap renormalization can be seen as a correction to the
            % subbands, but it is more complicated than temperature, as every k
            % point and carrier density has different corrections; this effect is taken
            % into account by correcting quasi-Fermi levels a posteriori. This is a first
            % order correction, not self-consistent (it is like computing just 1 iteration of
            % a self-consistent loop). To perform a real self-consistent calculation, one should
            % re-compute Lindhard screening with the new Fermi levels and the renormalization
            % integrals, for each iteration.
            
            % This DeE is the gap renormalization, due to many-body effect, 
            % shifted by its k=0 value; this is for electrons only.
            % By this way, to compute the total renormalization, the 
            % bidimensional gap should be added.
            % This is useful since now we are taking into account 
            % renormalization in the evaluation of quasi-Fermi levels.
            DeE=(DelEc0-DelEc0(:,1)*ones(1,length(kgrid)))*hbar/qel;
            [EFc] = fminbnd(@(EFc) f_rho(kgrid,SBC+DeE,EFc,kBTev,eDensity,C0+mesh.Eg,mode), EFco,EFco+.1,OPTIONS);
            
            % Renormalization for holes
            RenGapH = @(k,theta) HF_diagH(k,theta);
            for iv = 1:mesh.nvb
                PPWvv=PPWvvV{iv};
                PPEv=PPEvV{iv};
                for iq = 1:length(kgrid)
                    K1 = kgrid(iq);
                    RenormIntegrand=RenGapH(VKRHO,VKTHETA);
                    DelGapIntegrand=wkrho*RenormIntegrand*wktheta';
                    RR = DelGapIntegrand*const_HF;
                    DelEv0(iv,iq) = RR;
                end
                DelEv(iv,:)=spline(kgrid,DelEv0(iv,:),kgridf);
            end %iv
            
            % This DeH is the gap renormalization, due to many-body effect, 
            % shifted by its k=0 value; this is for holes only.
            % by this way, to compute the total renormalization, the 
            % bidimensional gap should be added.
            % This is useful since now we are taking into account 
            % renormalization in the evaluation of quasi-Fermi levels.
            DeH=(DelEv0-DelEv0(:,1)*ones(1,length(kgrid)))*hbar/qel;
            EFvo=EFv;
            [EFv] = fminbnd(@(EFv) f_rho(kgrid,SBV+DeH,EFv,kBTev,hDensity,-V0,mode),EFvo,EFvo+.05,OPTIONS);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            for ic = 1:mesh.ncb
                for iv = 1:mesh.nvb
                    itra=iv+(ic-1)*mesh.nvb;
                    DelGap(itra,:)=DelEv(iv,:)+DelEc(ic,:);
                end
            end
            
%                    'fine ren H', keyboard
            
            Ren(idh,:)=DelGap(1,:);
            
            [eDensity1,hDensity1] = f_charge_buca(mesh,mode,kgrid,SBC,SBV,EFc,EFv,T,C0+mesh.Eg,-V0);
            
%      'dopo REN', keyboard            
            
            if mode.ieh_equal==0
                if ide==1
                    EFvv(idh)=EFv;
                    hD1(idh)=hDensity1;
                end
                if idh==1
                    EFcv(ide)=EFc;
                    eD1(ide)=eDensity1;
                end
            else
                eD1(idh)=eDensity1;
                hD1(idh)=hDensity1;
                EFvv(idh)=EFv;
                EFcv(idh)=EFc;
            end
            
        else
            [eDensity1,hDensity1] = f_charge_buca(mesh,mode,kgrid,SBC,SBV,EFc,EFv,T,C0+mesh.Eg,-V0);
            
            if mode.ieh_equal==0
                if ide==1
                    EFvv(idh)=EFv;
                    hD1(idh)=hDensity1;
                end
                if idh==1
                    EFcv(ide)=EFc;
                    eD1(ide)=eDensity1;
                end
            else
                eD1(idh)=eDensity1;
                hD1(idh)=hDensity1;
                EFvv(idh)=EFv;
                EFcv(idh)=EFc;
            end
            DelGap=0;
            Ren(idh,:)=0;
        end
        
        ECV_tot=DelGap*FatRen+ECV_tot0; % applying gap renormalization
        
        M2=M2d./(hbar*ECV_tot).^2; % dipole momentum for gain / DeltaN
        M20=M2d./(hbar*ECV_tot0).^2; % dipole momentum for zero carriers
        M2es=M2esd./(hbar*ECV_tot).^2; % dipole momentum for spontaneous emission
        
        % Actually nv should be computed differently; what is computed here
        % is "1-nv". Refer for example to 1997 Bava Debernardi report, p. 3
        % These re-definitions of the bands take into account the quasi-Fermi
        % levels corrected due to the many-body gap renormalization.
        if(mode.iren==1)
            DeEf=(DelEc-DelEc(:,1)*ones(1,length(kgridf)))*hbar/qel;
            DeHf=(DelEv-DelEv(:,1)*ones(1,length(kgridf)))*hbar/qel;
        else
            DeHf=0;
            DeEf=0;
        end
        SBCf=SBCf0+DeEf;
        SBVf=SBVf0+DeHf;
        %
        nv=1./(1+exp((SBVf-EFv)/kBTev)); % this is (1-nv) in 1997Bava, (2)-(5)
        nc=1./(1+exp((SBCf-EFc)/kBTev)); % this is nc in 1997Bava, (2)-(5)
        % This is the Fermi factor for gain in 1997Bava, eq. (2)
        Fer=[];
        for ke=1:nbc
            Fer=[Fer; nv-1+ones(nbv,1)*nc(ke,:)];
        end
        
        % This is the Fermi factor for spontaneous emission in 1997Bava,
        % eq. (4)
        Fersp=[];
        for ke=1:nbc
            for kh=1:nbv
                Fersp=[Fersp; nv(kh,:).*nc(ke,:)];
            end
        end
        
        % Computing spontaenous recombination rate: the one integrated in
        % the whole lambda (or energy) range; for this reason, this is out
        % of the energy loop
        lambdaSPn=M2es.*Fersp.*ECV_tot.^3;
        Rspi = Cost_Rsp*sum(lambdaSPn*kdkf);
        Rsp(idh,ide)=Rspi;
        
        for indE=1:length(omegavet)
            omega=omegavet(indE);
            lanm=2*pi*c_light./omega*1e9;
            gintegrand=0;
            
            DeltaE=ECV_tot-omega; % tutto in 1/s
            DeltaV=DeltaE;
            [lm,fi0]=min(abs(DeltaE),[],2);
            for kb=1:length(fi0)
                if fi0(kb)>1 & fi0(kb)<length(kgridf)
                    DeltaE(kb,fi0(kb))=mean(abs(DeltaE(kb,fi0(kb)+[-1 0 1])));
                end
            end
            %
            %==============================================================
            % In the following, a broadening model previously developed by
            % Landsberg in 1966 is proposed. This has been simplified by
            % Martin in 1977Martin_SSC.
            % This model can be thought as some sort of alternative to the
            % non-Markovian lineshape, but it should also include the
            % effect of carrier-carrier scattering. Moreover, it does not
            % exhibit bad oscillations. Still, it is not optimal, as it
            % does not work, in its current implementation, for low
            % densities, as it would provide 0 accordingly to a simplified
            % version of 1987Zielinski_JQE. See also Crosslight manuals,
            % with reference to "Landsberg" or "Landsberg broadening".
            DeltaETemp=Pargain(indT).DeltaE_Temp;
            %
            % The choice of writing E1, EF removing the 2D gap comes from
            % reading 1977Martin_SSC, Fig. 1, where all data (EF, EFe,
            % EFh...) is provided such that each energy is referred to the
            % bottom of its well (EFe with respect to the conduction band
            % edge, EFh with respect to the valence band edge)
            if iland>=0
                E1=ECV_tot*hbar/qel;
                E0=ECV_tot0*hbar/qel+DeltaETemp;
                Eg2D=E1(:,1);
                Eg2D0=E0(:,1);
                EF=EFc+EFv-Eg2D0;
                E1=E1-Eg2D*ones(1,length(kgridf));
                %
                % Coefficients from 1977Martin_SSC, eq. (2)
                a0=1;
                a1=-2.229;
                a2=1.458;
                a3=-0.229;
                
	    y1=32.21;
            y2=-43.82;
            y3=15.07;

           if iLaNew==1
   	    a1=y2/y1;
            a2=y3/y1;
            a3=0;
           end         
           
            x=E1./(EF*ones(1,length(kgridf)));
            Factor=a0+a1.*x+a2.*x.^2+a3.*x.^3;
            fie=find(x>1);
           if iLaNew==1
            Factor(fie)=exp(-5.12.*x(fie).^2+6.35)./y1; 
           end              

            fie=find(x<0);            
            Factor(fie)=1;
            
            
%          'factor;m', keyboard
            cond=x>0 & x<1;        
%            cond=x>0;        
            
%            Flast=Factor.*(cond);     
            Flast=Factor;                     
                
                
                
                % Setting the minimum value of Factor (Flast)
                % instead than setting it to zero.
                fiF=find(Flast<=.01);
                if length(fiF)<prod(size(Flast))
                    Flast(fiF)=.1;
                end
            else
                Flast=1;
            end
            
            gammakS=gammak.*Flast;
%            'ga', keyboard
            %
            %==============================================================
            % This is an empirical estimate of the carrier dependence of
            % the lineshape. This is somehow similar to 1989Asada_JQE, eq.
            % (43), with some phenomenological parameter. This has been
            % introduced by Debernardi to eliminate strange oscillations
            % occurring when including non-Markovian effects.
            % Alternative, more rigorous models for the inclusion of
            % carrier-carrier scattering in the lineshape broadening can be
            % found in 1989Asada_JQE, or with the keyboard "Landsberg
            % broadening", see 1977Martin_SSC and 1966Landsberg_PSS
            if strcmp(iline,'nMark')  % gammak phenomenological correction
                ExpGamN=(EFc+EFv-mesh.Eg)*qel/kBT;
                gammakN=gammakS.*(1+exp(ExpGamN/mode.Expgammak));
                ganm=ganm0.*(1+exp(ExpGamN/mode.Expgammak));
                if indE==1
%                'verG', keyboard
                end
                %
                ExpGam=(DeltaE)*hbar/kBT;
            end
            ganm=ganm0;
            
            %==============================================================
            % lineshape
            % ganm
            Lor=1./line_shape(DeltaE,gammakS,ganm,iline);
            %
            % All quantities with "0" ending refer to zero density of
            % carriers; this is useful to subtract from DeltaN its value
            % with no carriers.
            DeltaE0=ECV_tot0-omega; % tutto in 1/s
            Lor0=1./line_shape(DeltaE0,gammakS,ganm,iline);
            %
            %==============================================================
            % Computation of integrand functions
            %==============================================================
            % All quantities named "In" are meant to be "Integrand"
            % functions in spherical k domain (so, differential k*dk). The
            % integrand can be found in 1988AhnD_JQE, eq. (8). The
            % contributions are: the square absolute value of the dipole
            % momentum, the Fermi factor, and the lineshape (Lorentzian,
            % non-Markovian or Landsberg).
            if iland==1 % "pure" Landsberg model: cutting when gain is less than 0
                % The multiplication times the "condition" function "cond"
                % which is different from 0 only for 0 <= x <= 1, is
                % corresponding to limiting the integrating bounds from 0 
                % to EFv+EFc, in a similar way to 1977Martin_SSC, eq. (3)
                In0=-real(Lor0).*M20.*(cond); % integrand computed for 0 carriers for computing DeltaN/N
                In_att=M2.*Fer.*Lor.*(cond); % integrand computing with carriers
            else % "modified" Landsberg model, and other lineshapes
                In0=-real(Lor0).*M20; % integrand computed for 0 carriers for computing DeltaN/N
                In_att=M2.*Fer.*Lor; % integrand computing with carriers
            end
            % Integrand for DeltaN/N and gain
            In=(In_att-In0)*kdkf;
            InI=(In_att-In0);                      
            
%            'Integrand Anders', keyboard
            
            % Integrand for spontaneous recombination
            In_es=k0(indE)*c_light/Nb*((M2.*Fersp.*Lor)*kdkf);
            GG=sum(In(mode.ntrans));
            
            %==============================================================
            % Integrals: final products
            %==============================================================
            % Gain, 1/cm
            G(idh,ide,indE)=k0(indE)*imag(GG)/100;
            
            % Relative refractive index change Delta_n/n, dimensionless
            Dep(idh,ide,indE)=-real(GG);
            
            % Spontaneous emission rate, 1/s
            Es(idh,ide,indE)=imag(sum(In_es(mode.ntrans)));
            
            %==============================================================
            % Debug plots
            %==============================================================
            iplot=0;
            if iplot==1
                figure, plot(kgrid,dip), pausak
                figure, plot(kgrid,lor,kgrid,imag(lor)), pausak
                figure, plot(kgrid,fer), pausak
                figure, plot(kgrid,dele), title('Del E'), pausak
                figure, plot(kgrid,imag(gintegrand)), title('Integrand')
                pausak
            end
            %==============================================================
        end  % energy (lambda) loop
%            'qui vedo', keyboard
            if mode.ieh_equal==0            
             EFsplit(idh,ide)=EFvv(idh)+EFcv(ide)+hbar.*DelGap(1,1)./qel-DeltaETemp;
            else
             EFsplit(idh,ide)=EFvv(idh)+EFcv(idh)+hbar.*DelGap(1,1)./qel-DeltaETemp;
            end
            
%            EFsplit(idh,ide)=EFvv(idh)+EFcv(ide)+hbar.*DelGap(1,1)./qel-DeltaETemp;

      if(mode.LUT==0)
            EFsave(idh,ide)=EFvv(idh)+EFcv(ide)+hbar.*DelGap(1,1)./qel-DeltaETemp;
            ECVtotsave(idh,ide)=ECV_tot(1,1)+DelGap(1,1)*hbar./qel-DeltaETemp;
            
%            'qui vedo', keyboard
        end
    end % hDensity loop
end % eDensity loop
%            'qui vedo', keyboard

