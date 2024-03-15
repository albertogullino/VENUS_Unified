%==========================================================================
function [eDensity,hDensity,eDenV,hDenV] = f_charge_buca(mesh,mode,kgrid,SB1,SB2,EFc,EFv,T,C0,V0)
%==========================================================================
%
qel = 1.6021766208e-019; % electron charge, C
kB = 1.3806488e-23; % Boltzmann constant, J/K
ipri=0;
%==========================================================================
eDensity = 0;
for ib=1:mesh.ncb;
    engy = SB1(ib,:); % eV
    fi=find(engy>=C0);
    fpes=ones(size(kgrid));
    if mode.flagLimitMaxDensity==1
        fpes(fi)=0;
    end
    f = fpes./(1.0D0 + exp((engy-EFc)/(kB*T/qel)));
    eDensityd = 2/(2*pi)^2*(2*pi)*trapz(kgrid,f.*kgrid);
    eDensity = eDensity + eDensityd;
    eDenV(ib) = eDensityd;
end
if ipri==1
    fprintf('EDENSITY,     (1/cm^2): %e\n',eDensity * 1.0D-4)
    fprintf('EFc               (eV): %e\n',EFc)
end
%==========================================================================
hDensity = 0;
for ib=1:mesh.nvb;
    engy = SB2(ib,:); % eV
    fi=find(engy>=V0);
    fpes=ones(size(kgrid));
    if mode.flagLimitMaxDensity==1
        fpes(fi)=0;
    end
    f = fpes./(1.0D0 + exp((engy-EFv)/(kB*T/qel)));
    hDensityd = 2/(2*pi)^2*(2*pi)*trapz(kgrid,f.*kgrid);
    hDensity = hDensity + hDensityd;
    hDenV(ib) =hDensityd;
end
if ipri==1
    fprintf('HDENSITY,     (1/cm^2): %e\n',hDensity * 1.0D-4)
    fprintf('EFv               (eV): %e\n',EFv)
end
%==========================================================================
