%==========================================================================
function [eDensity,hDensity,eDenV,hDenV,We,Wh,Wed,Whd] = ...
    f_charge_WFk(mesh,Ban,EFc,EFv,kBTev)
%==========================================================================
%
ipri=0;
kgrid=Ban.kgrid;
SB1=Ban.SBC;
SB2=Ban.SBV;
XVC=Ban.XVC;
XVV=Ban.XVV;
C0=mesh.C0+mesh.Eg;
V0=-mesh.V0;
%==========================================================================
eDensity = 0;
We=0;
Wed=0;
Ded=0;
for ib=1:mesh.ncb;
    engy = SB1(ib,:); % eV
    fi=find(engy>=C0);
    fpes=ones(size(kgrid));
    fpes(fi)=0;
    eX=exp((engy-EFc)/(kBTev));
    f = fpes./(1.0D0 + eX);
    wi=XVC(:,ib).^2;
    fd = eX.*f.^2;
    eDensityd = 2/(2*pi)^2*(2*pi)*trapz(kgrid,f.*kgrid);
    eDer = 2/(2*pi)^2*(2*pi)*trapz(kgrid,fd.*kgrid);
    eDensity  = eDensity + eDensityd;
    We=We+eDensityd*wi;
    Ded=Ded+eDer;
    Wed=Wed+eDer*wi;
    eDenV(ib) = eDensityd;
end
Wed=Wed/Ded;

if ipri==1
    fprintf('EDENSITY,     (1/cm^2): %e\n',eDensity * 1.0D-4)
    fprintf('EFc               (eV): %e\n',EFc)
end
%==========================================================================
hDensity = 0;
Wh=0;
Whd=0;
Dhd=0;
kdk=kgrid*diff(kgrid(2:3));
kdk(end)=kdk(end)/2;
kdk=kdk/pi;
for ib=1:mesh.nvb;
    engy = SB2(ib,:); % eV
    fi=find(engy>=V0);
    fpes=ones(size(kgrid));
    fpes(fi)=0;
    wi = (squeeze(XVV(:,ib,:)));
    eX=exp((engy-EFv)/(kBTev));
    f = fpes./(1.0D0 + eX);
    fd = eX.*f.^2;
    hDensityd = 2/(2*pi)^2*(2*pi)*trapz(kgrid,f.*kgrid);
    eDer = 2/(2*pi)^2*(2*pi)*trapz(kgrid,fd.*kgrid);
%     wh_charge=wi*(f.*kdk)';
%     wh_charged=wi*(fd.*kdk)';
    wh_charge=wi*(f.*kdk);
    wh_charged=wi*(fd.*kdk);
    Wh=Wh+wh_charge;
    Whd=Whd+wh_charged;
    Dhd = Dhd + eDer;
    hDensity = hDensity + hDensityd;
    hDenV(ib) =hDensityd;
end
Whd=Whd/Dhd;
if ipri==1
    fprintf('HDENSITY,     (1/cm^2): %e\n',hDensity * 1.0D-4)
    fprintf('EFv               (eV): %e\n',EFv)
end
%==========================================================================