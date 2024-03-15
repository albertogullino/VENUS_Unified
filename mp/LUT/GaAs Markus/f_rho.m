%==========================================================================
function [out] = f_rho(kgrid,SB,EF,kBTev,Density_target,Ebar,mode)
%==========================================================================
Density = 0; nb = size(SB,1);
for ib = 1:nb;
    engy = SB(ib,:); % eV
    fi=find(engy>=Ebar);
    %'engy',keyboard
    fpes=ones(size(kgrid));
    if mode.flagLimitMaxDensity==1
        fpes(fi)=0;
    end
    eX=exp((engy-EF)/(kBTev));
    f = fpes./(1.0D0 + eX);
    Density = Density + 2/(2*pi)^2*(2*pi)*trapz(kgrid,f.*kgrid);
    %'controllo ', keyboard
end

% fprintf('EDENSITY,     (1/cm^2): %e\n',Density * 1.0D-4)
out = abs((Density_target - (Density*1.0D-4))/Density_target);
%==========================================================================