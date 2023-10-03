
function [S11,S21] = f_ComputeSlab(n2,n1,n3,L,vlambda,theta,Polarization)

vk0=2*pi./vlambda;

nin_sintheta=n1.*sin(theta);
kxinc=vk0.*nin_sintheta;

Clight=299.792458; % speed of light (um/ps)
mu0=4*pi/10; % vacuum magnetic permeability (H/um)
Z0=mu0*Clight; % free-space impedance (ohm)

vkz1=f_psqrt((vk0.*n1).^2-kxinc.^2);
vkz2=f_psqrt((vk0.*n2).^2-kxinc.^2);
vkz3=f_psqrt((vk0.*n3).^2-kxinc.^2);

if(strcmp(Polarization,'TE'))
    Zinf_1=Z0*vk0./vkz1;
    Zinf_2=Z0*vk0./vkz2;
    Zinf_3=Z0*vk0./vkz3;
    S11p=(vkz1-vkz2)./(vkz1+vkz2);
    S11s=(vkz2-vkz3)./(vkz2+vkz3);
elseif(strcmp(Polarization,'TM'))
    Zinf_1=Z0*vkz1./(vk0*n1.^2);
    Zinf_2=Z0*vkz2./(vk0*n2.^2);
    Zinf_3=Z0*vkz3./(vk0*n3.^2);
    S11p=(vkz2/n2.^2-vkz1/n1.^2)./(vkz2/n2.^2+vkz1/n1.^2);
    S11s=(vkz3/n3^2-vkz2/n2.^2)./(vkz3/n3.^2+vkz2/n2.^2);
else
    error('Sbagliata selezione della polarizzazione :-(')
end

S21p=2.*sqrt(Zinf_1).*sqrt(Zinf_2)./(Zinf_2+Zinf_1);
S12p=S21p; S22p=-S11p;

S21s=2.*sqrt(Zinf_2).*sqrt(Zinf_3)./(Zinf_3+Zinf_2);
S12s=S21s; S22s=-S11s;

S11=S11p+(S11s.*S21p.*S12p.*exp(-j.*2.*vkz2.*L))./(1-S22p.*S11s.*exp(-j.*2.*vkz2.*L));
S21=S21p.*S21s.*exp(-j.*vkz2.*L)./(1-S22p.*S11s.*exp(-j.*2.*vkz2.*L));

return