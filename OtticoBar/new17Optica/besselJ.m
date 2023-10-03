
function J = besselJ(nu,z)

NU=ones(length(z),1)*nu;
Z=z*ones(1,length(nu));
J=besselj(NU,Z);

return
