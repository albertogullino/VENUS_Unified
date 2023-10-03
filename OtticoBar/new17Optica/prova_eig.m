load pro3


M=expm(P);

  [V,autov]=eig(P);
  daut=diag(autov);
  ver=P*V-V*autov;

Mdef = V * diag(exp(daut)) / V;

map(log10(ver))

map(log10(M-Mdef))