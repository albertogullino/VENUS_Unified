Ex=P.E2xp{k};
Ef=P.Ef{k};
Ex0=mean(Ex,2);
Ef0=mean(Ef,2);

I=abs(Ex0.^2);
% si2=20^2
% I=exp(-2*xro'.^2/si2);
Es=xdx2*I;
En=xdx*I;
sig2=Es./En;
%'m2', keyboard


Esf=xdx2f*abs(Ef0.^2);
Enf=xdxf*abs(Ef0.^2);
sig2f=Esf./Enf;

dsi0=2*sqrt(sig2*2);
dsi0f=2*sqrt(sig2f*2);

M2=pi/(4*lambda*z)*dsi0.*dsi0f;

%k
%'M2 ver', keyboard