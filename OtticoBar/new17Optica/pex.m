rr=3;
n=1;
k=.1;
L=1;
lam=1;
k0=2*pi/lam;
ber=k0*rr;
Del=(n/rr)^2-1;

kt=k*ber;

t=sqrt(1-k^2);


M0=-j*ber*L*t*( [1 0; 0 -1] + 1/(2*t^2)* Del * [1 1; -1 -1] ); 

autov=eig(M0)
analitico=-j*k0*sqrt(1-(kt/k0)^2)*L;
T=expm(M0);
ain=[1 0]';

out=T*ain
out_ana=exp(analitico)

