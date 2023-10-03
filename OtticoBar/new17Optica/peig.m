n=3;

n2=4;

la=1;
k0=2*pi/la;

be=j*k0*n;
del=(n2^2-n^2)/(2*n^2);
M=[-(1+del) -del; del (1+del)+del];

Mv=be*M;

[V,E]=eig(M);