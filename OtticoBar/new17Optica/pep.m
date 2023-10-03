A1=[1 2; pi 4];
A1=A1/sum(sum(A1));

I=eye(2);

A=I+A1;

D=10;
d=.001;
Et=expm(D*A+d*A1);
E0=expm(D*A);
Ed=expm(d*A);

L=d*A1;

Eda=I+L;
Eprod=E0*Ed;
Eproda=E0*Eda;
Eproda1=E0+E0*L;

Et-Eprod
Et-Eproda
Eproda-Eproda1
