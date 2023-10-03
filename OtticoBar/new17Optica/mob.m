function [mue,muh,me,mh,tc]=momatc(x,T)

x0=.45;
c1=[10000 -22000 8000];
c2=[-720 1160 -255];
ch=[740 -970 370];

m1=[83 67]*1e-3;
m2=[-60 320]*1e-3;
mh=[140 620]*1e-3;

ct=[-13.22 12.7 1];


fi1=find(x<x0);
y1=polyval(c1,x(fi1));
fi2=find(x>=x0);
y2 =polyval(c2,x(fi2));
mue=[y1; y2];
muh =polyval(ch,x);

y1=polyval(m1,x(fi1));
fi2=find(x>=x0);
y2 =polyval(m2,x(fi2));
me=[y1; y2];
mh =polyval(mh,x);

tc =550./polyval(ct,x)*T^(-1.25);

