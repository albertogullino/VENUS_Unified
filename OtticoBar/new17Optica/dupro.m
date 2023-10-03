
mass0=9.1e-31;
mi=pi*4e-7;
eps0=8.8541e-12;
c=1/sqrt(mi*eps0);
Z0=120*pi;
q=1.6e-19;


calfa=1e10*Z0*q^3*lambda0^3/(mass0*2*pi*c)^2/(4*pi);
cenne=1e6*q^2*lambda0^2/(2*eps0*mass0*(2*pi*c)^2);


fiA=find(Dop>0);
palfae=calfa*abs(Dop(fiA)')./(mue(fiA).*me(fiA).^2.*real(nv(fiA,1).'))/(1+fre)^3;
pennee=cenne*abs(Dop(fiA)')./(me(fiA).*real(nv(fiA,1).'))/(1+fre)^3;
figure, plot(palfae)
figure, plot(pennee,'r')

fiA=find(Dop<0);
palfah=calfa*abs(Dop(fiA)')./(muh(fiA).*mh(fiA).^2.*real(nv(fiA,1).'))/(1+fre)^3;
penneh=cenne*abs(Dop(fiA)')./(mh(fiA).*real(nv(fiA,1).'))/(1+fre)^3;
figure, plot(palfah)
figure, plot(penneh,'r')

