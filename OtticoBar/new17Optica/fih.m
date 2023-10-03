load saca

figure, plot(zet,real(Ez),'r.',zet,imag(Ez),'c.')

ia=20;
ib=21;
fa=Ez(ia);
fap=diff(Ez(ia+[-1 0]))/diff(zet(ia+[-1 0]));

fb=Ez(ib);
fbp=diff(Ez(ib+[0 1]))/diff(zet(ib+[0 1]));
d=diff(zet([ia ib]));


zeh=linspace(zet(ia),zet(ib),50);

nm=1000;
k0=2*pi/lambda;
nme=linspace(1,5,nm);
for k=1:nm
 x0=k0*nme(k)*d;
 c0=cos(x0);
 s0=sin(x0);
 M(1,1)=-fa*s0;
 M(1,2)=fap*c0*d/x0;
 M(2,1)=-fa*c0*x0/d;
 M(2,2)=-fap*s0;
 N(1,1)=fb-fa*c0-fap*d*s0/x0;
 N(2,1)=fbp+fa*x0/d-fap*c0;
 dev=inv(M)*N;
 dva(k)=dev(1);
 dvb(k)=dev(2);
end

figure, plot(nme, abs(dva),nme, abs(dvb), nme, k0*nme*d)