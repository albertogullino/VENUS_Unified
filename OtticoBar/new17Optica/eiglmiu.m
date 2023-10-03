function [Ksi,lambdas,Fi,uLong,uLong0,Fa]=eiglmiu(Lvp0,uFunc,uF0,relPerm,hn,din,iraff);

d=din;

%'entro miu', pausak

if ~exist('iraff')
 iraff=0;
end


[Ksi,lambdas,Fi,uLong,uLong0,Fa]=eiglmiu_zero(Lvp0,uFunc,uF0,relPerm,hn,d);
Lvp0=lambdas;

if iraff==0
  return
end

fiat=find(uF0>0);

ra=sqrt(relPerm(fiat(ceil(end/2))));
rqw=ra;
lambda=lambdas*1e6;
k0=2e4*pi/lambda;
gpla=rqw*k0*imag(Ksi)/uLong;



 
 z=cumsum(ones(size(Fi)))*hn;
 
 %figure, plot(z,Fi,z(fiat),Fi(fiat),'r.'), pausak

g0=gpla;
gplaf=gpla;
ic=0;
stg=[0 1]*gpla/20;
for ic=1:length(stg)
gplaf=g0-stg(ic);
nI=gplaf/(2*k0);
nref=sqrt(relPerm)-j*uFunc*nI;
relC=nref.^2;
[Ksi1,lambdas1,Fi1,uLong1,uLong01,Fa1]=eiglmiu_zero(Lvp0,uFunc,uF0,relC,hn,d);
 Lvp0=lambdas1;
 err=k0*rqw*imag(Ksi1)/uLong1;
 gpla1(ic)=err;
 %kpla1(ic)=imag(Ksi1);
 %lpla1(ic)=lambdas1;
 lambda=lambdas1*1e6;
 k0=2e4*pi/lambda; 
end



co=polyfit(stg,gpla1,1);
zeg=roots(co);

%figure, plot(stg,gpla1,zeg,polyval(co,zeg),'wo'), pausak

gver=g0-zeg;
gplaf=gver;
nI=gplaf/(2*k0);
nref=sqrt(relPerm)-j*uFunc*nI;
relC=nref.^2;
%relC=relPerm-j*uF0*nI*2*rqw;
[Ksi,lambdas,Fi1,uLong,uLong0,Fa]=eiglmiu_zero(Lvp0,uFunc,uF0,relC,hn,d);

 %figure, plot(z,Fi,z,Fi1,'r.'), pausak
 %g=k0*rqw*imag(Ksi1)/uLong1;
 lambda=lambdas*1e6;
 k0=2e4*pi/lambda;  

Ksi=j*(gplaf*uLong)/(k0*rqw);
%' fine', keyboard

return
%stp=[0 stg];
%figure, plot(stp,gpla1), pausak

fiat=find(uF0>0);

 figure, plot(z,relPerm,z,uF0*relPerm(fiat(10)),'r.'), 
 
 