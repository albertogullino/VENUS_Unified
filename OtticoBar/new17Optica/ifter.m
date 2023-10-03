load safi

fiat=find(uF0>0);
ra=sqrt(relPerm(fiat(end/2)));
fiav=find(relPerm==ra.^2);

%uF0(fiat(1)-[1 2])=1;
%uF0(fiat(end))=0;
%uFunc=uF0;


[Ksi,lambdas,Fi,uLong,uLong0,Fa]=eiglmiu_zero(Lvp0*1e-6,uFunc,uF0,relPerm,hn,d);
lambda=lambdas*1e6;
 gpla=2e4*pi*rqw/lambda*imag(Ksi)/uLong
k0=2e4*pi/lambda;


 
 z=cumsum(ones(size(Fi)))*hn;
 
 figure, plot(z,Fi,z(fiat),Fi(fiat),'r.'), pausak

err=100;
g0=gpla;
gplaf=gpla;
ic=0;
stg=[0 1]*gpla/20;
for ic=1:length(stg)
gplaf=g0-stg(ic);
nI=gplaf/(2*k0);
nref=sqrt(relPerm)-j*uF0*nI;
relC=nref.^2;
[Ksi1,lambdas1,Fi1,uLong1,uLong01,Fa1]=eiglmiu_zero(Lvp0*1e-6,uFunc,uF0,relC,hn,d);
 err=2e4*pi*rqw/lambda1*imag(Ksi1)/uLong1;
 gpla1(ic)=err;
 %if ic>1
 %if err*gpla1(ic-1)<0
 % break
 %end
 %end
lambda=lambdas1*1e6;
k0=2e4*pi/lambda; 
end



co=polyfit(stg,gpla1,1);
zeg=roots(co);

figure, plot(stg,gpla1,zeg,polyval(co,zeg),'wo'), pausak

gver=g0-zeg;
gplaf=gver;
nI=gplaf/(2*k0);
nref=sqrt(relPerm)-j*uF0*nI;
relC=nref.^2;
%relC=relPerm-j*uF0*nI*2*rqw;
[Ksi1,lambdas1,Fi1,uLong1,uLong01,Fa1]=eiglmiu_zero(Lvp0*1e-6,uFunc,uF0,relC,hn,d);
 err=2e4*pi*rqw/lambda1*imag(Ksi1)/uLong1

 figure, plot(z,Fi,z,Fi1,'r.'), pausak

gver
return
%stp=[0 stg];
%figure, plot(stp,gpla1), pausak

fiat=find(uF0>0);
fiav=find(relPerm==ra.^2);

 figure, plot(z,relPerm,z,uF0*relPerm(fiat(10)),'r.'), 
 
 