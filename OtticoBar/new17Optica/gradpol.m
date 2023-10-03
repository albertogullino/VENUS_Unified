function [x,f,xn]=gradpol(D,x1,x2,cox,n,cxp,imat);

if length(nargin)<7
 imat=0;
end 
 
np=n+1;
xp=linspace(0,D,np)';
x=diff(xp(1:np));
xc=linspace(x1,x2,np)';
xp=xc+diff(xc(1:2))/2;
xn=xp(1:n);
if length(cox)==1
 if imat==0
  f=nAlGaAs(cox,xn);
 else 
  f=real(nIngaas(cox,xn));
 end 
else
 f=polyval(cox,xn);
end
if length(cxp)~=0
 p=polyval(cxp,xn);
 f=f+i*p;
end
%p, f, pausak

%xp=linspace(0,X,100)';
%yp=linspace(f1,f2,100)';
%figure, plot(xp,yp,x,f,'*')

