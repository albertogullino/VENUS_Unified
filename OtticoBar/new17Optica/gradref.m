function [x,f,xd]=gradref(D,x1,x2,n);
np=n+1;
xp=linspace(0,D,np)';
xp1=linspace(0,D,n)';
x=diff(xp(1:np));
f=x1+(x2-x1)/D*xp1;
xd=-10*ones(size(x));
