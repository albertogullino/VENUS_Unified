% Impulso rettangolare. centrato in mu; supporto=T; max=1
function [y]= rect(x,mu,T)
n=length(x);
y=zeros(1,n);
i=find(abs(x-mu)<=T/2);
m=length(i);
y(i)=ones(1,m);
return