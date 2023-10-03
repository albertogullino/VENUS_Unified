function f=filtA(y,so)

ay=abs(y);
fiva=find(ay>1e-5);
ay=ay(fiva);
me=median(ay);
ma=max(ay);
f=y;
 fiva=find(abs(f)<1e-5);
 f(fiva)=1e-10;
if ma/me>so
 fin=find(abs(y)>so*me);
 f(fin)=1e-10;
% figure, plot(abs(f)), hold on, plot(abs(y),'r'),  'filt', pausak
end
