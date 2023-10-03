function [rp]=peAu(Amr)
F=abs(Amr);
l=length(F);
pu=3:l/2-1;
pud=2:l/2;
ff=F(pu);
df=diff(F(pud));
Fz=df(1:end-1).*df(2:end);
fized=find(Fz<0);
if length(fized)==1
 fized=[fized length(pu)];
end
PU=pu(fized(1:2));
%figure, plot(F), hold on, plot(PU,F(PU),'ro'), pausak, close
rp=PU(2)-PU(1);
rp=PU(1);