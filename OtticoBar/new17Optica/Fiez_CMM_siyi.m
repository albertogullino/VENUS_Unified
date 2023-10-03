function [ga,la]=fiez_CMM(lam0,L_i,n_i,iat,icav,ifp,dlam);


PO=[909    83   672   711];
if ifp==-10
figure, plot(L_i,'o'), hold on, plot(iat,L_i(iat),'rp')
a=axis;
a(3)=0;
a(4)=max(L_i)+10;
axis(a)
pausak
end

lav=lam0+linspace(-1,1,3)*dlam*lam0;
%lav=lam0-.01;

for k=1:length(lav)
lai=lav(k);
[ei,g0]=Flam_gamL(L_i,n_i,iat,icav,lai);
gov(k)=g0;
eiv(k)=ei;
end

%'qui y', keyboard

g00=g0;
alf=imag(eiv);
coa=polyfit(lav,alf,1);
LAM=roots(coa);

gv=real(eiv);
cog=polyfit(lav,gv,2);
GTH=polyval(cog,LAM);
G0=GTH;
cog0=polyfit(lav,gov,2);
Gga=polyval(cog0,LAM);

if ifp==-10
h=figure, 
set(h,'pos',PO)
subplot(211)
plot(lav,eiv,LAM,GTH,'wo'), title('real part'), grid
subplot(212)
plot(lav,imag(eiv),LAM,0,'wo'), title('imag part'), grid
pausak

%'inizio iterazione', pausak
end

G00=G0;
Gprec=G0;
G0=0;
iter=1;
Gve(iter)=Gprec;
Lve(iter)=LAM;
while abs(G0/Gprec-1)>1e-3 & iter<2

lav=LAM*(1+linspace(-2,1,2)*10^-(3+iter));
lav=LAM;
g0=Gprec*1.1;
%g0=Gprec*.9;

for k=1:length(lav)
lai=lav(k);
%' qui prima U', keyboard
[ei]=FlamU(L_i,n_i,iat,lai,g0);
eiv1(k)=ei;
end
eiv=eiv1+g0;

G0=real(eiv);
iter=iter+1;
Gp=Gprec;
Gprec=G0;
Gve(iter)=Gprec;
Lve(iter)=LAM;

G0=Gp;
%Gprec=G0;
abs(G0/Gprec-1);
%pausak
end


Gx=[Gga g0];
cg=polyfit(Gx,Gve,1);
x=linspace(min(Gx),max(Gx),10);
cg0=cg;
cg0(1)=cg0(1)-1;
gg0=roots(cg0);
Gstim=gg0;
ga=Gstim;
la=LAM;


if ifp==-10
figure, plot(x,polyval(cg,x),gg0,gg0,'ro')
pausak
end