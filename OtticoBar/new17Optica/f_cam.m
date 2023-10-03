function [f1o,ztot,nzo]=f_cam(z)

%'in fmu', keyboard
global fsto L_i n_i rr rfd rfu  iLP ifp iff Lam0 ifun iorta ibast


lambdai=Lam0+real(z);
nim=imag(z);

freq=0;
fiqw=find(fsto==-1);
iff0=1;
%lav=lambda;
ishow=0;

if iorta==0
 [Kte,Ktm,f1o,zio,nzo]=gaz_mtu(fiqw,nim,L_i,n_i,rr,rfd,rfu,lambdai,freq,0,iLP,ifp,iff0);
else 
% [Kte,Ktm,f1o,zio,nzo]=gaz_mtor(fiqw,nim,L_i,n_i,rr,rfd,rfu,lambdai,freq,0,iLP,ifp,iff0);
 [Kte,Ktm,f1o,zio,nzo]=gaz_mtoru(fiqw,nim,L_i,n_i,rr,rfd,rfu,lambdai,freq,0,iLP,ifp,iff0,ibast);
end
%[Kte,Ktm,f1o,zio,nzo]=gaz_mtu(fiqw,nim,L_i,n_i,rr,rfd,rfu,lambdai,freq,0,iLP,ifp,iff0);

ztot=cumsum(zio);

return

if ifp==-10
%'qui f_cam', keyboard
Ez=f1o;
z0=z;
indz=nzo;


fiq=fiqw;
fiq=1:length(Ez);
I=abs(Ez.^2);
%fiq=find(real(indz)==rr & imag(indz)==imag(z0));


if length(fiq)>0
 In=I(1,:)/max(I(1,fiq));
 figure, semilogy(ztot,In*max(indz(fiq)),ztot,real(indz))
% figure, plot(ztot,In*max(indz(fiq)),ztot,real(indz))
else
 In=I(1,:)/max(I(1,:));
 figure, semilogy(ztot,In*max(indz(:)),ztot,real(indz))
end
a=axis;
a(3)=1e-5;
a(4)=4;
axis(a)
pausak
%'in f_cam', keyboard
%keyboard
end