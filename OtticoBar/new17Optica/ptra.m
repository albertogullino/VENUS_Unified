clear all
close all

z=linspace(0,5,1001);
L=3;
rr=5;
n=3.4;

lam=3.15;
k0=2*pi/lam;
Ga=(rr-n)/(rr+n);
Tr=sqrt(1-Ga^2);
V=[1 -Ga; -Ga 1]/Tr;

fim=find(z<=L);
fid=find(z>L);

ap=exp(-j*k0*rr*z(fim));
ar=Ga*exp(-j*k0*rr*(2*L-z(fim)));
Vs=(ap+ar);
Is=(ap-ar);


%figure, plot(z(fim),abs(Vs),z(fim),abs(Is))

at=V*[ap(end); ar(end)];

Vt=[Vs sqrt(rr/n)*sum(at)*exp(-j*k0*n*(z(fid)-L))];
It=[Is sqrt(n/rr)*diff(-at)*exp(-j*k0*n*(z(fid)-L))];


It1=[Is   sqrt(n/rr)*Tr*ap(end)*exp(-j*k0*n*(z(fid)-L))];
Vt1=[Vs  sqrt(rr/n)*Tr*ap(end)*exp(-j*k0*n*(z(fid)-L))];

figure, plot(z,abs(Vt),'.',z,abs(It),'.'), 
title(' V-I Trasmissione')
pausak
figure, plot(z,abs(Vt1),'.',z,abs(It1),'.'), 
title(' V-I Scattering')
pausak

figure, plot(z,abs(Vt),'.',z,abs(Vt1),'o'), 
title(' V Trasmissione-Scattering')
pausak
figure, plot(z,abs(It),'.',z,abs(It1),'o'), 
title(' I Trasmissione-Scattering')
pausak

figure, plot(z,Vt.*conj(It),'. '), pausak
figure, plot(z,Vt1.*conj(It1),'.'), pausak