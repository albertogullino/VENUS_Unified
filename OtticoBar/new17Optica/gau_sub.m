iplo=1;

w2=w^2;
NP=50;
NPkx=35;
NPky=35;
%NPkx=200;
%NPky=250;
x=linspace(-a*w,a*w,NP);
y=x;

E=exp(-x.^2/(2*w2))'*exp(-y.^2/(2*w2));

if iplo==1
figure, plot(x,E), pausak
end

ks=1.5;
kx=linspace(-a/w*ks,a/w*ks,NPkx);
fi=find(kx==0);
kx(fi)=1e-10;
ky=linspace(-a/w*ks,a/w*ks,NPky);
fi=find(ky==0);
ky(fi)=1e-10;
kym=repmat(ky,NPkx,1);
kxm=repmat(kx',1,NPky);
kt2=kxm.^2+kym.^2;
kt=sqrt(kt2);
Ex=exp(-kt2*w2/2);
Vie=-w2*kym./kt.*exp(-kt2*w2/2);
Vim=w2*kxm./kt.*exp(-kt2*w2/2);

%Te1=w2*exp(-kt2*w2/2);
%Tm1=w2*exp(-kt2*w2/2);
if iplo==1
figure, plot(kx,Vie), pausak
figure, plot(kx,Vim), pausak
end
% ricostruzione
'prima di GA', keyboard
for ikx=1:length(kx)
 for iky=1:length(ky)
  Ga=gaHCG2(kx(ikx),ky(iky),pg);
  Vr=Ga*[Vie(ikx,iky); Vim(ikx,iky)];
  Vre(ikx,iky)=Vr(1);
  Vrm(ikx,iky)=Vr(2);
  Ge(ikx,iky)=Ga(1,1);
  Gm(ikx,iky)=Ga(2,2);
  Gem(ikx,iky)=Ga(1,2);
  Gme(ikx,iky)=Ga(2,1);
 end
end


dkx=diff(kx(1:2));
dky=diff(ky(1:2));
for ix=1:length(x)
 for iy=1:length(y)
 ep=exp(-j*(kxm*x(ix)+kym*y(iy)));
 ee=kxm./kt.*ep;
 em=kym./kt.*ep;
 eex=-kym./kt.*ep;
 emx=kxm./kt.*ep; 
 Eiy(ix,iy)=sum(sum(Vie.*ee+Vim.*em))*dkx*dky/(2*pi);
 Eix(ix,iy)=sum(sum(Vie.*eex+Vim.*emx))*dkx*dky/(2*pi);
 Ery(ix,iy)=sum(sum(Vre.*ee+Vrm.*em))*dkx*dky/(2*pi);
 Erx(ix,iy)=sum(sum(Vre.*eex+Vrm.*emx))*dkx*dky/(2*pi); 
 end
end 


for ix=1:length(x)
 for iy=1:length(y)
 ep=exp(-j*(kxm*x(ix)+kym*y(iy)));
 ee=kxm./kt.*ep;
 em=kym./kt.*ep;
 eex=-kym./kt.*ep;
 emx=kxm./kt.*ep; 
 Eiye(ix,iy)=sum(sum(Vie.*ee))*dkx*dky/(2*pi);
 Eixe(ix,iy)=sum(sum(Vie.*eex))*dkx*dky/(2*pi);
 Eiym(ix,iy)=sum(sum(Vim.*em))*dkx*dky/(2*pi);
 Eixm(ix,iy)=sum(sum(Vim.*emx))*dkx*dky/(2*pi); 
 end
end 


if iplo==1

figure, plot(x,Eiy), hold on, plot(x,E,'.'), pausak
%figure, semilogy(x,Er), hold on, semilogy(x,E,'.'), pausak
figure, plot(x,Eiy-E), pausak
figure, plot(x,Eix), pausak
figure, plot(x,Ery), hold on, plot(x,Eiy,'.'), pausak
figure, plot(x,Erx), hold on, plot(x,Eiy,'.'), pausak
end
Iinc=sum(sum(abs(Eiy).^2+abs(Eix).^2));
Iref=sum(sum(abs(Erx.^2)+abs(Ery.^2)));
R_eq=Iref/Iinc;

return

lambda=1.616;
k0=2*pi/lambda;
tetai=asin(kt/k0);
phii=atan2(kym,kxm);