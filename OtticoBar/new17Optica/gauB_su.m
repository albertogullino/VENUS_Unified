ima=1;
NP=50;
NPkx=10000;
NPkx=50;
%NPkx=200;
%NPky=250;
r=linspace(0,a*w,NP);
x=r;

E=exp(-x.^2/(2*w2));
if ima==1
figure, plot(x,E), pausak
end

ks=.7;
kt=linspace(0,a/w*ks,NPkx+1);
kt=kt(2:end);
kx=kt;
Ex=exp(-kt.^2*w2/2);
Ve=w2*Ex;
Vm=-Ve;
 dk=diff(kt(1:2));

if ima==1
figure, plot(kx,Ve), pausak
end
clear Ex Er1x Er2x
for ix=1:length(r)
 ep=besselj(0,kx*r(ix)).*kt/2;
 Ex(ix)=sum((Ve-Vm).*ep)*dk;
end 

if ima==1
figure, plot(x,Ex,x,E,'.'),  pausak
figure, plot(x,abs(Ex-E),'.'),  pausak
end

VieX=Ve;
VimX=Vm;
alf=linspace(0,pi,31);
%alf=linspace(0,2*pi,51);
alf=alf(2:end);
kc0=exp(-j*alf);
kx0=real(kc0);
ky0=imag(kc0);

no=alf(end)/2;
for ikx=1:length(kt)
 kx=kt(ikx)*kx0;
 ky=kt(ikx)*ky0;
  Vi=[VieX(ikx); VimX(ikx)];
 for kal=1:length(alf) 
  [Ga,Tr]=gaHCG1(kx(kal),ky(kal));
  Vr=Ga*Vi;
  Vt=Tr*Vi;
%  'fine X'
  Vek(kal)=Vr(1);
  Vmk(kal)=Vr(2);
  Vtek(kal)=Vt(1);
  Vtmk(kal)=Vt(2);  
  kxv(kal,ikx)=kx(kal);
  kyv(kal,ikx)=ky(kal);
  Ge(kal,ikx)=Ga(1,1);
  Gm(kal,ikx)=Ga(2,2);
  Gem(kal,ikx)=Ga(1,2);
 end 
 
  fIe1=sum(Vek.*sin(alf).^2)*diff(alf(1:2))/no;
  fIm1=sum(Vmk.*cos(alf).^2)*diff(alf(1:2))/no;
  fIe2=sum(Vek.*cos(alf).^2)*diff(alf(1:2))/no;
  fIm2=sum(Vmk.*sin(alf).^2)*diff(alf(1:2))/no;  

  fte1=sum(Vtek.*sin(alf).^2)*diff(alf(1:2))/no;
  ftm1=sum(Vtmk.*cos(alf).^2)*diff(alf(1:2))/no;
  fte2=sum(Vtek.*cos(alf).^2)*diff(alf(1:2))/no;
  ftm2=sum(Vtmk.*sin(alf).^2)*diff(alf(1:2))/no;  


  Ver1(ikx)=fIe1;
  Vmr1(ikx)=fIm1;
  Vter1(ikx)=fte1;
  Vtmr1(ikx)=ftm1;
  
  fIe1a=sum(Vek.*sin(alf).*cos(alf))*diff(alf(1:2))/no;
  fIm1a=sum(Vmk.*cos(alf).*sin(alf))*diff(alf(1:2))/no;
  Ver1a(ikx)=fIe1a;
  Vmr1a(ikx)=fIm1a;
  
  fte1a=sum(Vtek.*sin(alf).*cos(alf))*diff(alf(1:2))/no;
  ftm1a=sum(Vtmk.*cos(alf).*sin(alf))*diff(alf(1:2))/no;
  Vter1a(ikx)=fte1a;
  Vtmr1a(ikx)=ftm1a;  
%' fine alfa ', keyboard    
  Ver2(ikx)=fIe2;
  Vmr2(ikx)=fIm2;  
  Vter2(ikx)=fte2;
  Vtmr2(ikx)=ftm2;    
end

% Ver1=Ver1s;
% Vmr1=Vmr1s ;

Ver1u=Ver1-Ver1a;
Vmr1u=Vmr1-Vmr1a;
Ver2u=Ver2+Ver1a;
Vmr2u=Vmr2+Vmr1a;

Vter1u=Vter1-Vter1a;
Vtmr1u=Vtmr1-Vtmr1a;
Vter2u=Vter2+Vter1a;
Vtmr2u=Vtmr2+Vtmr1a;
% ricostruzione

Nf=61;
phi=linspace(0,2*pi,Nf)';
ff=cos(2*phi);
fg=sin(2*phi);
o=ones(size(phi));
clear Ex1 Ex2 Ey1 Ey2
clear Erx1 Erx2 Ery1 Ery2
for ix=1:length(r)
 ep0=besselj(0,kt*r(ix)).*kt/2*dk;
 ep2=besselj(2,kt*r(ix)).*kt/2*dk;
 epex1=ff*ep2+o*ep0;
 epey1=fg*ep2;
 epmx1=ff*ep2-o*ep0;
 epmy1=fg*ep2; 
 epex2=fg*ep2;
 epey2=-ff*ep2+o*ep0;
 epmx2=fg*ep2;
 epmy2=-ff*ep2-o*ep0;  
 Ex1(ix,:)=(Ve*epex1.'+Vm*epmx1.');
 Ey1(ix,:)=(Ve*epey1.'+Vm*epmy1.');
 Ex2(ix,:)=(Ve*epex2.'+Vm*epmx2.');
 Ey2(ix,:)=(Ve*epey2.'+Vm*epmy2.'); 

 Erx1(ix,:)=(Ver1u*epex1.'+Vmr1u*epmx1.');
 Ery1(ix,:)=(Ver1u*epey1.'+Vmr1u*epmy1.');
 Erx2(ix,:)=(Ver2u*epex2.'+Vmr2u*epmx2.');
 Ery2(ix,:)=(Ver2u*epey2.'+Vmr2u*epmy2.'); 

 Etx1(ix,:)=(Vter1u*epex1.'+Vtmr1u*epmx1.');
 Ety1(ix,:)=(Vter1u*epey1.'+Vtmr1u*epmy1.');
 Etx2(ix,:)=(Vter2u*epex2.'+Vtmr2u*epmx2.');
 Ety2(ix,:)=(Vter2u*epey2.'+Vtmr2u*epmy2.'); 

end 

%ima=0;
if ima==1
X=(cos(phi)*r)';
Y=(sin(phi)*r)';
map(Ex1,X,Y), title(' X component Pol1'), pausak
map(Ey1,X,Y), title(' Y component Pol1'), pausak
map(Erx1,X,Y), title(' ref X component Pol1'), pausak
map(Ery1,X,Y), title(' ref Y component Pol1'), pausak
map(Etx1,X,Y), title(' tras X component Pol1'), pausak
map(Ety1,X,Y), title(' tras Y component Pol1'), pausak

map(Ex2,X,Y), title(' X component Pol2'), pausak
map(Ey2,X,Y), title(' Y component Pol2'), pausak
map(Erx2,X,Y), title(' ref X component Pol2'), pausak
map(Ery2,X,Y), title(' ref Y component Pol2'), pausak
map(Etx2,X,Y), title(' tras X component Pol2'), pausak
map(Ety2,X,Y), title(' tras Y component Pol2'), pausak
end

Iinc1=sum(sum(abs(Ex1.^2)+abs(Ey1.^2)));
Iinc2=sum(sum(abs(Ex2.^2)+abs(Ey2.^2)));
Ir1=sum(sum(abs(Erx1.^2)+abs(Ery1.^2)));
Ir2=sum(sum(abs(Erx2.^2)+abs(Ery2.^2)));
It1=sum(sum(abs(Etx1.^2)+abs(Ety1.^2)));
It2=sum(sum(abs(Etx2.^2)+abs(Ety2.^2)));
R1_eq=Ir1/Iinc1
R2_eq=Ir2/Iinc2
T1_eq=It1/Iinc1
T2_eq=It2/Iinc2

