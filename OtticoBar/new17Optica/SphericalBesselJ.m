function w=SphericalBesselJ(nu,z)

w=zeros(length(z),length(nu));
indc=find(nu==0);
indr=find(z==0);
w(indr,indc)=1;
indznz=find(z~=0);
zarg=angle(z(indznz));
zargmat=zarg*ones(1,length(nu));
w(indznz,:)=(sqrt(pi./(2*abs(z(indznz))))*ones(1,length(nu))).*besselj(nu+0.5,z(indznz)).*exp(-j*zargmat/2);

