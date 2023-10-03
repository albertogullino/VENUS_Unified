load zer

sti0=linspace(-1,1,30)/10;
sti=sti0*imag(z0);
zz=z0+j*sti;

for kvv=1:length(sti)
 z=zz(kvv);
 vz(kvv)=f_mulut(z);
end

figure, plot(sti0,real(vz), sti0,imag(vz))