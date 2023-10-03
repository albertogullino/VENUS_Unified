function rocf=sha_gen;

R=.8
%disp(' R'), keyboard
Delta=-.5;
Delta=0;

ce0=1+j*1;

npfia=21;
np00=11;

%fi0=linspace(0,pi/2,npfia)';
%fi=[fi0; pi-flipud(fi0); pi+fi0; 2*pi-flipud(fi0)];

fi=linspace(0,pi/2,npfia)';
roc=ce0+R*exp(j*fi);


fir1=linspace(0,1,np00+2);
fir=fir1(2:end-1)';
rocd=real(roc(1))+j*imag(roc(1))*fir;
rocu=[rocd; roc];
rocd=j*imag(roc(end))+real(roc(end))*fir;
rocu=[rocu; rocd];
rocf=[rocu; -conj(flipud(rocu))];
rocf=[rocf; -(rocf)];
rocf=[rocf; rocf(1)];

figure, plot(rocf), axis equal


