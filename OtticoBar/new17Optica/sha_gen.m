function rocf=sha_gen(Psh,fi0,cell);

%cce=[0 j*3+2  4+j*7 7 7*j];
%Rvetx=[.4  1      .5   2  .8];
%Rvety=[.4  1      .5   2  .8];
%Delta=[ 0 -.2     .5  .8 -.05];
%sh_type=[1 2       1   1   2];
%Psh{1}=sh_type;
%Psh{2}=Rvetx;
%Psh{3}=Rvety;
%Psh{4}=Delta;
%Psh{5}=cce;
%'shagen'
%keyboard

sh_type=abs(Psh{1}(cell));
Delta=Psh{4}(cell);
Rx0=Psh{2}(cell);
Ry0=Psh{3}(cell);
cce=Psh{5}(cell);
npfia=length(fi0);
if sh_type==2 | sh_type==3
 Rx=Rx0*(1-Delta);
 Ry=Ry0*(1-Delta);
 cc0=Rx+j*Ry;
 R=sqrt(Rx0*Ry0)*Delta;
% 'r'
% keyboard
   np00=11;
   fi=linspace(0,pi/2,fix(npfia/4))';
   roc=cc0+R*exp(j*fi);
   fir1=linspace(0,1,np00);
   fir=fir1(1:end-1)';
   rocd=real(roc(1))+j*imag(roc(1))*fir;
   rocu=[rocd; roc];
   rocd=j*imag(roc(end))+real(roc(end))*flipud(fir);
%   fi1q=1:length(rocd);
   rocu=[rocu; rocd];
   rocf=[rocu; -conj(flipud(rocu(2:end-1)))];
   rocf=[rocf; -rocf(1:end-1)];
   rocf=rocf+cce;
%   fi4q=3*length(rocd)+fi1q;
elseif sh_type==1 | sh_type==4
 Rx=Rx0;
 Ry=Ry0;
 if sh_type==1
  Delta=0;
  Rx=Rx0;
  Ry=Ry0;
 end
 fi=linspace(0,pi*2,npfia+1)';
% fi1q=find(fi<pi/2);
% fi4q=find(fi>3*pi/2);
 if Rx==Ry
  dap=0;
  R00=Rx;
 else
%  dap=(Rx-Ry)/sqrt(Rx*Ry)/2;
  rapax=Rx/Ry;
  dap=-(1-rapax)/(1+rapax)*(1+Delta);
  R00=Rx/(1+Delta+dap);
 end
%  rocf=cce+sqrt(Rx*Ry)*(1+dap*cos(2*fi)+Delta*cos(4*fi)).*exp(j*fi);
%  rocf=cce+R00*(1+dap*cos(2*fi)+Delta*cos(4*fi)).*exp(j*fi);
  NAz=3;
%  'qui da aggiustare per PG', keyboard
  NAz=2;
  if NAz==3
   Delta=0;
   Ry=Rx*.8;
   rapax=Rx/Ry;
   dap=-(1-rapax)/(1+rapax)*(1+Delta);
   R00=Rx/(1+Delta+dap);   
  end
  rocf=cce+R00*(1+dap*cos(NAz*fi)+Delta*cos(4*fi)).*exp(j*fi);
  
end

%figure, plot(rocf), axis equal, pausak
%figure, plot(abs(rocf)),  pausak
%figure, plot(real(rocf),imag(rocf))
%keyboard

