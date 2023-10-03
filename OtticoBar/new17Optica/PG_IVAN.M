close all
iLP=1;
Kp=0;
clear GGev GGmv Tev Tmv

nstratiu=3;

% 1                           |    1.0
%GaAs            |         0.06949 microns   |    3.53
%InGaAs         |         0.07963 microns   |    3.08
%GaAs            |         0.06949 microns   |    3.53
%InGaAs         |         0.07963 micron    |    3.08
%GaAs            |         0.06949 microns   |    3.53
%InGaAs         |         0.07963 microns   |    3.08
%GaAs            |         0.06949 microns   |    3.53
%Air                |         1                           |    1.0


Lplit=[];
nplit=[];
Lplut=[0.06949];
nplut=[3.53];
Lplt=[0.06949; 0.07963];
nplt=[3.53; 3.08];
freq=0;
rr=1;
rfu=1;


lav=linspace(0.75,1.25,200);

 for it=1:length(lav)
   lambda=lav(it);
   [GGe1,GGm1,Te1,Tm1]=gaemms_h(Kp,freq,lambda,Lplit,nplit,Lplt,...
              nplt,nstratiu,rfu,rr,iLP,Lplut,nplut,rr);
   Gv(it)=GGe1;

end

figure,
plot(lav,abs(Gv).^2),
