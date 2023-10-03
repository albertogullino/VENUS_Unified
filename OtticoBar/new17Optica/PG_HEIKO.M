iLP=0;
Kp=linspace(0,.3,301)';
clear GGev GGmv Tev Tmv

nstv=[0:5:30];
nstv=[0];
for ist=1:length(nstv)
nstratiu=nstv(ist)
%   [GGe1,GGm1,TTe1,TTm1]=gaemms(KK,freq,lambda,Lvtr,nvtr,...
%                     Lbt,nbt,nstratiu,rfu,rr,iLP,Luv,nuv);

%pausak
Lplit=[];
nplit=[];
Lplut=[];
nplut=[];
   [GGe1,GGm1,Te1,Tm1]=gaemms(Kp,freq,lambda,Lplit,nplit,Lplt,...
              real(nplt),nstratiu,rfu,rr,iLP1,Lplut,nplut,rqw);

  GGev(:,ist)=GGe1;
  GGmv(:,ist)=GGm1;
  Tev(:,ist)=Te1;
  Tmv(:,ist)=Tm1;
end
figure,
plot(Kp,abs(GGev)), hold on, plot(Kp,abs(GGmv),'--')
pausak
figure,
plot(Kp,abs(Tev)), hold on, plot(Kp,abs(Tmv),'--')
figure
plot(Kp,abs(Tev)./abs(Tmv))
