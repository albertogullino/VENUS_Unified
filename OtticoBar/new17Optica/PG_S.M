load pgam
close all
iLP=0;
Kp=linspace(0,.4,301)';
clear GGev GGmv Tev Tmv

nstv=[0:5:25];
%nstv=[100];

%Lplit=[];
%nplit=[];
%Lplut=[];
%Lplut(3:5)=0;
%nplut=[];
%Lplt=[0.069; 0.060];
%nplt=[3.044; 3.5155];
%lambda=2*real(sum(Lplt.*nplt))
%pausak
ireal=1;
if ireal==1
 nplt=real(nplt);
 nplut=real(nplut);
 nplit=real(nplit);
end


for ist=1:length(nstv)
nstratiu=nstv(ist)
%   [GGe1,GGm1,TTe1,TTm1]=gaemms(KK,freq,lambda,Lvtr,nvtr,...
%                     Lbt,nbt,nstratiu,rfu,rr,iLP,Luv,nuv);

%pausak

   [GGe1,GGm1,Te1,Tm1]=gaemms_h(Kp,freq,lambda,Lplit,nplit,Lplt,...
              nplt,nstratiu,rfu,rr,iLP1,Lplut,nplut,real(rqw));

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
plot(Kp,abs(Tmv).^2-abs(Tev).^2)
figure
semilogy(Kp,abs(GGev).^2-abs(GGmv).^2)
