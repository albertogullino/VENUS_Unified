RAD='rada';
IFP=-4;
ISavetutto=0;
Dp=0:40:400;
NP=11;
for knp=1:NP
  eval([' load ',RAD,num2str(knp)])
  ifp=IFP;
  isavetutto=ISavetutto;
 postg_sub
 [Ga,ig]=min(gsov/vg);
 F=fsov(ig)*lambda*1000;
 La=(1+fsov(ig))*lambda*1000;
 Gtot(knp)=Ga;
 Ftot(knp)=F;
 Ltot(knp)=La;
end 

h=figure,
set(h,'pos',[45   156   464   706])
subplot(311)
plot(Dp,Ftot,'linewidth',2)

subplot(312)
plot(Dp,Ltot,'linewidth',2)

subplot(313)
semilogy(Dp,Gtot,'o-','linewidth',2)
pausak

h=figure,
set(h,'pos',[45   156   464   706])
subplot(311)
plot(Dp,Ltot,'linewidth',2)
 ylabel(' Wavelength ')

subplot(312)
semilogy(Dp,Gtot,'o-','linewidth',2) 
 ylabel(' gain/qw  (1/cm) ')

 subplot(313)
 plot(Dp,Gtot,'o-','linewidth',2) 
%  xlabel(YL)
 ylabel(' gain/qw  (1/cm) ')