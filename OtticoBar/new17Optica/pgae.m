la=1;

lav=linspace(.5,1.5,200);

n1=2
n2=4
nm=[n2 n1];
lm=la/4./nm;
NP=5;
ni=1
no=1
rr=1
iLP=0
Lu=lm(1);
nu=nm(1);

fr=0
k=0

for kl=1:length(lav)

 lai=lav(kl);
 [GGe2,GGm2,TTe2,TTm2]=gaperd(k,fr,lai,[],[],lm,nm,NP,no,rr,iLP,Lu,nu,ni,0);
 Ga(kl)=GGe2;
end

figure, plot(lav,abs(Ga.^2))