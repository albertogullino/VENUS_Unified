
lav=linspace(.8,1.,200);


fr=0
k=0

for kl=1:length(lav)

 lai=lav(kl);
 [GGe2,GGm2,TTe2,TTm2]=gaperd(k,fr,lai,Lvbr,nvbr,Lbb,nbb,nstratid,rfd,rr,iLP,Luvb,nuvb,ring,0);
 
 Ga(kl)=GGe2;
end

figure, plot(lav,abs(Ga.^2))