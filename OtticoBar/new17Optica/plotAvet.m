   fa=figure;
   set(fa,'pos',[163         207        1650        771])
   
   subplot(211)
   puinp=puini(2:end)-1;
   plot(Ap,'g.-'), hold on, plot(ApQ,'w'),        plot(puinp,Ap(puinp),'ro'),
   subplot(212)
   plot(log10(abs(Ap)),'g.-'), hold on, plot(log10(abs(ApQ)),'w'),        plot(puinp,log10(abs(Ap(puinp))),'ro'),
   
return

hold on
   plot(Ap,'r.-'), hold on, plot(ApQ,'y'),        