 load checkCampiEH_1D_rr5

 h=figure;
 set(h,'pos',[383    35   771   901])
 subplot(211)
 plot(z/1000,log10(abs(V/V(1)).^2))
 title(['rr = 5:        Gth mio =',num2str(g0,5),'   Gth Alb =',num2str(Gth,5)])
 ylabel('|Campo elettrico|^2  normalizzato')
 hold on

 plot(zet+4,log10(abs(Ez/Ez(1)).^2),'r.'),
 subplot(212)
 plot(z/1000,log10(abs(I/I(1)).^2))
 title(['Lam mio =',num2str(la0,6),'  Lam Alb =',num2str(LambdaMode/1000,6)])
 ylabel('|Campo magnetico|^2  normalizzato')
 xlabel(' long. coordinate (um)')
 hold on
 plot(zet+4,log10(abs(Hz/Hz(1)).^2),'r.'),
pausak

 load checkCampiEH_1D_rr1
 
  h=figure;
  set(h,'pos',[383    35   771   901])
  subplot(211)
  plot(z/1000,log10(abs(V/V(1)).^2))
  title(['rr = 1:         Gth mio =',num2str(g0,5),'   Gth Alb =',num2str(Gth,5)])
  ylabel('|Campo elettrico|^2  normalizzato')
  hold on
 
  plot(zet,log10(abs(Ez/Ez(1)).^2),'r.'),
  subplot(212)
  plot(z/1000,log10(abs(I/I(1)).^2))
  title(['Lam mio =',num2str(la0,6),'  Lam Alb =',num2str(LambdaMode/1000,6)])
  ylabel('|Campo magnetico|^2  normalizzato')
  xlabel(' long. coordinate (um)')
  hold on
 plot(zet,log10(abs(Hz/Hz(1)).^2),'r.'),