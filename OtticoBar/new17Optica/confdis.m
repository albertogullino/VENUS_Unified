load old
 fcov=1000*lambda;
  h=figure, 
  subplot(211), plot(fcov*fou',aou','o'), grid
  subplot(212), semilogy(fcov*fou',2*gou'/vgconv,'o'), grid
  
  load rsa
  figure(h), 
    subplot(211), hold on, plot(fcov*fou',aou','--'), 
    subplot(212), hold on, semilogy(fcov*fou',2*gou'/vgconv,'--'), 
  
  load rsam
  figure(h), 
    subplot(211), hold on, plot(fcov*fou',aou','.-'), 
    subplot(212), hold on, semilogy(fcov*fou',2*gou'/vgconv,'.-'),   