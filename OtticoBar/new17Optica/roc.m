  Base=2;
  Lbuf=5;
  x=linspace(.3,1.5,20);
  r_lens=x;
  
  l_minore=Base-r_lens;
  angolo=Lbuf./l_minore;
  Rc=r_lens.*sqrt(1+1./angolo.^2)   
  
  figure, plot(x,Rc)