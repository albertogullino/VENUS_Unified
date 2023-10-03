  freq=0;
  Lvtr=[];
  nvtr=[];
  Lbt=[];
  nbt=[];
  Luv=[];
  nuv=[];
  fapeu=0;
  nstratiu=0;
  iLP=0;
  rfu=1;
  rr=3;
  ring=1.52;
  lambda=1.5;
  k0=2*pi/lambda;
  KKs=linspace(0,.35,30);
  
  [GGe1,GGm1,TTe1,TTm1]=gaperd(KKs,freq,lambda,Lvtr,nvtr,...
                     Lbt,nbt,nstratiu,rfu,rr,iLP,Luv,nuv,ring);
                     
figure, plot(KKs,GGe1,KKs,GGm1)   , pausak                  
figure, plot(KKs,abs(GGe1),KKs,abs(GGm1))