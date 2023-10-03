
  inval=find(Gsov<=0);
%  Gdu=Gsov.*M2v;
  Gdu=Gsov;
  Gdu(inval)=1e100;
  [du,idu]=sort(Gdu);
  ival=find(du<1e100);
  iso=idu(ival);
return
kco=0;
for kis=iso                                             
  kco=kco+1;					
  Plot.Ap{kco}=Plots.Ap{kis};		
  Plot.ApQ{kco}=Plots.ApQ{kis};		
  Plot.parmod{kco}=Plots.parmod{kis};	
  Plot.XP{kco}=Plots.XP{kis};		
  Plot.YP{kco}=Plots.YP{kis};		
  Plot.X{kco}=Plots.X{kis};			
  Plot.Y{kco}=Plots.Y{kis};			
  Plot.E2xo{kco}=Plots.E2xo{kis};		
  Plot.E2xp{kco}=Plots.E2xp{kis};		
  Plot.FF{kco}=Plots.FF{kis};			
  Plot.E2xp{kco}=Plots.E2xp{kis};		
  if iLP==0
   Plot.E2yo{kco}=Plots.E2yo{kis};		
   Plot.E2yp{kco}=Plots.E2yp{kis};		
   Plot.E2zo{kco}=Plots.E2zo{kis};		
   Plot.E2zp{kco}=Plots.E2zp{kis};		
  end 
  Plot.Ef{kco}=Plots.Ef{kis};			
  Plot.Cug{kco}=Plots.Cug{kis};		
  Plot.gou{kco}=Plots.gou{kis};		
  Plot.aou{kco}=Plots.aou{kis};		
  Plot.fou{kco}=Plots.fou{kis};		
  Plot.ze{kco}=Plots.ze{kis};			
  Plot.gg0{kco}=Plots.gg0{kis};		
  Plot.FF{kco}    =Plots.FF{kis};
  Plot.parmod{kco}=Plots.parmod{kis};

 end 					

 