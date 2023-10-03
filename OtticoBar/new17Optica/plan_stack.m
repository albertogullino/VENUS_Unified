function Ti=plan_stack(kr,Pstack)

 ini=Pstack.in; 
 fin=Pstack.fin; 
 n=Pstack.n;
 Li=Pstack.Li;
 rep=Pstack.rep;

 
 
 Pstacki=Pstack;

%' in plan_stak', keyboard

 for ilo=1:length(ini)
  ii=ini(ilo);
  uu=fin(ilo);
  pun=[ii:uu];
  nl=n(pun);
  Lil=Li(pun);
  repi=rep(pun);
  Pstacki.n=nl;
  Pstacki.Li=Lil;
  Pstacki.rep=repi;  
  Ti(:,:,ilo)=plan_stack0(kr,Pstacki);
% 'plan_stack', keyboard
 end 

%'fine plan_stack', keyboard