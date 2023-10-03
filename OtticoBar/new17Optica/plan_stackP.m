function Ti=plan_stackP(kr,Pstack)

 ini=Pstack.in; 
 fin=Pstack.fin; 
 n=Pstack.n;
 Li=Pstack.Li;
 rep=Pstack.rep;
 
 
 Pstacki=Pstack;

' in stak', keyboard
for klo=1:length(ini)

 for ilo=ini(klo):fin(klo)
  ii=ilo;
  uu=ii;
  pun=[ii:uu];
  nl=n(pun);
  Lil=Li(pun);
  repi=1;
  Pstacki.n=nl;
  Pstacki.Li=Lil;
  Pstacki.rep=repi;  
  Ti(:,:,ilo)=plan_stack0(kr,Pstacki);

 end 

end % klo
'fine Tr mir plan', keyboard
