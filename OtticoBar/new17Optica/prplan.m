 isemp_plan=1
 
if ifr==1 & pol==pvet(1)

 vfia=(shai'+ai(1,:)').*Li.*(1-anir)*isemp_plan;
  
 
  
  fia=find(vfia~=0);
  iloadd=fia+1;
 
  vsce=[1; fia; length(vfia)];
  dvsce=diff(vsce);
  iloo=[];
  ilins=[];
  for kloo=1:length(dvsce)-1
   if dvsce(kloo)<3
    iloo=[iloo vsce(kloo):vsce(kloo+1)];
   else
    iloo=[iloo vsce(kloo+1)];
   end
    ilins=[ilins iloo(end)];
  end 
 
  ilinsi=ilins+1;
  ilinsu=[ilins(2:end)-1 vsce(end)];
  PstackP.in=ilinsi; 
  PstackP.fin=ilinsu;
  PstackP.n=ni(1,1:end);
  PstackP.Li=Li;
  PstackP.rep=fmlsi(:,1);
  PstackP.nmo=numodiacc+1;
  PstackP.KK=KK;
  PstackP.rr=rr;
 end %ifr
 'per', keyboard
 if pol==pvet(1)
 Tstack=plan_stackP(kcav,PstackP);
 end
 'dopo stack', keyboard
 
 return

 
  Tdu=1;
  
 ilopl=0;
  for klo=1:length(iloo)
   Tdu=Tstor(:,:,klo)*Tdu;
   if length(find(klo==ilins))==1
    ilopl=ilopl+1;
    Tdu=Tstack(:,:,ilopl)*Tdu;
   end
  end 
 
