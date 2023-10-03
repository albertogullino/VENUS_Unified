 keyboard
 if iloo(1)>1
    ilopl=1;
    Tdu=Tstack(:,:,ilopl);
 else 
     ilopl=0;
    Tdu=1;
 end   
  
 klo=0; 
 rep=fmlsi(iloo,2);
 nstmir=fmlsi(iloo,1);
 while klo<length(iloo)
  klo=klo+1;
%  ind=iloo(klo);
  nsti=nstmir(klo)
  pausak
   if nsti==1
%  klo, pausak
    Tdu=Tstor(:,:,klo)*Tdu;
   else
%     klo, pausak
    Tmir=1;
    NPair=rep(klo);
    pumir=[klo:klo+nsti-1];
    klos=klo;
    for pum=pumir
%     ind=iloo(pum);
     ind=pum;
     Tmir=Tstor(:,:,ind)*Tmir;
    end
    Tmir=Tmir^NPair;
    klo=klo+nsti-1;
    Tdu=Tmir*Tdu;
%    'mir', keyboard
   end
   if length(find(iloo(klo)==fipla))==1
    ilopl=ilopl+1;
    Tdu=Tstack(:,:,ilopl)*Tdu;
   end
  end   
  