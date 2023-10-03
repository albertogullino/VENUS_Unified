pui=[1:24 45:68];
Gadu=Gad(pui);
   dRi=diag(Gadu)*Mn12(pui,pui)-Mn22(pui,pui);
   dPi=Mn21(pui,pui)-diag(Gadu)*Mn11(pui,pui);
   Geqm=inv(dRi)*dPi;
   Mns=Mn;
   Mis=Mi;
    if  ifp==-10 %& ifiez>0
     Geq=diag(Geqm);
     hpro=figure; plot(abs(Geq)), pausak
    end 