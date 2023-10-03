
   if iproga==1
     ' qui per gamma', keyboard

     Ga1=diag(Gas);
     s=size(Tc);
     l1=s(1)/2;
     l2=s(1)/2+1;
     l3=s(1);
     T11=Tc(1:l1,1:l1);
     T12=Tc(1:l1,l2:l3);
     T21=Tc(l2:l3,1:l1);
     T22=Tc(l2:l3,l2:l3);

     GAU=(T21+T22*Ga1)/(T11+T12*Ga1);
%     GAU1=(T21+T22*Ga1)*inv(T11+T12*Ga1);
     mapab(GAU), pausak

     D=(diag(GAU).^2);
     Dvte=D(1:end/2);
     Dvtm=D(1+end/2:end);
     Dte=reshape(Dvte,length(KK),numodi);
     Dtm=reshape(Dvtm,length(KK),numodi);
     figure, plot(KK,abs(Dte)), hold on, plot(KK,abs(Dtm),'.-'),
%     figure, plot(KK,abs(Dte),KK,abs(Dtm),'--'), hold on,
%     plot(Kc,abs(Dtec),'.',Kc,abs(Dtmc),'.'),
     ' qui per gamma dopo', keyboard


   end
