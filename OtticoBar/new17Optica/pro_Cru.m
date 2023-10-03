%load saCr
%load saC
load sa4

icrit(21:end)=0;

ficri=find(icrit==1);

% [Gacritn,Trcrit]=Gam_critU(Tstor,Ga1,Mcrit,ficri);
 [Gacritn,Tn]=Gam_critU(Tstor,Ga1,Mcrit,ficri);
 
map(log10(Gacritn)), pausak
map(log10(Tn)), pausak

return

  lcri=length(ficri);

 Mn=1;
 for k=1:lcri
  Mn=Tstor(:,:,k)*Mn;
  de=det(Mn);
 end
 s=size(Mn);
  l1=s(1)/2;
  l2=s(1)/2+1;
  l3=s(1);
  M11=Mn(1:l1,1:l1);
  M12=Mn(1:l1,l2:l3);
  M21=Mn(l2:l3,1:l1);
  M22=Mn(l2:l3,l2:l3);
iT=inv(M22);
s11=-iT*M21;
s12=iT;
s21=M11-M12*iT*M21;
s22=M12*iT;

   Gadu=Ga1;
   
   A=M21*Gadu+M22;
   B=M11*Gadu+M12;
   inva=inv(A);
   Gacritu=B*inva; 