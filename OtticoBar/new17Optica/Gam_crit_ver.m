 function Gacrit=gam_crit_ver(Tstor,Ga1,P,icrit);
 
 s=size(Tstor);
 Mn=1;
 for k=1:icrit
  Mn=Tstor(:,:,k)*Mn;
 end
 
  l1=s(1)/2;
  l2=s(1)/2+1;
  l3=s(1);
  M11=Mn(1:l1,1:l1);
  M12=Mn(1:l1,l2:l3);
  M21=Mn(l2:l3,1:l1);
  M22=Mn(l2:l3,l2:l3);

   Gadu=Ga1;
   
   A=M21*Gadu+M22;
   B=M11*Gadu+M12;
   Gacrit=B*inv(A);
   
   
'gacrit_ver', keyboard
%keyboard