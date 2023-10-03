 function [Gacrit,Trcrit]=gam_crit(Tstor,Ga1,P,ficri);
 
 TscaU
 
 
 s=size(Tstor);
 Mn=1;
 for k=1:ficri-1
  Mn=Tstor(:,:,k)*Mn;
 end
 
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

   
   Id=eye(length(Ga1));

%   Gacritu=s22+s21*inv(Id-Gadu*s11)*Gadu*s12;
   
   
   
   Gacrit=s22u+s21u * inv(Id-Gacritu*s11u) *Gacritu*s12u;

   iM=inv(M11);   
   Tcritu=M22-M21*iM*M12;   
   Gcritu=M21*iM;   
   Trcrit=Tcritu * inv(Id-s22u*Gcritu) *s21u;
   
   
'gacrit', 
%keyboard