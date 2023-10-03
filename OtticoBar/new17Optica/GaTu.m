function [Gae,Tre]=GaTu(Mn,Ga1);

  s=size(Mn);
  l1=s(1)/2;
  l2=s(1)/2+1;
  l3=s(1);
  M11=Mn(1:l1,1:l1);
  M12=Mn(1:l1,l2:l3);
  M21=Mn(l2:l3,1:l1);
  M22=Mn(l2:l3,l2:l3);

% iT=inv(M22);
% s11=-iT*M21;
% s12=iT;
% s21=M11-M12*iT*M21;
% s22=M12*iT;

   A=M21*Ga1+M22;
   B=M11*Ga1+M12;
   inva=inv(A);
   Gae=B*inva;
%   Tre=inv(M22);   
   Tre=inva;   
   ' TRE', keyboard
return   

   iM=inv(M11);   
   Tcritu=M22-M21*iM*M12;   
   Gcritu=M21*iM;   
   Trcrit=Tcritu * inv(Id-s22u*Gcritu) *s21u;   