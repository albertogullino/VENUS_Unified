% function [Gacrit,Trcrit]=gam_critU(Tstor,Ga1,Mcri,ficri);
% funzione modificata May 16 per riflettivitā come sommatoria contributi

function [Gacritu,Trcritu]=gam_TILT(Ga1,Mcri,ficrin,fmlsi,dovesono,ipol);

ficri=abs(ficrin);
iPv=sign(ficrin);

lc=length(ficri);
firep=find(fmlsi(ficri,2)>1);
fire=ficri(firep);
if length(fire)>0
 re=fmlsi(fire,2);
 puesp=[ficri(1:firep(1)-1)' repmat(fire',1,re(1)) ficri(firep(end)+1:end)'];
else
 puesp=ficri';
 
end

 %' sono dentro  gam_crit', keyboard
 
 Ga1_sa=Ga1;
% TscaU
 
 lcri=length(puesp);
 s=size(Ga1);
 Mn=1;
 Tr1=1;


kcrit=puesp;
if dovesono==2
 pacri=-1;
 kcrit=puesp(end:-1:1);
end

%'prima di crit', keyboard
klast=kcrit(end-1:end);
kcar=kcrit(1:end-1);
% for k=1:lcri
ick=0;
 for k=kcar
  ick=ick+1;
%  iP=iPv(ick);
  iP=iPv(k);
  kt=k;

  Mn=Mcri{kt,ipol};
    Ga1_s=Ga1;
%  kt, pausak  
   [Ga1,Tr1d]=GaScattHCG(Mn,Ga1,kt,iP);
   Tr1=Tr1*Tr1d;
%   Tr1=Tr1d;
   Trd(:,:,k)=Tr1d;
   Grc(:,:,k)=Ga1;
 end
 
  Trcritu=Tr1; 
  Gacritu=Ga1; 
  
 return
 
  ic=0;
  for k=klast
   ic=ic+1;
 %  iP=iPv(ick);
   iP=iPv(k);
   kt=k;
 
   Mn=Mcri{kt,ipol};
     Ga1_s=Ga1;
 %  kt, pausak  
    [Ei,Si]=GaScatt_Tilt(Mn,Ga1,kt,iP);
    if ic==1
     E=Ei;
    end
    Tr1=Tr1*Tr1d;
 %   Tr1=Tr1d;
 end
 

 
 sab11=E*sb11*E;
 sab12=E*sb12;
 sab21=sb21*E;
 sab22=sb22;
 
 
 Gu=Ga1;
 In=inv(eye(size(Gu))-Gu*s11u);
 Gae=s22u+s21u*In*Gu*s12u;
 
 Par=(eye(size(Gu))-s11u*Gu);
 Tre=inv(Par)*s12u;


 
 Tpre=1;

% for k=lcri:-1:1
%  Tri=Trd(:,:,k);
%  Trc(:,:,lcri-k+1)=Tri;
% end 
Trc=Trd;
Trcritu=Tr1; 
Gacritu=Ga1; 
% 'fine gacritu', keyboard 
return 

% 'fine for', keyboard
  s=size(Mn);
 if s(1)>1
  s=size(Mn);
  l1=s(1)/2;
  l2=s(1)/2+1;
  l3=s(1);
  M11=Mn(1:l1,1:l1);
  M12=Mn(1:l1,l2:l3);
  M21=Mn(l2:l3,1:l1);
  M22=Mn(l2:l3,l2:l3);
%iT=inv(M22);
%s11=-iT*M21;
%s12=iT;
%s21=M11-M12*iT*M21;
%s22=M12*iT;

   Gadu=Ga1;
   
   A=M21*Gadu+M22;
   B=M11*Gadu+M12;
   inva=inv(A);
   Gacritu=B*inva;
%   Trcritu=inv(M22);
   Trcritu=Tr1*inva;
else
 Gacritu=Ga1;
 Trcritu=Tr1;
end

%map(log10(Gacritu))
%' contro ', keyboard   
   
   return
   
   Id=eye(length(Ga1));

%   Gacritu=s22+s21*inv(Id-Gadu*s11)*Gadu*s12;
   
   
   
   Gacrit=s22u+s21u * inv(Id-Gacritu*s11u) *Gacritu*s12u;

   iM=inv(M11);   
   Tcritu=M22-M21*iM*M12;   
   Gcritu=M21*iM;   
   Trcrit=Tcritu * inv(Id-s22u*Gcritu) *s21u;
   
   
'gacrit', 
keyboard