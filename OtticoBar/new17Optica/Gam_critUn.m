% function [Gacrit,Trcrit]=gam_critU(Tstor,Ga1,Mcri,ficri);
function [Gacritu,Trcritu,Trc,Grc]=gam_critUn(Tstor,Ga1,Mcri,ficri,fmlsi);

lc=length(ficri);
fire=find(fmlsi(1:lc,2)>1);
if length(fire)>0
 re=fmlsi(fire,2);
 re=[re(1); re];
 fim=find(diff(re)~=0);
 fim=[1 fim length(fire)];
 fi1=fire(1);
 if fi1>1
  puesp=[1:fi1-1];
 else
  puesp=[];
 end
 for k=1:length(fim)-1
  ind=fim(k):fim(k+1);
  pu=fire(ind);
  pua=repmat(pu,fmlsi(pu(1),2),1);
  puesp=[puesp pua'];
 end
 fie=fire(end);
 if fie<lc
  puesp=[puesp fie+1:lc];
 end
else
 puesp=ficri;
end

% ' sono dentro ', keyboard
 
 Ga1_sa=Ga1;
% TscaU
 
 lcri=length(puesp);
 s=size(Tstor);
 Mn=1;
 Tr1=1;
 for k=1:lcri
  kt=puesp(k);

  Mn=Tstor(:,:,kt);
    Ga1_s=Ga1;
   [Ga1,Tr1d]=GaTu(Mn,Ga1);
   Tr1=Tr1*Tr1d;
   Trd(:,:,k)=Tr1;
   Grc(:,:,k)=Ga1;
 end
 Tpre=1;

 for k=lcri:-1:1
  Tri=Trd(:,:,k);
%  Tiv=Tri*Tpre;
  Trc(:,:,lcri-k+1)=Tri;
 end 
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
%keyboard