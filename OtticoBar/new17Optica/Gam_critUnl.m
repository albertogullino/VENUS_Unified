% function [Gacrit,Trcrit]=gam_critU(Tstor,Ga1,Mcri,ficri);
 function [Gacritu,Trcritu,Trc,Grc]=gam_critUnl(Tstor,Ga1,Mcri,ficri,fmlsi,ifiez,ficrip);

%'dentro', keyboard

if ifiez>=1 | length(ficrip)==1

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

else  %ifiez
%'sono in OK', keyboard
 finor=1:ficrip(1)-1;
 fiend=ficrip(1)-1;
% calcolo gamma_in in modo efficiente, nel modo solito
 jsa=0;
 Mt=1;
 while jsa<fiend
  jsa=jsa+1
%  pausak
  irep=fmlsi(jsa,2); 
  if irep==1
   Mp=Tstor(:,:,jsa);  
%   'no rep'
  else
   nrep=fmlsi(jsa,1); 
   Mdu=1; 
   jsa=jsa-1;
%   ' rep', keyboard
   
   for jrep=1:nrep
    jsa=jsa+1
    Mn=Tstor(:,:,jsa);  
    Mdu=Mn*Mdu;
%    'jsa in rep', pausak
   end
   Mp=Mdu^irep;
%   ' potenza', pausak
  end 
  Mt=Mp*Mt;
 end
%'cont', keyboard
 Trc=0;
 Grc=0;
 [Gai,Tr1]=GaTu(Mt,Ga1);
 Ga_s=Ga1;
 Ga1=Gai;
 lcrit=ficrip(1):ficrip(2);
 Tr1=1;
 for kt=lcrit
  Mn=Tstor(:,:,kt);
  [Ga1,Tr1d]=GaTu(Mn,Ga1);
  Tr1=Tr1*Tr1d;
%  Trd(:,:,k)=Tr1;
%  Grc(:,:,k)=Ga1;
 end   
 Gacritu=Ga1;
 Trcritu=Tr1;
end

%'fine', keyboard