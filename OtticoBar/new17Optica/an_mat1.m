
%' anmat1'
global ldapu

ldap=ldapu;
%if iLP==1
% npias=numodi;
%else
% npias=2*numodi;
%end
npias=(length(ldap)-1);
%npias, keyboard

clear Pias1
M1=max(max(abs(Kplot1)));
for ir=1:npias
 puri=ldap(ir)+1:ldap(ir+1);
 for ic=1:npias
  puco=ldap(ic)+1:ldap(ic+1);
  m=mean(abs(Kplot1(puri,puco)));
  if m/M1<1e-6
   Pias1(ir,ic)=0;
  else
   Pias1(ir,ic)=1;
  end
  if ic==ir
   Pias1(ir,ic)=1;
  end
 end
end

clear iaccp
clear iaccpd

for ir=1:npias
 iacc=find(Pias1(ir,:)==1);
 if length(iacc)~=npias & length(iacc)~=0
  iaccpd(ir,1:length(iacc))=iacc;
 end
end

if exist('iaccpd')
 [su,iso]=sort(iaccpd(:,1));
 iaccpd=iaccpd(iso,:);
 fiv=find([1; diff(iaccpd(:,1))]~=0);
 iaccpdu=iaccpd(fiv,:);
 sa=size(iaccpdu);
 iaccp=[];
 for ka=1:sa(1)
  if length(find(iaccpdu(ka,:)~=0))>0
   iaccp=[iaccp; iaccpdu(ka,:)];
  end
 end
else
 iaccp=[];
end

  nktot=length(Kplot1);
  VPU=1:nktot^2;
  MPU=reshape(VPU,nktot,nktot);
%  ' mat1 ', keyboard

if length(iaccp)==0
 PUm=[];
else
 for irid=1:length(iaccp(:,1))
 iadu=iaccp(irid,:);
 iadu1=iadu(find(iadu>0));
 iaccdu=iadu1;
 pMdu=[];
 for ir=iaccdu
  for ic=iaccdu
   if ic>0 & ir>0
    Mdu=MPU;
    Mdu(ldapu(ir)+1:ldapu(ir+1),ldapu(ic)+1:ldapu(ic+1))=NaN;
    pMpd=find(isnan(Mdu)==1);
%    size(pMpd), pausak
    pMdu=[pMdu; pMpd];
   end
  end
 end
  PUm{irid}=sort(pMdu);
 end
end

if length(PUm)==0
 PUm=1;
end
