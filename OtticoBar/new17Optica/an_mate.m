if iLP==1
 npias=numodi;
else
 npias=2*numodi;
end
clear Pias1
M1=max(max(abs(Kplot1)));
for ir=1:npias
 puri=ldap(ir)+1:ldap(ir+1);
 for ic=1:npias
  puco=ldap(ic)+1:ldap(ic+1);
  ' dentro', pausak
  m=mean(abs(Kplot1(puri,puco)));
  if m/M1<1e-6
   Pias1(ir,ic)=0;
  else
   Pias1(ir,ic)=1;
  end
 end
end

clear iacc
clear iaccd iaccp

for ir=1:npias
 iacc=find(Pias1(ir,:)==1);
 if length(iacc)~=npias & length(iacc)~=0
  iaccd(ir,1:length(iacc))=iacc;
 end
end


if exist('iaccd')
 [su,iso]=sort(iaccd(:,1));
 iaccd=iaccd(iso,:);
 fiv=find([1; diff(iaccd(:,1))]~=0);
 iaccpdu=iaccd(fiv,:);
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
iacc=iaccp;
