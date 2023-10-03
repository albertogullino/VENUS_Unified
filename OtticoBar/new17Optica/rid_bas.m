
% riduco la base
firid=find(shailoop==5);
if length(firid>0)
 irid_bas=1;
else
 irid_bas=0;
end

if irid_bas==1
 asce=axtot+aytot;
 sa=size(asce);
 if min(sa)>1
  [du1,ips1]=max(asce);
  [du2,ips2]=max(du1);
  irim=ips1(ips2);
  icom=ips2;
 else
  [du1,icom]=max(asce);
  irim=1;
 end
 Kplot1=MKosp{icom};
 Kplot2=MKosm{icom};
 siK=size(Kplot1);
 if length(siK)>2
  Kplot1=reshape(Kplot1(:,:,irim),siK(1:2));
 end
 siK=size(Kplot2);
 if length(siK)>2
  Kplot2=reshape(Kplot2(:,:,irim),siK(1:2));
 end
 Dm=diag(Kplot1+Kplot2);
 if ifp>=1
  figure, plot(Dm),
 end

 lk=length(KK);
 pu=1:lk;
 pre=.05;
 Pus=[];
 ldap=[0];
 if iLP==1
  fiLP=1;
 else
  fiLP=2;
 end
 for nm=1:numodi*fiLP
  pui=pu+(nm-1)*lk;
  Du=abs(Dm(pui));
  Du=Du/max(Du);
  der=diff([0; Du]);
  fiu=find(Du-pre>0 & der>0);
  pua=fiu(1):lk;
% plot(pu,Du,pu(pua),Du(pua),'r.')
% pausak

  Pus=[Pus pui(pua)];

  ldap=[ldap ldap(end)+length(pua)];
 end

 if ifp>=1
  lD=length(Dm);
  put=1:lD;
  figure,
  plot(put,Dm,put(Pus),Dm(Pus),'r.')
  [length(Dm) length(Pus)]
  map(Kplot2)
  map(Kplot2(Pus,Pus))
 end
else
 Pus=1:si2(1);
 numoa=0:numodi-1;
 ld=[numoa numodi+numoa]*nk1max;
 ldap=[ld ld+lKA];

end
si2=[length(Pus) length(Pus)];

 KAp=KAp(Pus,Pus);
 KAm=KAm(Pus,Pus);

for icom=shailoop

 Kp1=MKosp{icom};
 if iztm==1
  Kp2=MKoszp{icom};
 end
 siK=size(Kp1);
 if length(siK)>2
  for ikm=1:siK(3)
   Kplot1=reshape(Kp1(:,:,ikm),siK(1:2));
   Kpp1(:,:,ikm)=Kplot1(Pus,Pus);
   if iztm==1
    Kplot1=reshape(Kp2(:,:,ikm),siK(1:2));
    Kpp2(:,:,ikm)=Kplot1(Pus,Pus);
   end
  end
 else
   Kpp1=Kp1(Pus,Pus);
   if iztm==1
    Kpp2=Kp2(Pus,Pus);
   end
 end
 MKosp{icom}=Kpp1;
 if iztm==1
  MKoszp{icom}=Kpp2;
 end

 Kp1=MKosm{icom};
 if iztm==1
  Kp2=MKoszm{icom};
 end
 siK=size(Kp1);
 if length(siK)>2
  for ikm=1:siK(3)
   Kplot1=reshape(Kp1(:,:,ikm),siK(1:2));
   Kpp1(:,:,ikm)=Kplot1(Pus,Pus);
   if iztm==1
    Kplot1=reshape(Kp2(:,:,ikm),siK(1:2));
    Kpp2(:,:,ikm)=Kplot1(Pus,Pus);
   end
  end
 else
   Kpp1=Kp1(Pus,Pus);
   if iztm==1
    Kpp2=Kp2(Pus,Pus);
   end
 end
 MKosm{icom}=Kpp1;
 if iztm==1
  MKoszm{icom}=Kpp2;
 end

end




% set global IdeOo pMu0u pMu0u1 lKA nk1max nures sim0

if irid_bas==1

  lKA=length(Pus);
  if iLP==0
   nures=numodi*4;
   sim0=nk1max*2;
  else
   nures=numodi*2;
   sim0=nk1max*2;
  end

  IdeOo=diag(ones(2*lKA,1));

  ld=ldap(1:end-1);
  nk1vd=diff(ldap);
  nk1v=[nk1vd nk1vd];
  ldu=[ld ld+lKA];
  pMu0u=[];
  for kdu=1:nures
   pMu0u=[pMu0u ldu+(kdu-1)*2*lKA*nk1v(kdu)];
  end

  pMu0u1=[0 2*lKA^2; lKA lKA+2*lKA^2];


  ld1=1:nk1max;
  if iLP==0
   ldud0=[ld1];
   ldud=[ldud0 ldud0+lKA];
  else
   ldud0=[ld1 ld1+lKA];
   ldud=ldud0;
  end

  pMc=[];
  for kdu=ldud
   pMc=[pMc ldud+(kdu-1)*lKA*2];
  end

  ldud=ld1;
  pMei=[];
  for kdu=ldud
   pMei=[pMei ldud+(kdu-1)*lKA];
  end

else


 if iLP==0
  lKA=nk1max*numodi*2;
  nures=numodi*4;
  sim0=nk1max*2;
 else
  lKA=nk1max*numodi;
  nures=numodi*2;
  sim0=nk1max*2;
 end

 IdeOo=diag(ones(2*lKA,1));


  numoa=0:numodi-1;
  ld=[numoa numodi+numoa]*nk1max;
  ldu=[ld ld+lKA];
  pMu0u=[];
  for kdu=1:nures
   pMu0u=[pMu0u ldu+(kdu-1)*2*lKA*nk1max];
  end

  pMu0u1=[0 2*lKA^2; lKA lKA+2*lKA^2];


  ld1=1:nk1max;
  if iLP==0
   ldud0=[ld1];
   ldud=[ldud0 ldud0+lKA];
  else
   ldud0=[ld1 ld1+lKA];
   ldud=ldud0;
  end
  pMc=[];
  for kdu=ldud
   pMc=[pMc ldud+(kdu-1)*lKA*2];
  end

  ldud=ld1;
  pMei=[];
  for kdu=ldud
   pMei=[pMei ldud+(kdu-1)*lKA];
  end

end

Pusasf=Pus;
disp(' RID_BAS'),
[nures/2*length(KK) length(Pus)]
%keyboard
%disp(' IdeOo')
%keyboard

