
%'% riduco la base'
%keyboard
firid=find(shailoop==5);
if length(firid>0)
 irid_bas=1;
else
 irid_bas=0;
end
 irid_bas=1;
% irid_bas=0;

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
 pre=.01;
 Pus=[];
 ldap=[0];
 if iLP==1
  fiLP=1;
  npam=1;
 else
  fiLP=2;
  npam=2;
 end
 if fiLP==2 & is_even(numodi)==0
  npam=1;
 end

 for nm=1:npam:numodi*fiLP
  pui=pu+(nm-1)*lk;
  Du=abs(Dm(pui));
  Du=Du/max(Du);
  der=diff([0; Du]);
  fiu=find(Du-pre>0 & der>0);
  pua=fiu(1):lk;
  if npam==1
   Pus=[Pus pui(pua)];
   ldap=[ldap ldap(end)+length(pua)];
  else
   Pus=[Pus pui(pua) pui(pua)+lk];
   ldad=ldap(end)+length(pua);
   ldap=[ldap ldad ldad+length(pua)];
  end
 end
% ifp=1
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

 lKA=si2(1);
 Pus=1:lKA;
 numoa=0:numodi-1;
% ld=[numoa numodi+numoa]*nk1max;
% ldap=[ld ld+lKA];
 ld=[numoa numodi+numoa]*nk1max;
 ldap=[ld ld(end)+nk1max];

end

Pusa=1:si2(1);
Pust=[Pus Pus+si2(1)];

% set global IdeOo pMu0u pMu0u1 lKA nk1max nures sim0

lKAn=length(Pus);


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
 IdeOon=diag(ones(2*lKAn,1));


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


disp(' RID_BAS'),
[nures/2*length(KK) length(Pus)]
%keyboard
%disp(' IdeOo')
%keyboard

