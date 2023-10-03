clear
close all
load vec1

clear KTempp
%figure, plot(adis/kcav0,yiT(:,end)), pausak
C=cumsum(yiT(:,end));

figure,
plot(ro_inT,T(1:end-1,1)),
hold on, plot(adis/kcav0,C(end)-C+Tdis(end,end),'.'), pausak
   s=size(yiT);
   isuT=[1:s(1)];
%   for iz=1:s(2)
h=figure;
pu=1:32;
   for iz=s(2)

    Ktoiiep=0;
    Ktoiiem=0;
    indis=0;
    for ndis=isuT
     indis=indis+1;
     if iresK==0
      ndis1=ndis;
     else
      ndis1=ndis+(iz-1)*length(zcut);
     end
     [adis(ndis1)/kcav0, yiT(ndis,iz)]
     map(Kosp(:,:,ndis1),pu,pu,h), pausak
     Ktoiiep=Ktoiiep+Kosp(:,:,ndis1)*yiT(ndis,iz);
     Ktoiiem=Ktoiiem+Kosm(:,:,ndis1)*yiT(ndis,iz);
    end  %ndis
    KTempp(:,:,iz)=Ktoiiep;
    KTempm(:,:,iz)=Ktoiiem;
   end
%  ' qui KTE', keyboard
