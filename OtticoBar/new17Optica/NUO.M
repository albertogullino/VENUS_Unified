clear
close all
load nuo1
pz=11;
C=cumsum(yiT(:,pz));

figure,
plot(ro_inT,T(1:end-1,1)),
sy=size(yiT);
pu=(1:sy(1))+(pz-1)*sy(1);
hold on, plot(adis(pu)/kcav0,C+Tdis(end,end),'.'), pausak

clear KTempp
h=figure;
pu=1:32;
   s=size(yiT);
   isuT=[1:s(1)];
%   for iz=1:s(2)
   for iz=s(2)
    Ktoiiep=0;
    Ktoiiem=0;
    indis=0;
    for ndis=fliplr(isuT)
     indis=indis+1;
     if iresK==0
      ndis1=ndis;
     else
      ndis1=ndis+(iz-1)*s(1);
     end
     [adis(ndis1)/kcav0, yiT(ndis,iz)]
%     map(Kosp(:,:,ndis1)), pausak
     map(Kosp(:,:,ndis1),pu,pu,h), pausak
     Ktoiiep=Ktoiiep+Kosp(:,:,ndis1)*yiT(ndis,iz);
     Ktoiiem=Ktoiiem+Kosm(:,:,ndis1)*yiT(ndis,iz);
    end  %ndis
    KTempp(:,:,iz)=Ktoiiep;
    KTempm(:,:,iz)=Ktoiiem;
   end
%  ' qui KTE', keyboard
