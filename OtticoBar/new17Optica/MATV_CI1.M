' matv_ci1'
%iexK=0;
%if exist('KAp')
% KAps=KAp;
% KAms=KAm;
% iexK=1;
%end

 mbv1=[-2+nubesi:nubesu+2];
% bes=besselj(mbv1,KK*a);

 nup1=numodi;
 Kiie=zeros(npk*nup1,npk*nup1);
 Kzer=Kiie;
 Kije=zeros(npk*nup1,npk*nup1);
 Kde=zeros(1,npk*nup1);
 Kd1e=zeros(1,npk*nup1);
 Kiim=zeros(npk*nup1,npk*nup1);
 Kijm=zeros(npk*nup1,npk*nup1);
 Kdm=zeros(1,npk*nup1);
 Kd1m=zeros(1,npk*nup1);

 Kdz=zeros(1,npk*nup1);
 Kiiz=zeros(npk*nup1,npk*nup1);
%%%%%%%%%%%%%%%%%%%%%%%

 clear  Ktoiiei Ktoijei Ktoiimi Ktoijmi
 clear  Ktoiizi Kosp Koszp Kosm Koszm

 for ndis=lxivet
  aloc1=adis(ndis);
%  aloc2=aloc1/kcav0;
  aloc=aloc1;
  bes=besselj(mbv1,KK*aloc1);

   ifidis=find(ifield<0);
   afiel=aitot(ifidis);
%   'afiel VE', keyboard
  if istrumix==1 & igainshape==0
   if (ifii(ndis)<=0 | length(find(afiel==aloc/kcav0))>0)
%  if ifii(ndis)<0
    forma
    scu=size(xcu);
    pug=1:scu(1);
    if exist('Cug')==1
     sC=size(Cug.x);
     sC2=sC(2);
    else
     sC2=0;
    end
    puc=[1:scu(2)]+sC2;
    Cug.x(pug,puc)=xcu;
    Cug.y(pug,puc)=ycu;
    Cug.z(pug,puc)=zcu;
   end
  end
 ip=0;


 for imu=pimu+1
   ipom=(ceil((imu-1)/pasnu));
   ik1v=Pusa([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;
 mu=mbv1(imu);
 jmu=imu;
   if mu==0
    imu0=1;
   end
   if mu~=-20
    fbe=[1 1 1 1];
   else
    if p==1
     fbe=[0 0 1 0];
    elseif p==-1
     fbe=[1 0 0 0];
    end
   end
    for ip1=ik1v
     ip=ip1+(ipom-1)*npk;
%    for ip1=1:npk
%     ip=ip+1;
     Q=KK(ip1)*aloc1;
     P1=bes(ip1,jmu+1)^2+bes(ip1,jmu-1)^2;
     P2=bes(ip1,jmu)*(bes(ip1,jmu+2)+bes(ip1,jmu-2));
     Kdia=.5*(P1-P2);
     Kde(ip)=Kdia*fbe(1);
     Kdm(ip)=Kdia*fbe(3);
     P1=bes(ip1,jmu+1)^2-bes(ip1,jmu-1)^2;
     P2=bes(ip1,jmu)*(bes(ip1,jmu+2)-bes(ip1,jmu-2));
     Kdia=.5*(P1-P2);
     Kd1e(ip)=Kdia*fbe(2);
     Kd1m(ip)=Kdia*fbe(4);
%     for ip2=ip1+1:npk
     for ip2=ip1+1:ik1v(end)
      Q1=KK(ip2)*aloc1;
      P1=bes(ip2,jmu-1)*bes(ip1,jmu-2)+bes(ip2,jmu+1)*bes(ip1,jmu);
      P2=bes(ip1,jmu-1)*bes(ip2,jmu-2)+bes(ip1,jmu+1)*bes(ip2,jmu);
      Kadia=1/(Q1^2-Q^2)*(Q*P1-Q1*P2);
      Kiie(ip,ip+ip2-ip1)=Kadia*fbe(1);
      Kiim(ip,ip+ip2-ip1)=Kadia*fbe(3);
      P1=-bes(ip2,jmu-1)*bes(ip1,jmu-2)+bes(ip2,jmu+1)*bes(ip1,jmu);
      P2=-bes(ip1,jmu-1)*bes(ip2,jmu-2)+bes(ip1,jmu+1)*bes(ip2,jmu);
      Kadia=1/(Q1^2-Q^2)*(Q*P1-Q1*P2);
      Kije(ip,ip+ip2-ip1)=Kadia*fbe(2);
      Kijm(ip,ip+ip2-ip1)=Kadia*fbe(4);
     end
    end
 end %mu

    Ktoiiei(ndis,:,:)=(Kiie+Kiie'+diag(Kde))*aloc^2/4;
    Ktoijei(ndis,:,:)=-(Kije+Kije'+diag(Kd1e))*aloc^2/4;
    Ktoiimi(ndis,:,:)=(Kiim+Kiim'+diag(Kdm))*aloc^2/4;
    Ktoijmi(ndis,:,:)=-(Kijm+Kijm'+diag(Kd1m))*aloc^2/4;

   if iztm==1
    ip=0;
    for imu=pimu
    ipom=(ceil(imu/pasnu));
    ik1v=Pusa([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;
    mu=mbv(imu);
    jmu=imu+1;
%    if mu~=0
    if mu~=-20
      fbe=[1 1 1 1];
     else
      if p==1
       fbe=[0 0 1 0];
      elseif p==-1
       fbe=[1 0 0 0];
      end
     end
%      for ip1=1:npk
%       ip=ip+1;
      for ip1=ik1v
       ip=ip1+(ipom-1)*npk;
       Q=KK(ip1)*aloc1;
       P1=bes(ip1,jmu)^2;
       P2=bes(ip1,jmu-1)*bes(ip1,jmu+1);
       Kdia=.5*(P1-P2);
       Kdz(ip)=Kdia*fbe(3);
%       for ip2=ip1+1:npk
       for ip2=ip1+1:ik1v(end)
        Q1=KK(ip2)*aloc1;
        P1=bes(ip2,jmu)*bes(ip1,jmu-1);
        P2=bes(ip1,jmu)*bes(ip2,jmu-1);
        Kadia=1/(Q1^2-Q^2)*(Q*P1-Q1*P2);
        Kiiz(ip,ip+ip2-ip1)=Kadia*fbe(3);
       end
      end
    end %mu

    Ktoiizi(ndis,:,:)=(Kiiz+Kiiz'+diag(Kdz))*aloc^2;
   end % iztm


 end  %ndis

%  isu=prod(size(T));

  if igainshape==0
   isu=find(aytot==aiat(1));
  else
   isu=[1:lxi];
   if exist('xiT')
   isuT=[1:length(xiT)];
   end
  end
%  'isu Vett', keyboard

    si=size(Kiie);
    for ipui=1:lxi
     Ktoiie=reshape((Ktoiiei(ipui,:,:)),si);
     Ktoije=reshape((Ktoijei(ipui,:,:)),si);
     Ktoiim=reshape((Ktoiimi(ipui,:,:)),si);
     Ktoijm=reshape((Ktoijmi(ipui,:,:)),si);
     if isi==0
      KEEp=diag(ZEv.*KKv)*(Ktoiie);
      KEMp=diag(ZEv.*KKv)*(Ktoije);
      KMEp=diag(ZMv.*KKv)*(Ktoijm);
      KMMp=diag(ZMv.*KKv)*(Ktoiim);
     else
      diae=diag(sqrt(ZEv.*KKv));
      diam=diag(sqrt(ZMv.*KKv));
      diae(1,1)=diae(1,1)*fatde1;
      diam(1,1)=diam(1,1)*fatde1;
      KEEp=diae*(Ktoiie)*diae;
      KEMp=diae*(Ktoije)*diam;
      KMEp=diam*(Ktoijm)*diae;
      KMMp=diam*(Ktoiim)*diam;
     end
     Kosp(:,:,ipui)=[KEEp KEMp; KMEp KMMp];
     Kosm(:,:,ipui)=[KEEp KEMp; KMEp KMMp];
%     Kosm(:,:,ipui)=[KEEp -KEMp; -KMEp KMMp];
     if iztm==1
      Ktoiiz=reshape((Ktoiizi(ipui,:,:)),si);
      if isi==0
       KMMpz=diag(.5*ZMv.*KKv.*KKv./bev)*(Ktoiiz)*diag(KKv./bev);
      else
       diam=sqrt(ZMv.*KKv);
       diam(1)=diam(1)*fatde1;
       KMMpz=diag(.5*diam.*KKv./bev)*(Ktoiiz)*diag(KKv./bev.*diam);
      end
      Koszp(:,:,ipui)=[Kzer Kzer; Kzer KMMpz];
     else
      Koszp(:,:,ipui)=0;
     end
     Koszm=Koszp;
    end


%  if sha==shav.a(1) & igau==0
   ipui=find(ayi==aiat(1));
   if length(ipui)>0
    si2=2*si;
    KAp=reshape((Kosp(:,:,ipui)),si2);
    KAm=reshape((Kosm(:,:,ipui)),si2);
    if iztm==1
     KAzp=reshape((Koszp(:,:,ipui)),si2);
     KAzm=reshape((Koszp(:,:,ipui)),si2);
    end
   end


   if ifp>=3
    figure,
    surf(KAp), shading('interp'), view(0,90), colorbar, pausak
   end

   KApsub=KAp;
   KAmsub=KAm;
   if iztm==1
    KAzpsub=KAzp;
    KAzmsub=KAzm;
   end

%  end

Kospsub=Kosp;
Kosmsub=Kosm;
Koszpsub=Koszp;
Koszmsub=Koszm;


%'  prima di temp Vett ', keyboard

if igainshape==1
 if igau==4
  s=size(yiT);
  for iz=1:s(2)
   Ktoiie=0;
   Ktoije=0;
   Ktoiim=0;
   Ktoijm=0;
   Ktoiiz=0;
   indis=0;
   for ndis=isuT
    indis=indis+1;
     Kdum=reshape(Ktoiiei(ndis,:,:),si);
     Ktoiie=Ktoiie+Kdum*yiT(indis,iz);
     Kdum=reshape(Ktoijei(ndis,:,:),si);
     Ktoije=Ktoiie+Kdum*yiT(indis,iz);
     Kdum=reshape(Ktoiimi(ndis,:,:),si);
     Ktoiim=Ktoiim+Kdum*yiT(indis,iz);
     Kdum=reshape(Ktoijmi(ndis,:,:),si);
     Ktoijm=Ktoiie+Kdum*yiT(indis,iz);
   end  %ndis
     if isi==0
      KEEp=diag(ZEv.*KKv)*(Ktoiie);
      KEMp=diag(ZEv.*KKv)*(Ktoije);
      KMEp=diag(ZMv.*KKv)*(Ktoijm);
      KMMp=diag(ZMv.*KKv)*(Ktoiim);
     else
      KEEp=diae*(Ktoiie)*diae;
      KEMp=diae*(Ktoije)*diam;
      KMEp=diam*(Ktoijm)*diae;
      KMMp=diam*(Ktoiim)*diam;
     end
   KTemp(:,:,iz)=[KEEp KEMp; KMEp KMMp];
   indis=0;
   if iztm==1
    for ndis=isuT
     indis=indis+1;
     Kdum=reshape(Ktoiizi(ndis,:,:),si);
     Ktoiiz=Ktoiiz+Kdum*yiT(indis,iz);
    end  %ndis

      if isi==0
       KMMpz=diag(.5*ZMv.*KKv.*KKv./bev)*(Ktoiiz)*diag(KKv./bev);
      else
       diam=sqrt(ZMv.*KKv);
       diam(1)=diam(1)*fatde1;
       KMMpz=diag(.5*diam.*KKv./bev)*(Ktoiiz)*diag(KKv./bev.*diam);
      end
    KTempz(:,:,iz)=[Kzer Kzer; Kzer KMMpz];
   else
    KTempz(:,:,iz)=[Kzer Kzer; Kzer Kzer];
   end
  end
 end  %igau==4

%' fine KTemp vett ', keyboard
end

' fine matv_ci1 ', keyboard



%clear KAp KAm Kosp Kosm
%
%if iexK==1
% KAp=KAps;
% KAm=KAms;
%end

