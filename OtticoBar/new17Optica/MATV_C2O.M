
fime=find(adis<0);

if length(fime)==1
[adisd,iso]=sort(ailu_mean*kcav0);
fid=find(diff([0; adisd])>0);
adis=adisd(fid);
clear fmul
for ka=1:length(adis)
 fi=find(adis(ka)==adisd);
 fmul(ka,1)=length(fi);
end

lxivet=1:length(adis);
dx_mean=dilu_mean(iso(fid)).*fmul;
fattme=dx_mean/sum(dx_mean);

' entro ci vett mediato ', keyboard

ZEKK=ZEv(Pusas0).*KKv(Pusas0);
ZMKK=ZMv(Pusas0).*KKv(Pusas0);
 lPU=length(Pus0);
% lPU=Pusas0(end);
 Kset=zeros(lPU,lPU);
 Vset=zeros(1,lPU);

 mbv1=[-2+nubesi:nubesu+2];
% bes=besselj(mbv1,KK*a);

 nup1=numodi;

% Kiie=zeros(npk*nup1,npk*nup1);
% Kzer=Kiie;
% Kije=zeros(npk*nup1,npk*nup1);
% Kde=zeros(1,npk*nup1);
% Kd1e=zeros(1,npk*nup1);
% Kiim=zeros(npk*nup1,npk*nup1);
% Kijm=zeros(npk*nup1,npk*nup1);
% Kdm=zeros(1,npk*nup1);
% Kd1m=zeros(1,npk*nup1);
% Kdz=zeros(1,npk*nup1);
% Kiiz=zeros(npk*nup1,npk*nup1);

%%%%%%%%%%%%%%%%%%%%%%%
%'qui', keyboard

 clear  Ktoiiei Ktoijei Ktoiimi Ktoijmi
 clear  Ktoiizi Kosp Koszp Kosm Koszm
' ver vet', keyboard


 for ndis=lxivet

 Kiie=Kset;
 Kzer=Kiie;
 Kije=Kset;
 Kde=Vset;
 Kd1e=Vset;
 Kiim=Kset;
 Kijm=Kset;
 Kdm=Vset;
 Kd1m=Vset;

 Kdz=Vset;
 Kiiz=Kset;

  aloc1=adis(ndis);
%  aloc2=aloc1/kcav0;
  alocver=aloc1/kcav0;
  aloc=aloc1;
  bes=besselj(mbv1,KK*aloc1);

  ip=0;


 for imu=pimu+1
% for imu=pimu
   for iarm=1:farm
    ipom=(ceil((imu-2)/pasnu));
    ipom1=(ceil((farm*(imu-2)+iarm-1)/pasnu));
%    ipom=(ceil((imu)/pasnu));
%    ipom1=(ceil((farm*(imu)+iarm-1)/pasnu));
%    [imu ipom ipom1]
%    keyboard
    ik1v=Pusas0([ldap(ipom1)+1:ldap(ipom1+1)])-(ipom-1)*npk;


%    if exist('puk_ret') & alocver>ra0
%     fist=find(diff(Pusas0)>1)+1;
%     fist=fist(1);
%     ip=fist-1;
%     ik1v=puk_ret-fist;
%    end


    if alocver>ra0
     if exist('puk_ret')
      if imu==pimu(1)+1
       fist=find(diff(Pusas0)>1)+1;
       fist=fist(1);
       ip=fist-1;
       ik1v=puk_ret-fist;
      end
     else
      ik1v=[];
     end

    end




%   ipom=(ceil((imu-1)/pasnu));
%   ik1v=Pusas0([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;
 mu=mbv1(imu);
 jmu=imu;
% mu=mbv1(imu+1);
% jmu=imu+1;
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
%     ip=ip1+(ipom-1)*npk;
%    for ip1=1:npk
     ip=ip+1;
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
     fip2=find(ik1v==ip1);
     ipc=ip;
     for ip2=ik1v(fip2+1:end)
%     for ip2=ip1+1:ik1v(end)
      ipc=ipc+1;
      Q1=KK(ip2)*aloc1;
      P1=bes(ip2,jmu-1)*bes(ip1,jmu-2)+bes(ip2,jmu+1)*bes(ip1,jmu);
      P2=bes(ip1,jmu-1)*bes(ip2,jmu-2)+bes(ip1,jmu+1)*bes(ip2,jmu);
      Kadia=1/(Q1^2-Q^2)*(Q*P1-Q1*P2);
%'      [imu ip1 ip2 ip], pausak '
%      [imu ip1 ip2 ip], pausak
%'      [ip  ip+ip2-ip1 ip1 ip2], pausak '
%      [ip  ip+ip2-ip1 ip1 ip2], pausak
%      ipc=ip+ip2-ip1;
      Kiie(ip,ipc)=Kadia*fbe(1);
      Kiim(ip,ipc)=Kadia*fbe(3);
      P1=-bes(ip2,jmu-1)*bes(ip1,jmu-2)+bes(ip2,jmu+1)*bes(ip1,jmu);
      P2=-bes(ip1,jmu-1)*bes(ip2,jmu-2)+bes(ip1,jmu+1)*bes(ip2,jmu);
      Kadia=1/(Q1^2-Q^2)*(Q*P1-Q1*P2);
      Kije(ip,ipc)=Kadia*fbe(2);
      Kijm(ip,ipc)=Kadia*fbe(4);
     end
    end
   end % farm
 end %mu

    Ktoiiei(ndis,:,:)=(Kiie+Kiie'+diag(Kde))*aloc^2/4;
    Ktoijei(ndis,:,:)=-(Kije+Kije'+diag(Kd1e))*aloc^2/4;
    Ktoiimi(ndis,:,:)=(Kiim+Kiim'+diag(Kdm))*aloc^2/4;
    Ktoijmi(ndis,:,:)=-(Kijm+Kijm'+diag(Kd1m))*aloc^2/4;

%' qui fin ', keyboard

   if iztm==1
    ip=0;
    for imu=pimu
 %   ipom=(ceil(imu/pasnu));
 %   ik1v=Pusas0([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;
   for iarm=1:farm
    ipom=(ceil((imu-1)/pasnu));
    ipom1=(ceil((farm*(imu-1)+iarm-1)/pasnu));
%    ipom=(ceil((imu-1)/pasnu));
%    ipom1=(ceil((farm*(imu-1)+iarm-1)/pasnu));
    ik1v=Pusas0([ldap(ipom1)+1:ldap(ipom1+1)])-(ipom-1)*npk;

    if alocver>ra0
     if exist('puk_ret')
      if imu==pimu(1)
       fist=find(diff(Pusas0)>1)+1;
       fist=fist(1);
       ip=fist-1;
       ik1v=puk_ret-fist;
      end
     else
      ik1v=[];
     end

    end


%    if exist('puk_ret') & alocver>ra0
%     ik1v=puk_ret;
%     if imu==pimu(1)+1
%      ip=ik1v(1)-1;
%     end
%    end
%    if imu==pimu(1) & alocver>ra0
%     fist=find(diff(Pusas0)>1)+1;
%     pu=fist:length(Pusas0);
%     ik1v=Pusas0(pu)-(ipom-1)*npk;
%     ip=fist-1;
%    end


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
       ip=ip+1;
%       ip=ip1+(ipom-1)*npk;
       Q=KK(ip1)*aloc1;
       P1=bes(ip1,jmu)^2;
       P2=bes(ip1,jmu-1)*bes(ip1,jmu+1);
       Kdia=.5*(P1-P2);
       Kdz(ip)=Kdia*fbe(3);
%       for ip2=ip1+1:npk
     fip2=find(ik1v==ip1);
     ipc=ip;
     for ip2=ik1v(fip2+1:end)
      ipc=ipc+1;
%       for ip2=ip1+1:ik1v(end)
        Q1=KK(ip2)*aloc1;
        P1=bes(ip2,jmu)*bes(ip1,jmu-1);
        P2=bes(ip1,jmu)*bes(ip2,jmu-1);
        Kadia=1/(Q1^2-Q^2)*(Q*P1-Q1*P2);
%       ipc=ip+ip2-ip1;
        Kiiz(ip,ipc)=Kadia*fbe(3);
       end
      end
     end %iarm
    end %mu

    Ktoiizi(ndis,:,:)=(Kiiz+Kiiz'+diag(Kdz))*aloc^2;
   end % iztm

%' end size ', keyboard

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
      Madd=0;
      Maddz=0;

    si=size(Kiie);
%    for ipui=1:lxi
%     alocver=adis(ipui)/kcav0;
    siK=size(Ktoiiei);
    for ipui=1:siK(1)
%     if alocver>ra0
%      Madd=Madd0;
%      Maddz=Maddz0;
%     else
%      Madd=0;
%      Maddz=0;
%     end

     Ktoiie=reshape((Ktoiiei(ipui,:,:)),si);
     Ktoije=reshape((Ktoijei(ipui,:,:)),si);
     Ktoiim=reshape((Ktoiimi(ipui,:,:)),si);
     Ktoijm=reshape((Ktoijmi(ipui,:,:)),si);
     if isi==0
      KEEp=diag(ZEKK)*(Ktoiie);
      KEMp=diag(ZEKK)*(Ktoije);
      KMEp=diag(ZMKK)*(Ktoijm);
      KMMp=diag(ZMKK)*(Ktoiim);
     else
      diae=diag(sqrt(ZEKK));
      diam=diag(sqrt(ZMKK));
      diae(1,1)=diae(1,1)*fatde1;
      diam(1,1)=diam(1,1)*fatde1;
      KEEp=diae*(Ktoiie)*diae;
      KEMp=diae*(Ktoije)*diam;
      KMEp=diam*(Ktoijm)*diae;
      KMMp=diam*(Ktoiim)*diam;
     end

     Kospdu=[KEEp segem*KEMp; segem*KMEp KMMp]+Madd;
     Kosp(:,:,ipui)=Kospdu;
     Kosmdu=[KEEp segem*KEMp; segem*KMEp KMMp]+Madd;
     Kosm(:,:,ipui)=Kosmdu;
%     Kosm(:,:,ipui)=[KEEp -KEMp; -KMEp KMMp];
     if iztm==1
      Ktoiiz=reshape((Ktoiizi(ipui,:,:)),si);
       Bva=KKv./bev;
       Bva=Bva(Pusas0);
      if isi==0
       KMMpz=diag(.5*ZMKK.*Bva)*(Ktoiiz)*diag(Bva);
      else
       diam=sqrt(ZMKK);
       diam(1)=diam(1)*fatde1;
       KMMpz=diag(.5*diam.*Bva)*(Ktoiiz)*diag(Bva.*diam);
      end
      Koszpdu=[Kzer Kzer; Kzer KMMpz]+Maddz;
      Koszp(:,:,ipui)=Koszpdu;
      Koszm(:,:,ipui)=Koszpdu;
     else
      Koszp(:,:,ipui)=0;
     end
     Koszm=Koszp;
%     pausak
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
  if exist('KAp')
   KApsub=KAp;
   KAmsub=KAm;
   if iztm==1
    KAzpsub=KAzp;
    KAzmsub=KAzm;
   end
  end



Kosp_me=mean(Kosp,3);
Kosm_me=mean(Kosm,3);
Koszp_me=mean(Koszp,3);
Koszm_me=mean(Koszm,3);

lxivet=lxivet_sa;
adis=adis_sa;
clear KAp KAm Kosp Kosm Koszp Koszm

' fine matve_cim', keyboard


%clear KAp KAm Kosp Kosm
%
%if iexK==1
% KAp=KAps;
% KAm=KAms;
%end
end
