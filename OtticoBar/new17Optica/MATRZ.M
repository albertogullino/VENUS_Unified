%' entro du ', keyboard
%iexK=0;
%if exist('KAp')
% KAps=KAp;
% KAms=KAm;
% iexK=1;
%end

 mbv1du=[-1 0];
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

 clear  Ktoiizi Koszp Koszm

 for ndis=lxivet
  aloc1=adis(ndis);
%  aloc2=aloc1/kcav0;
  aloc=aloc1;
  bes=besselj(mbv1du,KK*aloc1);

   ifidis=find(ifield<0);
   afiel=aitot(ifidis);
%  if istrumix==1 & igainshape==0
%   if (ifii(ndis)<=0 | length(find(afiel==aloc/kcav0))>0)
  if istrumix==1 & igainshape==0 & DOE==0
   ifidis=find(ifield<0);
   afiel=radii.a(ifidis);
   Prele=10;
   if exist('Rel')
    if Rel>0
     Prele=abs(aloc/kcav0-Rel);
    end
   end
   if shaC==7
    ReS=PAle{1}.Rel*kcav;
    ReSe=PAle{1}.Rm_rel*kcav;
    condi=(ifii(ndis)<0 | aloc==ReS | aloc==ReSe );
   else
    condi=(ifii(ndis)<0 |  (ndis==1 | ndis==length(adis)) | (Prele<.01) | lxi<5);
   end
%   shaC
%   condi
%   keyboard
%   if iplan==1
    condi=0;
%   end
   if condi==1
%   if (ifii(ndis)<0 |  (ndis==1 | ndis==length(adis)) | (Prele<.01) | lxi<5)
   sht=sha_ish(fis(fish));
%   'afiel VE', keyboard
   if mean(sht(lxivet))==7
    sha_plo=7;
   end

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

   if iztm==1
    ip=0;
    for imu=2
    mu=0;
    jmu=2;
      for ip1=1:npk
       ip=ip+1;
       Q=KK(ip1)*aloc1;
       P1=bes(ip1,jmu)^2;
       P2=-bes(ip1,jmu-1)^2;
       Kdia=.5*(P1-P2);
       Kdz(ip)=Kdia;
       for ip2=ip1+1:npk
        Q1=KK(ip2)*aloc1;
        P1=bes(ip2,jmu)*bes(ip1,jmu-1);
        P2=bes(ip1,jmu)*bes(ip2,jmu-1);
        Kadia=1/(Q1^2-Q^2)*(Q*P1-Q1*P2);
        Kiiz(ip,ip+ip2-ip1)=Kadia;
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
     if iztm==1
      Ktoiiz=reshape((Ktoiizi(ipui,:,:)),si);
      if isi==0
       KMMpz=diag(.5*ZMv.*KKv.*KKv./bev)*(Ktoiiz)*diag(KKv./bev);
      else
       diam=sqrt(ZMv.*KKv);
       diam(1)=diam(1)*fatde1;
       KMMpz=diag(.5*diam.*KKv./bev)*(Ktoiiz)*diag(KKv./bev.*diam);
      end
      Koszp(:,:,ipui)=[KMMpz];
     else
      Koszp(:,:,ipui)=0;
     end
     Koszm=Koszp;
    end


%  if sha==shav.a(1) & igau==0
   ipui=find(ayi==aiat(1));
   if length(ipui)>0
    si2=si;
    if iztm==1
     KAzp=reshape((Koszp(:,:,ipui)),si2);
     KAzm=reshape((Koszp(:,:,ipui)),si2);
    end
   end


   if ifp>=3
    figure,
    surf(KAzp), shading('interp'), view(0,90), colorbar, pausak
   end

   if exist('KAp')
    if iztm==1
     KAzpsub=KAzp;
     KAzmsub=KAzm;
    end
   end

%  end

Koszpsub=Koszp;
Koszmsub=Koszm;
%   ' Kto scal fint',pausak
