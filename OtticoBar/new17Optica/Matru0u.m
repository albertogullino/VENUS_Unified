fatde1=1;
if ~exist('yi')
 yi=ones(lxi,1);
end

 mbv1=[-1 0];

 nup1=1;
 Kiie=zeros(npk*nup1,npk*nup1);
 Kzer=Kiie;

 Kdz=zeros(1,npk*nup1);
 Kiiz=zeros(npk*nup1,npk*nup1);
%%%%%%%%%%%%%%%%%%%%%%%

 clear  Ktoiizi Koszm Koszp  Kdz

 for ndis=1:lxi
  aloc1=adis(ndis);
%  aloc2=aloc1/kcav0;
  aloc=aloc1;
  bes=besselj(mbv1,KK*aloc1);


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

    Ktoiizi(:,:,ndis)=(Kiiz+Kiiz'+diag(Kdz))*yi(ndis)*aloc^2;
   end % iztm

 end  %ndis

 ' ver', keyboard

  si=size(Kiie);
  if nirid~=0
   isu=[1:lxi-nirid];
   isui=lxi-nirid;
  else
   isu=[1:lxi];
  end
  if lxi>1
   if length(isu)>1
    if iztm==1
     Ktoiiz=reshape(sum(Ktoiizi(isu,:,:)),si);
    end
   else
    if iztm==1
     Ktoiiz=reshape((Ktoiizi(isu,:,:)),si);
    end
   end
   if iztm==1
    if isi==0
     KMMpz=diag(0.5*ZMv.*KKv.*KKv./bev)*(Ktoiiz)*diag(KKv./bev);
    else
     diam=sqrt(ZMv.*KKv);
     diam(1)=diam(1)*fatde1;
     KMMpz=diag(0.5*diam.*KKv./bev)*(Ktoiiz)*diag(KKv./bev.*diam);
    end
    KAzp=[KMMpz];
   else
    KAzp=0;
   end

   if nirid~=0
%    for inir=1:nirid
%     ipui=isui+inir;
    for inir=1:length(ipuos)
     ipui=ipuos(inir);
     disp(' matru0u '), keyboard
     if iztm==1
      Ktoiiz=reshape((Ktoiizi(ipui,:,:)),si);
      if isi==0
       KMMpz=diag(.5*ZMv.*KKv.*KKv./bev)*(Ktoiiz)*diag(KKv./bev);
      else
       diam=sqrt(ZMv.*KKv);
       diam(1)=diam(1)*fatde1;
       KMMpz=diag(.5*diam.*KKv./bev)*(Ktoiiz)*diag(KKv./bev.*diam);
      end
      Koszp(:,:,inir)=[KMMpz];
     else
      Koszp(:,:,inir)=0;
     end

    end
   end
  else
    inir=1;
    ipui=ipuos(inir);
     if iztm==1
      Ktoiiz=reshape((Ktoiizi(ipui,:,:)),si);
      if isi==0
       KMMpz=diag(.5*ZMv.*KKv.*KKv./bev)*(Ktoiiz)*diag(KKv./bev);
      else
       diam=sqrt(ZMv.*KKv);
       diam(1)=diam(1)*fatde1;
       KMMpz=diag(.5*diam.*KKv./bev)*(Ktoiiz)*diag(KKv./bev.*diam);
      end
      KAzp=[KMMpz];
     else
      KAzp=0;
      KAzm=0;
     end

     if iztm==1
      Koszp(:,:,inir)=KAzp;
     else
      Koszp(:,:,inir)=0;
     end

  end

  KAzm=KAzp;
  Koszm=Koszp;

  if ifp>=3
   figure,
   surf(KAzm), shading('interp'), view(0,90), colorbar, pausak
  end

si=size(Kiie);
