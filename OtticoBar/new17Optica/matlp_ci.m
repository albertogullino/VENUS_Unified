%' entro ci ', keyboard


iexK=0;
if exist('KAp')
 KAps=KAp;
 KAms=KAm;
 iexK=1;
end

 mbv1=[-1+nubesi:nubesu+1];

 nup1=numodi;
 Kiie=zeros(npk*nup1,npk*nup1);
 Kde=zeros(1,npk*nup1);
 clear  Ktoiiei Kovri

 for ndis=lxivet
  aloc=adis(ndis);
  if istrumix==1
   if ifii(ndis)<0
    forma
    scu=size(xcu);
    pug=1:scu(1);
    if exist('Cug.x')==1
     sC2=size(Cug.x);
     sC2=sC(2);
    else
     sC2=0;
    end
    puc=[1:scu(2)]+sC2;
    Cug.x(pug,puc)=xcu;
    Cug.y(pug,puc)=ycu;
    Cug.z(pug,puc)=zcu;
   end

%   if ifii(ndis)<0
%    forma
%    icosha=icosha+1;
%    Cug.x(:,icosha)=xcu';
%    Cug.y(:,icosha)=ycu';
%    Cug.z(:,icosha)=zcu';
%   end
  end
  bes=besselj(mbv1,KK*aloc);
  ip=0;
%  for imu=2
%  mu=mbv1(imu);
%  jmu=imu;
  for imu=pimu
  mu=mbv(imu);
  jmu=imu+1;
    for ip1=1:npk
     ip=ip+1;
     Q=KK(ip1)*aloc;
     P1=bes(ip1,jmu)^2;
     P2=bes(ip1,jmu-1)*bes(ip1,jmu+1);
     Kdia=.5*(P1-P2);
     Kde(ip)=Kdia;
     for ip2=ip1+1:npk
      Q1=KK(ip2)*aloc;
      P1=bes(ip2,jmu)*bes(ip1,jmu-1);
      P2=bes(ip1,jmu)*bes(ip2,jmu-1);
      Kadia=1/(Q1^2-Q^2)*(Q*P1-Q1*P2);
      Kiie(ip,ip+ip2-ip1)=Kadia;
     end
    end
  end %mu
  Ktoiiei(ndis,:,:)=(Kiie+Kiie'+diag(Kde))*aloc^2/2;
  Kti(ndis,:,:)=(Kiie+Kiie'+diag(Kde))*aloc^2/2;
  Kovri(ndis,:,:)=(Kiie+Kiie'+diag(Kde));
    ' scal', pausak
 end  %ndis

  si=size(Kiie);
  if lxi>1
   if length(isu)>1
    Ktoiie=reshape(sum(Ktoiiei(isu,:,:)),si);
    if isi==0
     dia=diag(ZEv.*KKv);
     dia1=1;
     KA=dia*Ktoiie;
    else
     dia=diag(sqrt(ZEv.*KKv));
     dia1=dia;
     KA=dia*Ktoiie*dia1;
    end
   else
    Ktoiie=reshape((Ktoiiei(isu,:,:)),si);
    if isi==0
     dia=diag(ZEv.*KKv);
     dia1=1;
     KA=dia*Ktoiie;
    else
     dia=diag(sqrt(ZEv.*KKv));
     dia1=dia;
     KA=dia*Ktoiie*dia1;
    end
   end
   if nirid~=0
    for inir=1:length(ipuos)
     ipui=ipuos(inir);
     if isi==0
      Kos(:,:,inir)=diag(ZEv.*KKv)*reshape(Ktoiiei(ipui,:,:),si);
     else
      dia=diag(sqrt(ZEv.*KKv));
      Kos(:,:,inir)=dia*reshape(Ktoiiei(ipui,:,:),si)*dia;
     end
     KoN(:,:,inir)=reshape(Kovri(ipui,:,:),si);
    end
   end
  else
   Ktoiie=reshape(Ktoiiei,si);
   if isi==0
    KA=diag(ZEv.*KKv)*Ktoiie;
   else
    dia=diag(sqrt(ZEv.*KKv));
    KA=dia*Ktoiie*dia;
   end
   if nirid~=0
    Kos=KA;
    KoN=reshape(Kovri(1,:,:),si);
   end
  end


  if ifp>=1
   figure,
   surf(KA), shading('interp'), view(0,90), colorbar, pausak
  end
  Znor=ZEv;

 KAp=KA;
 KAm=KA;
 if exist('Kos')
  Kosp=Kos;
  Kosm=Kos;
 else
  Kosp=KA;
  Kosm=KA;
 end

KApsub=KAp;
KAmsub=KAm;
Kospsub=Kosp;
Kosmsub=Kosm;

clear KAp KAm Kosp Kosm

if iexK==1
 KAp=KAps;
 KAm=KAms;
end

if istrumix==0
 if igau==4
  s=size(yiT);
  for iz=1:s(2)
   Ktoiie=0;
   indis=0;
   for ndis=isuT
    indis=indis+1;
    Kdum=reshape(Kti(ndis,:,:),si);
    Ktoiie=Ktoiie+Kdum*yiT(indis,iz);
   end  %ndis
   KTemp(:,:,iz)=dia*Ktoiie*dia1;
  end
 end  %igau==4
end




