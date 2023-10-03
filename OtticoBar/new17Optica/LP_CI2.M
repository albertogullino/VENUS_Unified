ZEKK=ZEv(Pusas).*KKv(Pusas);
iexK=0;
%'ci'
%keyboard
if exist('KAp')
 KAps=KAp;
 KAms=KAm;
 iexK=1;
end

if istrumix==1
 nirid=length(lxivet);
 ipuos=lxivet;
end

 mbv1=[-1+nubesi:nubesu+1];

 nup1=numodi;
  lPU=length(Pus);
  Kset=zeros(lPU,lPU);
% Kiie=zeros(npk*nup1,npk*nup1);
% Kde=zeros(1,npk*nup1);
 Kiie=Kset;
 Kde=zeros(1,lPU);

 clear  Ktoiiei Kovri Kti

 for ndis=lxivet
  aloc=adis(ndis);

  if istrumix==1 & igainshape==0
   ifidis=find(ifield<0);
   afiel=radii.a(ifidis);
%   'afiel LP', keyboard
   if (ifii(ndis)<=0 | length(find(afiel==aloc/kcav0))>0)
%   if ifii(ndis)<0
%    ' qui prima di forma ', keyboard
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

   for iarm=1:farm
    ipom=(ceil(imu/pasnu));
    ipom1=(ceil((farm*imu+iarm-1)/pasnu));
    ik1v=Pusas([ldap(ipom1)+1:ldap(ipom1+1)])-(ipom-1)*npk;
    [imu iarm]
%   'ik1v', keyboard
    mu=mbv(imu);
    jmu=imu+1;
    for ip1=ik1v
     ip=ip+1;
     Q=KK(ip1)*aloc;
     P1=bes(ip1,jmu)^2;
     P2=bes(ip1,jmu-1)*bes(ip1,jmu+1);
     Kdia=.5*(P1-P2);
     Kde(ip)=Kdia;
%     for ip2=ip1+1:npk
     ipc=ip;
     fip2=find(ik1v==ip1);
     for ip2=ik1v(fip2+1:end)
%     for ip2=ip1+1:ik1v(end)
      ipc=ipc+1;
      Q1=KK(ip2)*aloc;
      P1=bes(ip2,jmu)*bes(ip1,jmu-1);
      P2=bes(ip1,jmu)*bes(ip2,jmu-1);
      Kadia=1/(Q1^2-Q^2)*(Q*P1-Q1*P2);
%      Kiie(ip,ip+ip2-ip1)=Kadia;
      Kiie(ip,ipc)=Kadia;
     end
    end
   end % farm

  end %mu
  Ktoiiei(ndis,:,:)=(Kiie+Kiie'+diag(Kde))*aloc^2/2;    % conta il Q*Del Q
  Kti(ndis,:,:)=(Kiie+Kiie'+diag(Kde))*aloc^2/2;
  Kovri(ndis,:,:)=(Kiie+Kiie'+diag(Kde));
 end  %ndis

%' qui K ', keyboard
  si=size(Kiie);
%  isu=prod(size(T));

  if igainshape==0
   isu=find(aytot(:,1)==aiat(1));
  else
   isu=[1:lxi];
   if exist('xiT')
   isuT=[1:length(xiT)];
   end
  end
%  'isu LP', keyboard

  if lxi>1
   if length(isu)>1

%    Kdu=0;
%    for nsu=isu
%     Kdu=Kdu+Ktoiiei(nsu,:,:)*yiN(nsu);
%    end
%    Ktoiie=reshape(Kdu(1,:,:),si);
    if isi==0
     dia=diag(ZEKK);
     dia1=1;
%     KA=dia*Ktoiie;
    else
     dia=diag(sqrt(ZEKK));
     dia1=dia;
%     KA=dia*Ktoiie*dia1;
    end
%
%%    if exist('yiN_ag') & ianti_gui==1
%     Kdu=0;
%     for nsu=isu
%      Kdu=Kdu+Ktoiiei(nsu,:,:)*yiN_ag(nsu);
%     end
%     Ktoiie=reshape(Kdu(1,:,:),si);
%     if isi==0
%      dia=diag(ZEv.*KKv);
%      dia1=1;
%      KA_ag=dia*Ktoiie;
%     else
%      dia=diag(sqrt(ZEv.*KKv));
%      dia1=dia;
%      KA_ag=dia*Ktoiie*dia1;
%     end
%    end


   elseif length(isu)==1
    Ktoiie=reshape((Ktoiiei(isu,:,:)),si);
    if isi==0
     dia=diag(ZEKK);
     dia1=1;
     KA=dia*Ktoiie;
    else
     dia=diag(sqrt(ZEKK));
     dia1=dia;
     KA=dia*Ktoiie*dia1;
    end
   end

   if nirid~=0
    for inir=1:length(ipuos)
     ipui=ipuos(inir);
     if isi==0
      Kos(:,:,inir)=diag(ZEKK)*reshape(Ktoiiei(ipui,:,:),si);
     else
      dia=diag(sqrt(ZEKK));
      Kos(:,:,inir)=dia*reshape(Ktoiiei(ipui,:,:),si)*dia;
     end
     KoN(:,:,inir)=reshape(Kovri(ipui,:,:),si);
    end
   end

  else
   Ktoiie=reshape(Ktoiiei,si);
   if isi==0
    KA=diag(ZEKK)*Ktoiie;
   else
    dia=diag(sqrt(ZEKK));
    KA=dia*Ktoiie*dia;
   end
%  'qi sper', keyboard

   if nirid~=0
    Kos=KA;
    KoN=reshape(Kovri(1,:,:),si);
   end
  end


  if ifp>=1 & exist('KA')
   figure,
   surf(KA), shading('interp'), view(0,90), colorbar, pausak
  end
  Znor=ZEv(Pusas);

  if exist('KA')
   KAp=KA;
   KAm=KA;

   KApsub=KAp;
   KAmsub=KAm;
  end

   if exist('Kos')
    Kosp=Kos;
    Kosm=Kos;
   else
%    Kosp=KA;
%    Kosm=KA;
    Kosp=0;
    Kosm=0;
   end

   Kospsub=Kosp;
   Kosmsub=Kosm;
%  'dopo isu LP: fine K', keyboard

clear KAp KAm Kosp Kosm

if iexK==1
 KAp=KAps;
 KAm=KAms;
end

%' prima di Temp LP ', keyboard

if igainshape==1
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

%' prima di DOE LP ', keyboard
if DOE==1
   Ktoiie=0;
   for ndis=lxivet
    Ktoiie=Ktoiie+Wdis(ndis)*reshape(Kti(ndis,:,:),si);
   end  %ndis
   if isi==0
    Kosp=diag(ZEKK)*Ktoiie;
   else
    dia=diag(sqrt(ZEKK));
    Kosp=dia*Ktoiie*dia;
   end
   Kosm=Kosp;
   if ifp==-10
   map(Kosp), pausak
   Wd=[Wdis; Wdis];

   Wd=reshape(Wd,prod(size(Wd)),1);
   Ad=adis/kcav0;
   Ad=[Ad; Ad];
   Ad=reshape(Ad,prod(size(Wd)),1);
   Ad=[0; Ad; Ad(end)+5];
%   Wd=[Wd; Wd(end); Wd(end)];
   Wd=[Wd; -1; -1];
   figure, plot(Ad,(Wd+1)/2,'r')
   ' dopo di DOE LP ', keyboard
   end
end


%' fine ci ', keyboard

