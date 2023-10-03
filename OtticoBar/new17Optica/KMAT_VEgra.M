igainshape=1;  %non scrive reticolo per disegno
fatde1=1;
imu0=0;
if ifp~=-4
disp(' ')
disp(' ')
end

%ZEv=ZEv*1.2;
%'%ZMv=ZMv*1.2; ', keyboard

ZEKK=ZEv(Pusas0).*KKv(Pusas0);
ZMKK=ZMv(Pusas0).*KKv(Pusas0);

clear A B AB

alon=1;

pimu=2:pasnu:length(mbv)-1;
%'entro ve2', keyboard
meun=1;

 mbv1=[-2+nubesi:nubesu+2];
   if iany==2 | ianys==2
    ianysd=1;
   else
    ianysd=0;
   end

%'kmatve', pausak



 sha_grat_lens

raver=rv/kcav;
for ktro=1:length(rad_any)
 [du,imi]=min(abs(raver-rad_any(ktro)));
 po_radi(ktro)=imi;
end

adis=raver(po_radi);
pura=1:length(adis);


 if ifalso<=0
  lKK=length(KK);
  lrv=length(rv);
  z=KK*rv;

  clear besr besz

  for iz=1:length(mbv)
   nu=mbv(iz);
   bes0=real(besselj(nu,z));
   besr(:,:,iz)=bes0;
  end

%  for iz=1:length(mbv)
%   nu=mbv(iz);
%   bes0=besselj(mbv,-1)
%  end

la=length(adis);
npk1=length(Pusas0);
npk1=ldap(2);
 Fi1=zeros(nubes,nubes,npk1,npk1,la);
 Fi2=Fi1;
 Fi3=Fi1;
 Fi4=Fi1;
 Fi5=Fi1;
 Fi6=Fi1;
%'qui ', keyboard

 rvp=[rv rv(lrv)+diff(rv(lrv-1:lrv))];
 R1=rvp/alon;
 Rp=R1(1:lrv);

 if igint==1
  RdR=sgimp.*wi/alon.*Rp;
 else
  dR=diff(R1);
  RdR=sgimp.*Rp.*dR;
 end

    if ifp>1
      figure
    end

 ifig=0;
     ik1c=0; 
     ik2c=0; 
imuc=0;
  for imu=pimu
  imuc=imuc+1;

 %  ipom=(ceil(imu/pasnu));
 %  ik1v=Pusas0([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;
   for iarm=1:farm
    ipom=(ceil((imu-1)/pasnu));
    ipom1=(ceil((farm*(imu-1)+iarm-1)/pasnu));
    ik1v=Pusas0([ldap(ipom1)+1:ldap(ipom1+1)])-(ipom-1)*npk;

  jmu=imu-1;
  mu=mbv(imu);
   fimu=find(pimu>=imu);
inuc=0;

   for inu=pimu
   inuc=inuc+1;

%    ipom=(ceil(inu/pasnu));
%    ik2v=Pusas0([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;
    for iarn=1:farm
     ipom=(ceil((inu-1)/pasnu));
     ipom1=(ceil((farm*(inu-1)+iarn-1)/pasnu));
     ik2v=Pusas0([ldap(ipom1)+1:ldap(ipom1+1)])-(ipom-1)*npk;
     ik2c=0; 

   jnu=inu-1;
   nu=mbv(inu);
%    [jmu jnu], pausak
  if inu>=imu

   if (nu+mu)/2-fix((nu+mu)/2)==0
%    for ik1=1:npk
%     for ik2=1:npk
     ik1c=0; 
     for ik1=ik1v
     ik1c=ik1c+1;
     ik2c=0; 
      for ik2=ik2v
      ik2c=ik2c+1;
%      pausak
      F1=.5*besr(ik2,:,inu-1).*besr(ik1,:,imu-1).*RdR.*A(:,jmu,jnu)';
      Fid=cumsum(F1);
      Fi1(imuc,inuc,ik1c,ik2c,pura)=Fid(po_radi);
      F2=.5*besr(ik2,:,inu+1).*besr(ik1,:,imu+1).*RdR.*A(:,jmu,jnu)';
      Fid=cumsum(F2);
      Fi2(imuc,inuc,ik1c,ik2c,pura)=Fid(po_radi);
      F3=.5*besr(ik2,:,inu-1).*besr(ik1,:,imu+1).*RdR.*B(:,jmu,jnu)';
      Fid=cumsum(F3);
      Fi3(imuc,inuc,ik1c,ik2c,pura)=Fid(po_radi);
      F4=.5*besr(ik2,:,inu+1).*besr(ik1,:,imu-1).*RdR.*B(:,jmu,jnu)';
      Fid=cumsum(F4);
      Fi4(imuc,inuc,ik1c,ik2c,pura)=Fid(po_radi);
      if iztm==1
       F5d=besr(ik2,:,inu).*besr(ik1,:,imu).*RdR;
       F5=F5d.*A(:,jmu,jnu)';
       Fid=cumsum(F5);
       Fi5(imuc,inuc,ik1c,ik2c,pura)=Fid(po_radi);
       F6=F5d.*B(:,jmu,jnu)';
       Fid=cumsum(F6);
       Fi6(imuc,inuc,ik1c,ik2c,pura)=Fid(po_radi);
      end
     end %ik2
    end  %ik1
    if ifp>=2
     if iztm==1
 % tutti   plot(r,F1,r,F2,r,F3,r,F4,r,F5,r,F6)
       plot(r,F1,r,F2,'.')
     else
 % tutti   plot(r,F1,r,F2,r,F3,r,F4)
      plot(r,F1,r,F2,'.')
     end
     pausak
    end
    if ifig==1
     figure
     surf(reshape(Fi1(jmu,jnu,:,:,1),npk1,npk1)), shading('interp'), view(0,90), colorbar, pausak
     figure
     surf(reshape(Fi2(jmu,jnu,:,:,1),npk1,npk1))
     shading('interp'), view(0,90), colorbar, pausak
     figure
     surf(reshape(Fi3(jmu,jnu,:,:,1),npk1,npk1))
     shading('interp'), view(0,90), colorbar, pausak
     figure
     surf(reshape(Fi4(jmu,jnu,:,:,1),npk1,npk1))
     shading('interp'), view(0,90), colorbar, pausak
     close all
    end
   end  %if

  else
    for krad=pura
     Fi1(imuc,inuc,:,:,krad)=reshape(Fi1(inuc,imuc,:,:,krad),npk1,npk1).';
     Fi2(imuc,inuc,:,:,krad)=reshape(Fi2(inuc,imuc,:,:,krad),npk1,npk1).';
     Fi3(imuc,inuc,:,:,krad)=reshape(Fi4(inuc,imuc,:,:,krad),npk1,npk1).';
     Fi4(imuc,inuc,:,:,krad)=reshape(Fi3(inuc,imuc,:,:,krad),npk1,npk1).';
     Fi5(imuc,inuc,:,:,krad)=reshape(Fi5(inuc,imuc,:,:,krad),npk1,npk1).';
     Fi6(imuc,inuc,:,:,krad)=reshape(Fi6(inuc,imuc,:,:,krad),npk1,npk1).';
    end
  end

     end% farn
   end   %nu
     end% farm
  end    %mu
 end %ifalso

%'fine Fi', keyboard

 lPU=length(Pus0);
 Kset=zeros(lPU,lPU);

for krad=1:length(adis)

% Kpp=zeros(npk*nup1,npk*nup1);
 Kpp=Kset;
 Kpm=Kpp;
 Ksp=Kpp;
 Ksm=Kpp;
 KAz1=Kpp;
 KBz1=Kpp;
 Kzer=Kpp;
%' qui ret', keyboard




 if ifalso<=0
 pmat=1:npk;
 imuc=0;
  for jmu=pimu
  imuc=imuc+1;

 imu=jmu-1;
   mu=mbv(jmu);
   if mu==0
    fnor=4*pi;
   else
    fnor=2*pi;
   end
%   pmu=pmat+npk*(imu-1)/pasnu;
   for iarm=1:farm
    ipom=(ceil((imu)/pasnu));
    ipom1=(ceil((farm*(imu)+iarm-1)/pasnu));
    ik1v=Pusas0([ldap(ipom1)+1:ldap(ipom1+1)])-(ipom-1)*npk;
    pmu=[ldap(ipom1)+1:ldap(ipom1+1)];
    lpm=length(pmu);

%   ipom=(ceil(imu/pasnu));
%   pmu=[ldap(ipom)+1:ldap(ipom+1)];
%   ik1v=Pusas0([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;
   inuc=0;
   for jnu=pimu
      inuc=inuc+1;
      inu=jnu-1;
      nu=mbv(jnu);
%      ipom=(ceil(inu/pasnu));
%      pnu=[ldap(ipom)+1:ldap(ipom+1)];
%      ik2v=Pusas0([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;
   for iarn=1:farm
    ipom=(ceil((inu)/pasnu));
    ipom1=(ceil((farm*(inu)+iarn-1)/pasnu));
    ik2v=Pusas0([ldap(ipom1)+1:ldap(ipom1+1)])-(ipom-1)*npk;
    pnu=[ldap(ipom1)+1:ldap(ipom1+1)];
    lpn=length(pnu);
      S=Fi1(imuc,inuc,:,:,krad)+Fi2(imuc,inuc,:,:,krad);
      Kpp(pmu,pnu)=reshape(S(1,1,pmu-ldap(imuc),pnu-ldap(inuc)),lpm,lpn)/fnor;
%      '[imu inu]'
%      [imu inu], pausak
%      keyboard

%      Sr=reshape(S,npk1,npk1)/fnor;
%      Kpp(pmu,pnu)=Sr(ik1v,ik2v);
      S=Fi1(imuc,inuc,:,:,krad)-Fi2(imuc,inuc,:,:,krad);
      Kpm(pmu,pnu)=reshape(S(1,1,pmu-ldap(imuc),pnu-ldap(inuc)),lpm,lpn)/fnor;
%      Sr=reshape(S,npk1,npk1)/fnor;
%      Kpm(pmu,pnu)=Sr(ik1v,ik2v);
      S=Fi3(imuc,inuc,:,:,krad)+Fi4(imuc,inuc,:,:,krad);
      Ksp(pmu,pnu)=reshape(S(1,1,pmu-ldap(imuc),pnu-ldap(inuc)),lpm,lpn)/fnor;
%      Sr=reshape(S,npk1,npk1)/fnor;
%      Ksp(pmu,pnu)=Sr(ik1v,ik2v);
      S=Fi3(imuc,inuc,:,:,krad)-Fi4(imuc,inuc,:,:,krad);
      Ksm(pmu,pnu)=reshape(S(1,1,pmu-ldap(imuc),pnu-ldap(inuc)),lpm,lpn)/fnor;
%      Sr=reshape(S,npk1,npk1)/fnor;
%      Ksm(pmu,pnu)=Sr(ik1v,ik2v);
      if iztm==1
       S=Fi5(imuc,inuc,:,:,krad);
       KAz1(pmu,pnu)=reshape(S(1,1,pmu-ldap(imuc),pnu-ldap(inuc)),lpm,lpn)/fnor;
%       Sr=reshape(S,npk1,npk1)/fnor;
%       KAz1(pmu,pnu)=Sr(ik1v,ik2v);
       S=Fi6(imuc,inuc,:,:,krad);
       KBz1(pmu,pnu)=reshape(S(1,1,pmu-ldap(imuc),pnu-ldap(inuc)),lpm,lpn)/fnor;
%       Sr=reshape(S,npk1,npk1)/fnor;
%       KBz1(pmu,pnu)=Sr(ik1v,ik2v);
      end

    end
    end
   end
   end
 end %ifalso


  if isi==0
   diae=diag(ZEKK);
   diam=diag(ZMKK);
   diae1=1;
   diam1=1;
  else
   diae=diag(sqrt(ZEKK));
   diam=diag(sqrt(ZMKK));
   diae(1,1)=diae(1,1)*fatde1;
   diam(1,1)=diam(1,1)*fatde1;
   diae1=diae;
   diam1=diam;
  end


   for pol=-1:2:1


    KEEp=diae*(Kpp+pol*Ksp)*diae1;
    KEMp=diae*(Kpm+pol*Ksm)*diam1;
    KMEp=diam*(Kpm-pol*Ksm)*diae1;
    KMMp=diam*(Kpp-pol*Ksp)*diam1;

    if pol==-1
     KAmdu=[KEEp segem*KEMp; segem*KMEp KMMp];
     if ianysd==1
      Kaam=[KEEa segem*KEMa; segem*KMEa KMMa];
     end
    else
     KApdu=[KEEp segem*KEMp; segem*KMEp KMMp];
     if ianysd==1
      Kaap=[KEEa segem*KEMa; segem*KMEa KMMa];
     end
    end

   end  %pol


  if iztm==1
   if isi==1
    pol=1;
    diam=sqrt((ZMKK));
    prodz=j*diag(KKv./bev.*diam);
    prodzc=conj(prodz);
    prodz=diag(KKv./bev.*diam);
    prodz=prodz(Pusas0);
    prodzc=prodz;
    KMMpz=0.5*prodz*(KAz1+pol*KBz1)*prodzc;
    pol=-1;
    KMMmz=0.5*prodz*(KAz1+pol*KBz1)*prodzc;
   else
%   disp(' pordz')
%   keyboard
    pol=1;
    prodz=j*KKv./bev;
    prodz=prodz(Pusas0);
    prodz_con=prodz;
    prodz=KKv./bev;
    prodz=prodz(Pusas0);
    prodz_con=prodz;

%    du=diag(.5*ZMKK.*prodz)*(Ktoiiz+KAz1+pol*KBz1)*diag(conj(prodz_con));
    KMMpz=diag(.5*ZMKK.*prodz)*(KAz1+pol*KBz1)*diag(prodz_con);
    pol=-1;
    KMMmz=diag(.5*ZMKK.*prodz)*(KAz1+pol*KBz1)*diag(prodz_con);
   end
   KAzpdu=[Kzer Kzer; Kzer KMMpz];
   KAzmdu=[Kzer Kzer; Kzer KMMmz];
  else
   KAzpdu=0;
   KAzmdu=0;
  end

ndis=krad;
 Kospd (ndis,:,:)=KApdu*alon^2;
 Kosmd (ndis,:,:)=KAmdu*alon^2;
 Koszpd(ndis,:,:)=KAzpdu*alon^2;
 Koszmd(ndis,:,:)=KAzmdu*alon^2;
% disp (' qui uno '), keyboard

 end  % sui raggi

   siz=size(Kospd);
   if length(siz)==3
    si=siz(2:3);
   else
    si=siz;
   end
clear KANzp KANzm KANp KANm
    for iir=pura
     KANp(:,:,iir) = Kospd(iir,:,:);
     KANm(:,:,iir) = Kosmd(iir,:,:);
     KANzp(:,:,iir)=Koszpd(iir,:,:);
     KANzm(:,:,iir)=Koszmd(iir,:,:);
    end

 % disp(' sono in matruas ')
  Kplot1=reshape(Kospd(1,:,:),si);
  Kplot2=reshape(Kosmd(1,:,:),si);
% ' Kplot1', keyboard
 if ifp>0,
%  Kplot3=reshape(Koszpd(1,:,:),si);
%  Kplot4=reshape(Koszmd(1,:,:),si);
%  map(log10(abs(Kplot)))
  map((Kplot1))
   title(' Pol 1 ')
  map((Kplot2))
   title(' Pol -1 ')
  keyboard,
 end,

if mmvet(end)==mm
 clear Fi*
end 

%'fine retoc', keyboard