%' entro lmat_ve ', keyboard
fatde1=1;
imu0=0;
if ifp~=-4
disp(' ')
disp(' ')
end

%ZEv=ZEv*1.2;
%ZMv=ZMv*1.2;
clear A B AB

pimu=2:pasnu:length(mbv)-1;
%keyboard
meun=1;

 mbv1=[-2+nubesi:nubesu+2];
   if iany==2 | ianys==2
    ianysd=1;
   else
    ianysd=0;
   end

%'kmatve', pausak

if sha>1

 ipl1=1;
 ipl2=1;
 for ndis=1:lxi
 icirc0=0;

 aloc=adis(ndis);
 bloc=bdis(ndis);
 alon=min([aloc bloc]);
%  'alon ', keyboard
 a=aloc;
 b=bloc;
 dap=pdi(ndis);
 if sha==5
  Pshi=Psh{ndis};
  iring(ndis)=length(find(Pshi{1}==-1));
  cce=Pshi{5};
 end
% ('ICI 1'), keyboard
  if sha==6
   P=PV{ndis};
   ind_vet=ish(isor(ndis));
  end
%  ' passo shaC ', keyboard
   if shaC==7
    ReS=PAle{1}.Rel*kcav;
    ReSe=PAle{1}.Rm_rel*kcav;
    condi1=( abs(bloc/ReS-1)<.02 );
    condi2=( abs(bloc/ReSe-1)<.02 );
    condi=( condi1*ipl1 | condi2*ipl2);
    if condi1==1
     ipl1=0;
    end
    if condi2==1
     ipl2=0;
    end
    sha_plo=7;
%    sha=4;
    if length(find(shavet==4))>0 
     condi=1;
    end
   else
    condi=1;
   end

  if ifp==-10
%   bloc/kcav
   condi
   ' forma ', pausak
  end

  if igainshape==0
   if ifii(ndis)<=0 & igainshape==0 & condi==1

   forma
%   icosha=icosha+1;
%   lcu=length(xcu);
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


 ifalso=0;
% if sha==3 & bloc/aloc==1
%  ifalso=1;
% end
 if iany==1 & ianys==1
  ifalso=1;
 end

 if sha==3 & bloc/aloc~=1 & iany>1
  disp(' non ancora sistemato per anisotropia ellittica ')
  keyboard
 end

 if sha==2
  shape_sq
  sgimp=ones(size(rv));
 elseif sha==3 & ifalso==0
  shape_el
 elseif sha>3
   if sha==4
    load cce
    if abs(cce)~=0
     Pshi{1}=sha;
     Pshi{2}(1)=bloc;
     Pshi{3}(1)=aloc;
     Pshi{4}(1)=dap;
     Pshi{5}(1)=cce*kcav;
%' wui sce', keyboard
     s_sce
    else
     eval(shape)
    end
  else
        eval(shape)
   end
 end

  if ireturn==1
   %(JMO
  %'dopo shape ', keyboard
  'dopo shape '
  %JMO)
   return
  end

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

 Fi1=zeros(nubes+1,nubes+1,npk,npk);
 Fi2=Fi1;
 Fi3=Fi1;
 Fi4=Fi1;
 Fi5=Fi1;
 Fi6=Fi1;

 rvp=[rv rv(lrv)+diff(rv(lrv-1:lrv))];
 R1=rvp/alon;
 Rp=R1(1:lrv);

 if igint==1
  RdR=sgimp.*wi/alon.*Rp;
 else
  dR=diff(R1);
  RdR=sgimp.*Rp.*dR;
 end
% ifp=5

    if ifp>1
      figure
    end

 ifig=0;
  for imu=pimu
  jmu=imu-1;
  mu=mbv(imu);
   fimu=find(pimu>=imu);
   for inu=pimu
   jnu=inu-1;
   nu=mbv(inu);
%    [jmu jnu], pausak
  if inu>=imu

   if (nu+mu)/2-fix((nu+mu)/2)==0 | pasnu==1
    for ik1=1:npk
     for ik2=1:npk
      F1=.5*besr(ik2,:,inu-1).*besr(ik1,:,imu-1).*RdR.*A(:,jmu,jnu)';
      Fi1(jmu,jnu,ik1,ik2)=sum(F1);
      F2=.5*besr(ik2,:,inu+1).*besr(ik1,:,imu+1).*RdR.*A(:,jmu,jnu)';
      Fi2(jmu,jnu,ik1,ik2)=sum(F2);
      F3=.5*besr(ik2,:,inu-1).*besr(ik1,:,imu+1).*RdR.*B(:,jmu,jnu)';
      Fi3(jmu,jnu,ik1,ik2)=sum(F3);
      F4=.5*besr(ik2,:,inu+1).*besr(ik1,:,imu-1).*RdR.*B(:,jmu,jnu)';
      Fi4(jmu,jnu,ik1,ik2)=sum(F4);
      if iztm==1
       F5d=besr(ik2,:,inu).*besr(ik1,:,imu).*RdR;
       F5=F5d.*A(:,jmu,jnu)';
       Fi5(jmu,jnu,ik1,ik2)=sum(F5);
       F6=F5d.*B(:,jmu,jnu)';
       Fi6(jmu,jnu,ik1,ik2)=sum(F6);
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
     surf(KK,KK,reshape(Fi1(jmu,jnu,:,:),npk,npk)), shading('interp'), view(0,90), colorbar, pausak
     figure
     surf(KK,KK,reshape(Fi2(jmu,jnu,:,:),npk,npk))
     shading('interp'), view(0,90), colorbar, pausak
     figure
     surf(KK,KK,reshape(Fi3(jmu,jnu,:,:),npk,npk))
     shading('interp'), view(0,90), colorbar, pausak
     figure
     surf(KK,KK,reshape(Fi4(jmu,jnu,:,:),npk,npk))
     shading('interp'), view(0,90), colorbar, pausak
     close all
    end
   end  %if
   else
     Fi1(jmu,jnu,:,:)=reshape(Fi1(jnu,jmu,:,:),npk,npk).';
     Fi2(jmu,jnu,:,:)=reshape(Fi2(jnu,jmu,:,:),npk,npk).';
     Fi3(jmu,jnu,:,:)=reshape(Fi4(jnu,jmu,:,:),npk,npk).';
     Fi4(jmu,jnu,:,:)=reshape(Fi3(jnu,jmu,:,:),npk,npk).';
     Fi5(jmu,jnu,:,:)=reshape(Fi5(jnu,jmu,:,:),npk,npk).';
     Fi6(jmu,jnu,:,:)=reshape(Fi6(jnu,jmu,:,:),npk,npk).';
   end

   end   %nu
  end    %mu
 end %ifalso

%  'dopo fii ', keyboard

%' fii ', pausak


 nup1=fix(nubes/pasnu)+1;
 Kiie=zeros(npk*nup1,npk*nup1);
 Kany=zeros(npk*nup1,npk*nup1);
 Ktoa=zeros(npk*nup1,npk*nup1);
 Ktoad=zeros(npk*nup1,npk*nup1);
 Ktoa1=zeros(npk*nup1,npk*nup1);
 Ktoa0=zeros(npk*nup1,npk*nup1);
 Ktoa0v=zeros(npk*nup1,npk*nup1);
 dKtoa=0;
 iKa1=0;
 iKa2=0;
 Kany1=zeros(npk*nup1,npk*nup1);
 Kany2=zeros(npk*nup1,npk*nup1);
 Kzer=Kiie;
 Kije=zeros(npk*nup1,npk*nup1);
 Kde=zeros(1,npk*nup1);
 Kanyd=zeros(1,npk*nup1);
 Kanyd1=zeros(1,npk*nup1);
 Kanyd2=zeros(1,npk*nup1);
 Kd1e=zeros(1,npk*nup1);
 Kiim=zeros(npk*nup1,npk*nup1);
 Kijm=zeros(npk*nup1,npk*nup1);
 Kdm=zeros(1,npk*nup1);
 Kd1m=zeros(1,npk*nup1);

 Kdz=zeros(1,npk*nup1);
 Kiiz=zeros(npk*nup1,npk*nup1);



 if ifalso>=0

   bes=real(besselj(mbv1,KK*alon));


 ip=0;
 for imu=pimu
 mu=mbv(imu);
  if mu==0 & length(pimu)==1
   imu0=1;
  end
 jmu=imu+1;
 ja=jmu+1;
 % if mu~=0
    for ip1=1:npk
     ip=ip+1;
     Q=KK(ip1)*alon;
     P1=bes(ip1,jmu+1)^2+bes(ip1,jmu-1)^2;
     P2=bes(ip1,jmu)*(bes(ip1,jmu+2)+bes(ip1,jmu-2));
     Kdia=.5*(P1-P2);
     Kde(ip)=Kdia;
     Kdm(ip)=Kdia;
     P1=bes(ip1,jmu+1)^2-bes(ip1,jmu-1)^2;
     P2=bes(ip1,jmu)*(bes(ip1,jmu+2)-bes(ip1,jmu-2));
     Kdia=.5*(P1-P2);
     Kd1e(ip)=Kdia;
     Kd1m(ip)=Kdia;
     if ianysd==1
      Kanyd(ip)=0.5*(bes(ip1,ja)^2-bes(ip1,ja-1)*bes(ip1,ja+1));
 %     if ip1==2,
 %     [ja, Kanyd(ip)]
 %     keyboard, end
      if mu==1
       ja1=imu;
       Kanyd1(ip)=0.5*(bes(ip1,ja1)^2+bes(ip1,ja1-1)^2);
      end
      if mu==0
       Kanyd2(ip)=Kanyd(ip);
       Kanyd(ip)=0;
      end
     end
     for ip2=ip1+1:npk
      Q1=KK(ip2)*alon;
      P1=bes(ip2,jmu-1)*bes(ip1,jmu-2)+bes(ip2,jmu+1)*bes(ip1,jmu);
      P2=bes(ip1,jmu-1)*bes(ip2,jmu-2)+bes(ip1,jmu+1)*bes(ip2,jmu);
      Kadia=1/(Q1^2-Q^2)*(Q*P1-Q1*P2);
      Kiie(ip,ip+ip2-ip1)=Kadia;
      Kiim(ip,ip+ip2-ip1)=Kadia;
      P1=-bes(ip2,jmu-1)*bes(ip1,jmu-2)+bes(ip2,jmu+1)*bes(ip1,jmu);
      P2=-bes(ip1,jmu-1)*bes(ip2,jmu-2)+bes(ip1,jmu+1)*bes(ip2,jmu);
      Kadia=1/(Q1^2-Q^2)*(Q*P1-Q1*P2);
      Kije(ip,ip+ip2-ip1)=Kadia;
      Kijm(ip,ip+ip2-ip1)=Kadia;
      if ianysd==1
       P1=bes(ip2,ja)*bes(ip1,ja-1);
       P2=bes(ip1,ja)*bes(ip2,ja-1);
       Kaa=1/(Q1^2-Q^2)*(Q*P1-Q1*P2);
       if mu==0
        Kany2(ip,ip+ip2-ip1)=Kaa;
        iKa2=1;
       else
        Kany(ip,ip+ip2-ip1)=Kaa;
       end
       if mu==1
        ja1=imu;
        P1=bes(ip2,ja1)*bes(ip1,ja1-1);
        P2=bes(ip1,ja1)*bes(ip2,ja1-1);
        Kaa=1/(Q1^2-Q^2)*(Q*P1-Q1*P2);
        Kany1(ip,ip+ip2-ip1)=Kaa;
        iKa1=1;
 %      disp('non diag')
 %      keyboard
       end
      end
     end
    end
 %   if mu==1
 %     Kanyd=Kanyd;
 %     Kany=Kany;
 %     Kp=(Kany+Kany.'+diag(Kanyd))/4;
 %     Ktoad(:,npk+1:length(Kanyd))=Kp(:,1:length(Kanyd)-npk);
 %%     dKtoa=Ktoad.';
 %   end
 end %mu
    Ktoiie=(Kiie+Kiie'+diag(Kde))/4;
    Ktoije=-(Kije+Kije'+diag(Kd1e))/4;
    Ktoiim=(Kiim+Kiim'+diag(Kdm))/4;
    Ktoijm=-(Kijm+Kijm'+diag(Kd1m))/4;

  if iztm==1
   ip=0;
   for imu=pimu
   mu=mbv(imu);
   jmu=imu+1;
 %% if mu~=0
   if mu~=-20
     fbe=[1 1 1 1];
    else
     if p==1
      fbe=[0 0 1 0];
     elseif p==-1
      fbe=[1 0 0 0];
     end
    end
     for ip1=1:npk
      ip=ip+1;
      Q=KK(ip1)*alon;
      P1=bes(ip1,jmu)^2;
      P2=bes(ip1,jmu-1)*bes(ip1,jmu+1);
      Kdia=.5*(P1-P2);
      Kdz(ip)=Kdia*fbe(3);
      for ip2=ip1+1:npk
       Q1=KK(ip2)*alon;
       P1=bes(ip2,jmu)*bes(ip1,jmu-1);
       P2=bes(ip1,jmu)*bes(ip2,jmu-1);
       Kadia=1/(Q1^2-Q^2)*(Q*P1-Q1*P2);
       Kiiz(ip,ip+ip2-ip1)=Kadia*fbe(3);
      end
     end
   end %mu

   Ktoiiz=(Kiiz+Kiiz'+diag(Kdz));
  end % iztm

 else
    Ktoiie=0;
    Ktoije=0;
    Ktoiim=0;
    Ktoijm=0;
    Ktoiiz=0;

 end  % ifalso



%  disp(' pausa Kij0')
%  keyboard
 % pausak


 Kpp=zeros(npk*nup1,npk*nup1);
 Kpm=Kpp;
 Ksp=Kpp;
 Ksm=Kpp;
 KAz1=Kpp;
 KBz1=Kpp;


 if ifalso<=0
 pmat=1:npk;
  for jmu=pimu
   imu=jmu-1;
   mu=mbv(jmu);
   if mu==0
    fnor=4*pi;
   else
    fnor=2*pi;
   end
   pmu=pmat+npk*(imu-1)/pasnu;
   for jnu=pimu
      inu=jnu-1;
      nu=mbv(jnu);
      pnu=pmat+npk*(inu-1)/pasnu;
      S=Fi1(imu,inu,:,:)+Fi2(imu,inu,:,:);
      Kpp(pmu,pnu)=reshape(S,npk,npk)/fnor;
      S=Fi1(imu,inu,:,:)-Fi2(imu,inu,:,:);
      Kpm(pmu,pnu)=reshape(S,npk,npk)/fnor;
      S=Fi3(imu,inu,:,:)+Fi4(imu,inu,:,:);
      Ksp(pmu,pnu)=reshape(S,npk,npk)/fnor;
      S=Fi3(imu,inu,:,:)-Fi4(imu,inu,:,:);
      Ksm(pmu,pnu)=reshape(S,npk,npk)/fnor;
      if iztm==1
       S=Fi5(imu,inu,:,:);
       KAz1(pmu,pnu)=reshape(S,npk,npk)/fnor;
       S=Fi6(imu,inu,:,:);
       KBz1(pmu,pnu)=reshape(S,npk,npk)/fnor;
      end

    end
   end
 end %ifalso

    if ianysd==1
     Kp=(Kany+Kany.'+diag(Kanyd))/4;
     Ktoad(:,npk+1:length(Kanyd))=Kp(:,1:length(Kanyd)-npk);
     Ktoa=Ktoad;
     if iKa1==1
      Ktoa1=(Kany1+Kany1'+diag(Kanyd1))/4;
     else
      Ktoa1=0;
     end
     if iKa2==1
      Kp=(Kany2+Kany2'+diag(Kanyd2))/4;
      Ktoad(:,npk+1:length(Kanyd))=Kp(:,1:length(Kanyd)-npk);
      Ktoa0v=Ktoad*2;
      Ktoa0=Ktoad;
     else
      Ktoa0=0;
      Ktoa0v=0;
     end
    end



  if isi==0
   diae=diag(ZEv.*KKv);
   diam=diag(ZMv.*KKv);
   diae1=1;
   diam1=1;
  else
   diae=diag(sqrt(ZEv.*KKv));
   diam=diag(sqrt(ZMv.*KKv));
   diae(1,1)=diae(1,1)*fatde1;
   diam(1,1)=diam(1,1)*fatde1;
   diae1=diae;
   diam1=diam;
  end

    matpolp0=ones(size(diae));
    matpolp1=ones(size(diae));
    if imu0==1
%     fimam=1:length(KK);
%     matpolp0(fimam,fimam)=0;
     if length(ldapu)>0
      fimam=1:ldapu(2);
     else
      fimam=1:length(KK);
     end
     matpolp0(fimam,:)=0;
     matpolp0(:,fimam)=0;
    end

   for pol=-1:2:1
    if pol==-1
     matpolm=matpolp0;
     matpolp=matpolp1;
    else
     matpolp=matpolp0;
     matpolm=matpolp1;
    end

    KEEp=diae*(Ktoiie.*matpolp+Kpp+pol*Ksp)*diae1;
    KEMp=diae*(Ktoije.*matpolp+Kpm+pol*Ksm)*diam1;
    KMEp=diam*(Ktoijm.*matpolm+Kpm-pol*Ksm)*diae1;
    KMMp=diam*(Ktoiim.*matpolm+Kpp-pol*Ksp)*diam1;

    if ianysd==1
     pm1=pol-1;
     pp1=pol+1;
     KEEa=diae*(+Ktoa+Ktoa.'+dKtoa+pol*Ktoa1-pm1*Ktoa0-pm1*Ktoa0v.')*diae1;
     KEMa=diae*(+Ktoa-Ktoa.'-dKtoa+pol*Ktoa1-pm1*Ktoa0-pp1*Ktoa0v.')*diam1;
     KMEa=diam*(-Ktoa+Ktoa.'+dKtoa+pol*Ktoa1-pp1*Ktoa0-pm1*Ktoa0v.')*diae1;
     KMMa=diam*(-Ktoa-Ktoa.'-dKtoa+pol*Ktoa1-pp1*Ktoa0-pp1*Ktoa0v.')*diam1;
    end

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
    diam=sqrt((ZMv.*KKv));
    prodz=j*diag(KKv./bev.*diam);
    prodzc=conj(prodz);
    prodz=diag(KKv./bev.*diam);
    prodzc=prodz;
    KMMpz=0.5*prodz*(Ktoiiz+KAz1+pol*KBz1)*prodzc;
    pol=-1;
    KMMmz=0.5*prodz*(Ktoiiz+KAz1+pol*KBz1)*prodzc;
   else
%   disp(' pordz')
%   keyboard
    pol=1;
    prodz=j*KKv./bev;
    prodz_con=prodz;
    prodz=KKv./bev;
    prodz_con=prodz;
    KMMpz=diag(.5*ZMv.*KKv.*prodz)*(Ktoiiz+KAz1+pol*KBz1)*diag(prodz_con);
    pol=-1;
    KMMmz=diag(.5*ZMv.*KKv.*prodz)*(Ktoiiz+KAz1+pol*KBz1)*diag(prodz_con);
   end
   KAzpdu=[Kzer Kzer; Kzer KMMpz];
   KAzmdu=[Kzer Kzer; Kzer KMMmz];
  else
   KAzpdu=0;
   KAzmdu=0;
  end

%' fine ', pausak
  if icirc0==1
   lxivet=ndis;
   matve_ci
   KApdu=KApdu+KApsub/alon^2;
   KAmdu=KAmdu+KAmsub/alon^2;
   KAzpdu=KAzpdu+KAzpsub/alon^2;
   KAzmdu=KAzmdu+KAzmsub/alon^2;
  end


 Kospd (ndis,:,:)=KApdu*alon^2;
 Kosmd (ndis,:,:)=KAmdu*alon^2;
 Koszpd(ndis,:,:)=KAzpdu*alon^2;
 Koszmd(ndis,:,:)=KAzmdu*alon^2;
% disp (' qui uno '), keyboard

%' igainshape vett ', keyboard
  if igainshape==0
   if aiat(1)==ayi(ndis) & sha==shav.a(1)
    KAp=KApdu*alon^2;
    KAm=KAmdu*alon^2;
    KAzp=KAzpdu*alon^2;
    KAzm=KAzmdu*alon^2;
   end
  end

 end  % sui raggi

   siz=size(Kospd);
   if length(siz)==3
    si=siz(2:3);
   else
    si=siz;
   end

    for iir=1:lxi
     Kosp(:,:,iir) = Kospd(iir,:,:);
     Kosm(:,:,iir) = Kosmd(iir,:,:);
     Koszp(:,:,iir)=Koszpd(iir,:,:);
     Koszm(:,:,iir)=Koszmd(iir,:,:);
    end

   if sha==5
    ifiring=find(iring>0);
    ifiring0=find(iring==0);
    if length(ifiring)==1 & length(ifiring0)==1
     iir=ifiring;
     iir0=ifiring0;
     Kosp(:,:,iir) = Kospd(iir,:,:)- Kospd(iir0,:,:);
     Kosm(:,:,iir) = Kosmd(iir,:,:)- Kosmd(iir0,:,:);
     Koszp(:,:,iir)=Koszpd(iir,:,:)-Koszpd(iir0,:,:);
     Koszm(:,:,iir)=Koszmd(iir,:,:)-Koszmd(iir0,:,:);
    end
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

  clear Kospd Kosmd Koszpd Koszmd

 if iany>=2
  Kosanp=Kaap/2;
  Kosanm=Kaam/2;
 end

% disp(' qui  xiT  Vett'), keyboard
 %&&&&&&&&&&&&&&&

  if igainshape==1 & igainshapeT==0
   isu=[1:lxi];
    Kdu1=0;
    Kdu2=0;
    Kdu3=0;
    Kdu4=0;
    for nsu=isu
     Kdu1=Kdu1+Kosp(:,:,nsu)*yiN(nsu);
     Kdu2=Kdu2+Kosm(:,:,nsu)*yiN(nsu);
     Kdu3=Kdu3+Koszp(:,:,nsu)*yiN(nsu);
     Kdu4=Kdu4+Koszm(:,:,nsu)*yiN(nsu);
    end
    KAp=Kdu1;
    KAm=Kdu2;
    KAzp=Kdu3;
    KAzm=Kdu4;

    if ianti_gui==1
     Kdu1=0;
     Kdu2=0;
     Kdu3=0;
     Kdu4=0;
     for nsu=isu
      Kdu1=Kdu1+Kosp(:,:,nsu)*yiN_ag(nsu);
      Kdu2=Kdu2+Kosm(:,:,nsu)*yiN_ag(nsu);
      Kdu3=Kdu3+Koszp(:,:,nsu)*yiN_ag(nsu);
      Kdu4=Kdu4+Koszm(:,:,nsu)*yiN_ag(nsu);
     end
     KAp_ag=Kdu1;
     KAm_ag=Kdu2;
     KAzp_ag=Kdu3;
     KAzm_ag=Kdu4;
    end

  end

%  if igainshapeT==1
 if igainshape==2 & igainshapeT==0
  isu=[1:lxi/2];

    Kdu1=0;
    Kdu2=0;
    Kdu3=0;
    Kdu4=0;
    for nsu=isu
     Kdu1=Kdu1+Kosp(:,:,nsu)*yiN(nsu);
     Kdu2=Kdu2+Kosm(:,:,nsu)*yiN(nsu);
     Kdu3=Kdu3+Koszp(:,:,nsu)*yiN(nsu);
     Kdu4=Kdu4+Koszm(:,:,nsu)*yiN(nsu);
    end
    KAp=Kdu1;
    KAm=Kdu2;
    KAzp=Kdu3;
    KAzm=Kdu4;


    if ianti_gui==1
     Kdu1=0;
     Kdu2=0;
     Kdu3=0;
     Kdu4=0;
     for nsu=isu
      Kdu1=Kdu1+Kosp(:,:,nsu)*yiN_ag(nsu);
      Kdu2=Kdu2+Kosm(:,:,nsu)*yiN_ag(nsu);
      Kdu3=Kdu3+Koszp(:,:,nsu)*yiN_ag(nsu);
      Kdu4=Kdu4+Koszm(:,:,nsu)*yiN_ag(nsu);
     end
     KAp_ag=Kdu1;
     KAm_ag=Kdu2;
     KAzp_ag=Kdu3;
     KAzm_ag=Kdu4;
    end
   end
  if igainshapeT==1

   s=size(yiT);
   isuT=[1:s(1)];
   for iz=1:s(2)
    Ktoiiep=0;
    Ktoiiem=0;
    Kdu1=0;
    Kdu2=0;
    Kdu3=0;
    Kdu4=0;
    indis=0;
    for ndis=isuT
     indis=indis+1;
     if iresK==0
      ndis1=ndis;
     else
      ndis1=ndis+(iz-1)*s(1);
     end
     Ktoiiep=Ktoiiep+Kosp(:,:,ndis1)*yiT(ndis,iz);
     Ktoiiem=Ktoiiem+Kosm(:,:,ndis1)*yiT(ndis,iz);
     Kdu1=Kdu1+Kosp(:,:,ndis1)*yiT(ndis,iz);
     Kdu2=Kdu2+Kosm(:,:,ndis1)*yiT(ndis,iz);
     Kdu3=Kdu3+Koszp(:,:,ndis1)*yiT(ndis,iz);
     Kdu4=Kdu4+Koszm(:,:,ndis1)*yiT(ndis,iz);
    end  %ndis
    KTempp(:,:,iz)=Kdu1;
    KTempm(:,:,iz)=Kdu2;
    KTemzp(:,:,iz)=Kdu3;
    KTemzm(:,:,iz)=Kdu4;
   end
%  ' qui KTE', keyboard
  end

  if exist('KA')
%   KAp=KA;
%   KAm=KA;

   KApsub=KAp;
   KAmsub=KAm;
  end

% ' sono qui ', keyboard


else %sha==1

 lxivet=[1:lxi];
 matve_ci
 Kosp=Kospsub;
 Kosm=Kosmsub;
 Koszp=Koszpsub;
 Koszm=Koszmsub;
%' dopo matve_ci ' , keyboard

 if exist('KAp')
  KAp=KApsub;
  KAm=KAmsub;
  if iztm==1
   KAzp=KAzpsub;
   KAzm=KAzmsub;
  end
  Kplot1=KAp;
  Kplot2=KAm;
 else
  Kplot1=Kosp;
  Kplot2=Kosm;
 end

end  %sha





si=size(Kiie);
%si2=2*si;
Znor=[ZEv; ZMv];

%ifp=2
if ifp>=-3, disp(' Kmat_ve '),
 if ifp>1,
  keyboard,
 end,
end

% disp(' Kmat_ve '),
% keyboard
