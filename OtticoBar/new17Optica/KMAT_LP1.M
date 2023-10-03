%  'kmat 1 in ', keyboard
clear Kos KoN
ifalso=0;
pimu=1:pasnu:length(mbv);
meun=0;
clear A B AB

%' Entro K Prova ', pausak

if sha>1

 for ndis=1:lxi

  icirc0=0;
  aloc=adis(ndis);
  bloc=bdis(ndis);
  alon=min([aloc bloc]);
%  alon=aloc;
  a=aloc;
  b=bloc;
  dap=pdi(ndis);
  if sha==5
   Pshi=Psh{ndis};
   iring(ndis)=length(find(Pshi{1}==-1));
   cce=Pshi{5};
   if ifp>=0
    figure, plot(cce,'ro'), pausak
   end
  end

% ('ICI LP'), keyboard
  if sha==6
   P=PV{ndis};
  end
%  if sha==5
%   ityar=ceiv(ndis
%   nar=ceiv(ndis,2);
%   cce=ceiv(ndis,3:nar+2);
%  end
%'forma'
%keyboard
  if igainshape==0
   if ifii(ndis)<=0 & igainshape==0
%   ifidis=find(ifield<0);
%   afiel=aitot(ifidis);
%   'afiel', keyboard
%   if (ifii(ndis)<=0 | length(find(afiel==a))>0) & igainshape==0
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


  if sha==2
   shape_sq
   sgimp=ones(size(rv));
  elseif sha==3
   shape_el
  elseif sha>3
   if sha==4
    global cce
    if length(cce)==0
     load cce
    end
    if abs(cce)~=0
     Pshi{1}=sha;
     Pshi{2}(1)=bloc;
     Pshi{3}(1)=aloc;
     Pshi{4}(1)=dap;
     Pshi{5}(1)=cce*kcav;
     s_sce
    else
     eval(shape)
    end
   else
    eval(shape)
   end
  end

   lKK=length(KK);
   lrv=length(rv);
   z=KK*rv;
   clear besr

   for iz=1:length(mbv)
    nu=mbv(iz);
    bes0=besselj(nu,z);
    besr(:,:,iz)=real(bes0);
   end

  Fi1=zeros(nubes+1,nubes+1,npk,npk);
  Fi2=Fi1;
  if ipolar==0
   Fi3=Fi1;
   Fi4=Fi1;
  end

  rvp=[rv rv(lrv)+diff(rv(lrv-1:lrv))];
  R1=rvp/alon;
  Rp=R1(1:lrv);

  if igint==1
   RdR=wi/alon.*Rp;
  else
   dR=diff(R1);
   RdR=Rp.*dR;
  end
%  disp(' kmat_lp'), keyboard

  ifig=0;
   for imu=pimu
   ipom=(ceil(imu/pasnu));
   ik1v=Pusa([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;
   jmu=imu;
   mu=mbv(imu);
    for inu=pimu
    ipom=(ceil(inu/pasnu));
    ik2v=Pusa([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;
    jnu=inu;
    nu=mbv(inu);
    if (nu+mu)/2-fix((nu+mu)/2)==0 | pasnu==1
     for ik1=ik1v
      for ik2=ik2v
       F1=besr(ik2,:,inu).*besr(ik1,:,imu).*RdR.*A(:,jmu,jnu)';
       Fi1(jmu,jnu,ik1,ik2)=sum(F1);
       F2=besr(ik2,:,inu).*besr(ik1,:,imu).*RdR.*B(:,jmu,jnu)';
       Fi2(jmu,jnu,ik1,ik2)=sum(F2);
       if ipolar==0 & sha==4
        F3=besr(ik2,:,inu).*besr(ik1,:,imu).*RdR.*D(:,jmu,jnu)';
        Fi3(jmu,jnu,ik1,ik2)=sum(F3);
        F4=besr(ik2,:,inu).*besr(ik1,:,imu).*RdR.*C(:,jmu,jnu)';
        Fi4(jmu,jnu,ik1,ik2)=sum(F4);
       end
      end %ik2
     end  %ik1
     if ifp>=2
       if ipolar==0 & sha==4
        figure, plot(r,F1,'.-',r,F2,'.-',r,F3,'.-',r,F4,'.-')
       else
        figure, plot(r,F1,'.-',r,F2,'.-')
        ' integrale = '
        [sum(F1) sum(F2)]
       end
      pausak
     end
     if ifig==1
      figure
      surf(KK,KK,reshape(Fi1(jmu,jnu,:,:),npk,npk)), shading('interp'), view(0,90), colorbar, pausak
      figure
      surf(KK,KK,reshape(Fi2(jmu,jnu,:,:),npk,npk))
      shading('interp'), view(0,90), colorbar, pausak
%      close all
     end
    end  %if
    end   %nu
   end    %mu

  nup1=fix(nubes/pasnu)+1;
  Kiie=zeros(npk*nup1,npk*nup1);

  if ifalso>=0

   mbv1=[-1+nubesi:nubesu+1];
   bes=real(besselj(mbv1,KK*alon));

  Kde=zeros(1,npk*nup1);

  %

  ip=0;
  for imu=pimu
   ipom=(ceil(imu/pasnu));
   ik1v=Pusa([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;
  mu=mbv(imu);
  jmu=imu+1;
     for ip1=ik1v
      ip=ip1+(ipom1-1)*npk;
%     for ip1=1:npk
%      ip=ip+1;
      Q=KK(ip1)*alon;
      P1=bes(ip1,jmu)^2;
      P2=bes(ip1,jmu-1)*bes(ip1,jmu+1);
      Kdia=.5*(P1-P2);
      Kde(ip)=Kdia;
%      for ip2=ip1+1:npk
      for ip2=ip1+1:ik1v(end)
       Q1=KK(ip2)*alon;
       P1=bes(ip2,jmu)*bes(ip1,jmu-1);
       P2=bes(ip1,jmu)*bes(ip2,jmu-1);
       Kadia=1/(Q1^2-Q^2)*(Q*P1-Q1*P2);
       Kiie(ip,ip+ip2-ip1)=Kadia;
      end
     end
  end %mu


     Ktoiie=(Kiie+Kiie'+diag(Kde))/2;

  % disp(' pausa Kij0')
  % pausak
  else
   Ktoiie=0;
  end  %ifalso


  Kpp=zeros(npk*nup1,npk*nup1);
  Kpm=zeros(npk*nup1,npk*nup1);
  if ipolar==0 & sha==4
   Kppa=Kpp;
   Kpma=Kpp;
  end

%' dopo prima ', keyboard

  pmat=1:npk;
  fnor=2*pi;
  for imu=1:pasnu:nubes+1
   mu=mbv(imu);
   if mu==0
    fnor=4*pi;
   else
    fnor=2*pi;
   end
   pmu=pmat+npk*(imu-1)/pasnu;
   for inu=1:pasnu:nubes+1
    pnu=pmat+npk*(inu-1)/pasnu;
     S=Fi1(imu,inu,:,:);
     Kpp(pmu,pnu)=reshape(S,npk,npk)/fnor;
     S=Fi2(imu,inu,:,:);
     Kpm(pmu,pnu)=reshape(S,npk,npk)/fnor;
       if ipolar==0 & sha==4
        S=Fi3(imu,inu,:,:);
        Kppa(pmu,pnu)=reshape(S,npk,npk)/fnor;
        S=Fi4(imu,inu,:,:);
        Kpma(pmu,pnu)=reshape(S,npk,npk)/fnor;
       end
   end
  end

%' dopo fine ', keyboard

%    pol=1;
%    KApdu=(Ktoiie+Kpp+pol*Kpm);
%    pol=-1;
%    KAmdu=(Ktoiie+Kpp+pol*Kpm);

   if isi==0
    dia=diag(ZEv.*KKv*alon^2);
    pol=1;
    KApdu=dia*(Ktoiie+Kpp+pol*Kpm);
    pol=-1;
    Ktod=Ktoiie;
    if mbv(1)==0 & length(Ktod)>1
      Ktod(1:nk1max,1:nk1max)=0;
    end
    KAmdu=dia*(Ktod+Kpp+pol*Kpm);
       if ipolar==0
        K1=(Ktoiie+Kpp+Kpm);
        Ktod=Ktoiie;
        if mbv(1)==0 & length(Ktod)>1
         Ktod(1:nk1max,1:nk1max)=0;
        end
        K2=(Ktod+Kpp-Kpm);
        if sha==4
%        K3=Kppa;
%        K4=Kpma;
         K3=(Kppa-Kpma);
         K4=(Kppa+Kpma);
        else
         K3=zeros(size(K1));
         K4=K3;
        end
        dia2=diag([ZEv.*KKv; ZEv.*KKv]*alon^2);
        KApdua=dia2*[K1 K3; K4 K2];
       end
   else
    dia=diag(sqrt(ZEv.*KKv))*alon;
    pol=1;
    KApdu=dia*(Ktoiie+Kpp+pol*Kpm)*dia;
    pol=-1;
    Ktod=Ktoiie;
    if mbv(1)==0 & length(Ktod)>1
      Ktod(1:nk1max,1:nk1max)=0;
    end
    KAmdu=dia*(Ktod+Kpp+pol*Kpm)*dia;
       if ipolar==0
        K1=(Ktoiie+Kpp+Kpm);
        Ktod=Ktoiie;
        if mbv(1)==0 & length(Ktod)>1
         Ktod(1:nk1max,1:nk1max)=0;
        end
        K2=(Ktod+Kpp-Kpm);
        if sha==4
%        K3=Kppa;
%        K4=Kpma;
         K3=(Kppa-Kpma);
         K4=(Kppa+Kpma);
        else
         K3=zeros(size(K1));
         K4=K3;
        end
        dia2=diag(sqrt([ZEv.*KKv; ZEv.*KKv])*alon);
        KApdua=dia2*[K1 K3; K4 K2]*dia2;
       end
   end
       if ipolar==0
        KApdu=KApdua;
       end

%   '[icirc0 sha ifalso aloc rv(1)]'
%   [icirc0 sha ifalso aloc rv(1)]
%   pausak
   if icirc0==1
    lxivet=ndis;
    ' p1'
    pausak
    lp_ci1
    KAp0=KApdu;
    KApdu=KApdu+KApsub;
    KAmdu=KAmdu+KAmsub;
   end

  Kplot1=KApdu;
  Kplot2=KAmdu;
  if ifp>=1
    if icirc0==1
     figure,
     surf(KAp0), shading('interp'), view(0,90), colorbar,
     figure,
     surf(KApsub), shading('interp'), view(0,90), colorbar,
    end
    figure,
   surf(Kplot1), shading('interp'), view(0,90), colorbar,
   title(' Pol 1 ')
    figure,
   surf(Kplot2), shading('interp'), view(0,90), colorbar,
   title(' Pol -1 '), pausak
  end

%  'ndis'
%  pausak
  Kosp (:,:,ndis)=KApdu;
  Kosm (:,:,ndis)=KAmdu;
  Koszp(:,:,ndis)=KAmdu;
  Koszm(:,:,ndis)=KAmdu;

%' igainshape LP ', keyboard

  if igainshape==0
   if aiat(1)==ayi(ndis) & sha==shav.a(1)
    KAp=KApdu;
    KAm=KAmdu;
    KAzp=KApdu;
    KAzm=KAmdu;
   end
  end


 end  %lxi

% disp(' qui  xiT  LP'), keyboard

  if igainshape==1
   isu=[1:lxi];
   if exist('xiT')
    isuT=[1:length(xiT)];
   end
    Kdu=0;
    for nsu=isu
     Kdu=Kdu+Kosp(:,:,nsu)*yiN(nsu);
    end
    KA=Kdu;

    if ianti_gui==1
     Kdu=0;
     for nsu=isu
      Kdu=Kdu+Kosp(:,:,nsu)*yiN_ag(nsu);
     end
     KA_ag=Kdu;
    end
  end
  if exist('KA')
%   KAp=KA;
%   KAm=KA;

   KApsub=KAp;
   KAmsub=KAm;
  end

else  %sha=1

  lxivet=[1:lxi];

  lp_ci1
%'  dopo newlp_ci ', keyboard

        if ipolar==0
         K3=zeros(size(Kospsub));
         Kospsub=[Kospsub K3; K3 Kosmsub];
         sti=size(Kti);
         if length(sti)>2
%         stu=[sti(2:3) sti(1)];
%         stf=[sti(1) 2*sti(2:3)];
%         Kpa=reshape(Kti,stu);
%         keyboard
%          K3=zeros(size(Kpa));
%          Kpad=[Kpa K3; K3 Kpa];
%          Kti=reshape(Kpad,stf);
           K3=zeros(sti(2),sti(3));
           for ima=1:sti(1)
            Kdu=reshape(Kti(ima,:,:),sti(2),sti(3));
            if isi==0
             Kdu1=diag(ZEv.*KKv)*(Kdu);
            else
             diae=diag(sqrt(ZEv.*KKv));
             Kdu1=diae*(Kdu)*diae;
            end
            Ktiu(ima,:,:)=[Kdu1 K3; K3 Kdu1];
           end
           Kti=Ktiu;
%         ' Kti ', keyboard
         end
         K4=zeros(size(KApsub));
         KApsub=[KApsub K4; K4 KAmsub];
        end


   Kosp=Kospsub;
   Kosm=Kosmsub;

   Koszp=Kospsub;
   Koszm=Kosmsub;


  if exist('KA')
   KAp=KApsub;
   KAm=KAmsub;
   KAzp=KApsub;
   KAzm=KAmsub;
   Kplot1=KAp;
   Kplot2=KAm;
  else
   Kplot1=Kosp;
   Kplot2=Kosm;
  end


end  %sha

  if sha==5
    ifiring=find(iring>0);
    ifiring0=find(iring==0);
%    'ifiring'
%    keyboard
    if length(ifiring)==1 & length(ifiring0)==1
     iir=ifiring;
     iir0=ifiring0;
     Kosp(:,:,iir) =Kosp(:,:,iir) -Kosp(:,:,iir0) ;
     Kosm(:,:,iir) =Kosm(:,:,iir) -Kosm(:,:,iir0) ;
     Koszp(:,:,iir)=Koszp(:,:,iir)-Koszp(:,:,iir0);
     Koszm(:,:,iir)=Koszm(:,:,iir)-Koszm(:,:,iir0);
    end
%    'Kosp'
%    pausak
  end

 if exist('KApdu')
  si=size(KApdu);
 else
  si=size(Kiie);
 end
 si2=si;

 if itetm==2
  matru0u
 end



Ideltaap=0;
Ideltaam=0;
Ideltaa=0;
istatt=0;
flp=1-iLP1;
mr=1;
Ide=diag(ones(1,2*length(be)));

deltany=0;


%  ZEv=mr./sqrt(1-flp*KKv.^2*mr^2);
%  ZMv=mr.*sqrt(1-flp*KKv.^2*mr^2);

   if deltany~=0
    ZEva=mr./sqrt(1-flp*KKva.^2*mr^2);
    ZMva=mr.*sqrt(1-flp*KKva.^2*mr^2);
    if itetm==2
     ZEva=ZMva;
    end
    Idelteac=ZEva/2;
    Ideltmac=ZMva/2;
   end


%  Ideltad=([ZEv])/2;
%  Ideltadz=([ZEv.*KKv.^2./(1-KKv.^2)])/2;
%  KKt=0;
%  if ipolar==0
%   Ideltad=[Ideltad; Ideltad];
%  end
%
% Idelta=diag(Ideltad);
% Ideltaz=diag(Ideltadz);
 if ifp>=-1
  disp(' Kmat_lp '), keyboard
 end
%  disp(' Kmat_lp '), keyboard

