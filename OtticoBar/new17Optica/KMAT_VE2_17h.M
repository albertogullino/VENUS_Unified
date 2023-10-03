ifp=10
      clear besr1
      clear besr2

%'kmat_ve2', keyboard
if sha==6 & ifp==-10
% 'entro kvet2', keyboard
end
if sha==6 
 if isfield(Ps,'igraef_new')==1
  if Ps.igraef_new==2
   ireturn=2;
   return
  end
 end 
% 'entro kvet2', keyboard
end

% 'entro kvet2', keyboard
% 'entro kvet2', keyboard

conde=length(pros<0);
coext= (sha==6 | (conde>0 & shad==7));
if coext==1 & ired_ret==2
 puPa=1;
elseif coext==0 & ired_ret==2
 puPa=2;
end
if ired_ret==3 
 puPa=1;
end 


if exist('puPa')
 Pus0=Par{puPa}.Pus0;
 Pusas0=Par{puPa}.Pusas0;
 ldapu=Par{puPa}.ldapu;
 ldap=Par{puPa}.ldap;
 mbv=Par{puPa}.mbv;
 nubes=Par{puPa}.nubes;
 nubesi=Par{puPa}.nubesi;
 nubesu=Par{puPa}.nubesu;
end



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

pimu=2:pasnu:length(mbv)-1;
%pimu=mbv(2):pasnu:mbv(end-1);

'entro ve2', keyboard

meun=1;

 mbv1=[-2+nubesi:nubesu+2];
   if iany==2 | ianys==2
    ianysd=1;
   else
    ianysd=0;
   end

%'kmatve', pausak

if sha>1

 if sha==11
  sha_vortex
 end
global cce
 for ndis=1:lxi
 icirc0=0;

 aloc=adis(ndis);
 bloc=bdis(ndis);
 alon=min([aloc bloc]);

 a=aloc;
 b=bloc;
 dap=pdi(ndis);
 %' adis ', keyboard
 if sha==5
  Pshi=Psh{ndis};
  iring(ndis)=length(find(Pshi{1}==-1));
  cce=Pshi{5};
 end
 %('ICI 1'), keyboard
  if sha==6
   P=PV{ndis};
   ind_vet=ish(isor(ndis));
  end
  if igainshape==0
   
  
%   if ifii(ndis)<=0 & igainshape==0 & P.D>.1
   if ifii(ndis)<=0 & igainshape==0

%' prima di forma', keyboard
   forma

%   icosha=icosha+1;
%   lcu=length(xcu);
' dopo forma', keyboard
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
 elseif sha>3 & sha~=11
  eval(shape)
%  if sha==6 & igraef_new==2
  if sha==6 & igr_app==1 & iretpiano==1
%   iret=find(shavet(:,1)==6);
%   shavet(iret,1)=0;  
   ireturn=2;
%   ' passo ireturn=2', keyboard
   return
  end
 end
'dopo shape', keyboard
  if ireturn==1
   if ifirstgrat==0
    k_anip
   end
   if ifp~=-4
%    'dopo k_anip ', keyboard
   end
   return
  end

'dopo shag', keyboard
 if ifalso<=0
  lKK=length(KK);
  lrv=length(rv);
  z=KK*rv;

  clear besr besz

%' qui prima di bers', keyboard
%' qui prima di bers', keyboard
tic
iob=0;
if iob==1
  for iz=1:length(mbv)
   nu=mbv(iz);
   bes0=real(besselj(nu,z));
   besr(:,:,iz)=bes0;
  end
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
  for imu=pimu
 %  ipom=(ceil(imu/pasnu));
 %  ik1v=Pusas0([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;
   for iarm=1:farm
    ipom=(ceil((imu-1)/pasnu));
    ipom1=(ceil((farm*(imu-1)+iarm-1)/pasnu));
    ik1v=Pusas0([ldap(ipom1)+1:ldap(ipom1+1)])-(ipom-1)*npk;

     if iob==0
     'dentro iob', keyboard
%      clear besr1
      z1=KK(ik1v)*rv;
      imuvero=mbv(imu);
      dmu=imuvero-imu;
%      for iz=imuvero+[-1 0 1]
      for iz=imu+[-1 0 1]
       izver=iz+dmu;
       bes0=real(besselj(izver,z1));
       besr1(:,:,iz)=bes0;
      end   
     end
  jmu=imu-1;
  mu=mbv(imu);
   fimu=find(pimu>=imu);
   for inu=pimu
%    ipom=(ceil(inu/pasnu));
%    ik2v=Pusas0([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;
    for iarn=1:farm
     ipom=(ceil((inu-1)/pasnu));
     ipom1=(ceil((farm*(inu-1)+iarn-1)/pasnu));
     ik2v=Pusas0([ldap(ipom1)+1:ldap(ipom1+1)])-(ipom-1)*npk;

   jnu=inu-1;
   nu=mbv(inu);
%    [jmu jnu], pausak
     if iob==0
%     'dentro iob', keyboard
%      clear besr2
      z1=KK(ik2v)*rv;
%      for iz=inu+[-1 0 1]
      imuvero=mbv(inu);
      dnu=imuvero-mbv(inu);
%      for iz=imuvero+[-1 0 1]
      for iz=inu+[-1 0 1]
       izver=iz-1+dnu;       
       izver=iz+dnu;       
       bes0=real(besselj(izver,z1));      
%       bes0=real(besselj(iz-1,z1));
       besr2(:,:,iz)=bes0;
      end     
     end
  if inu>=imu

   if (nu+mu)/2-fix((nu+mu)/2)==0
%    for ik1=1:npk
%     for ik2=1:npk

     for ik1=ik1v
      for ik2=ik2v
       if iob==1      
        besr1=besr;
        besr2=besr;
        ik1p=ik1;
        ik2p=ik2;
       else
        ik1p=find(ik1==ik1v);
        ik2p=find(ik2==ik2v);
       end
       
        F1=.5*besr2(ik2p,:,inu-1).*besr1(ik1p,:,imu-1).*RdR.*A(:,jmu,jnu)';
        Fi1(jmu,jnu,ik1,ik2)=sum(F1);
        %'F2', keyboard
        F2=.5*besr2(ik2p,:,inu+1).*besr1(ik1p,:,imu+1).*RdR.*A(:,jmu,jnu)';
        Fi2(jmu,jnu,ik1,ik2)=sum(F2);
        F3=.5*besr2(ik2p,:,inu-1).*besr1(ik1p,:,imu+1).*RdR.*B(:,jmu,jnu)';
        Fi3(jmu,jnu,ik1,ik2)=sum(F3);
        F4=.5*besr2(ik2p,:,inu+1).*besr1(ik1p,:,imu-1).*RdR.*B(:,jmu,jnu)';
        Fi4(jmu,jnu,ik1,ik2)=sum(F4);
        if iztm==1
         F5d=besr2(ik2p,:,inu).*besr1(ik1p,:,imu).*RdR;
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

     end% farn
   end   %nu
     end% farm
  end    %mu
 end %ifalso

toc
'dopo Fiii', keyboard
%'dopo Fiii', keyboard
%'dopo Fiii', keyboard

 nup1=fix(nubes/pasnu)+1;
% ' nup1 ', keyboard
 lPU=length(Pus0);
 Kset=zeros(lPU,lPU);
 Vset=zeros(1,lPU);
 Kiie=Kset;
 Kiie=Kset;
 Kany=Kset;
 Ktoa=Kset;
 Ktoad=Kset;
 Ktoa1=Kset;
 Ktoa0=Kset;
 Ktoa0v=Kset;
 dKtoa=0;
 iKa1=0;
 iKa2=0;
 Kany1=Kset;
 Kany2=Kset;
 Kzer=Kiie;
 Kije=Kset;
 Kde=Vset;
 Kanyd=Vset;
 Kanyd1=Vset;
 Kanyd2=Vset;
 Kd1e=Vset;
 Kiim=Kset;
 Kijm=Kset;
 Kdm=Vset;
 Kd1m=Vset;

 Kdz=Vset;
 Kiiz=Kset;

% Kiie=zeros(npk*nup1,npk*nup1);
% Kany=zeros(npk*nup1,npk*nup1);
% Ktoa=zeros(npk*nup1,npk*nup1);
% Ktoad=zeros(npk*nup1,npk*nup1);
% Ktoa1=zeros(npk*nup1,npk*nup1);
% Ktoa0=zeros(npk*nup1,npk*nup1);
% Ktoa0v=zeros(npk*nup1,npk*nup1);
% dKtoa=0;
% iKa1=0;
% iKa2=0;
% Kany1=zeros(npk*nup1,npk*nup1);
% Kany2=zeros(npk*nup1,npk*nup1);
% Kzer=Kiie;
% Kije=zeros(npk*nup1,npk*nup1);
% Kde=zeros(1,npk*nup1);
% Kanyd=zeros(1,npk*nup1);
% Kanyd1=zeros(1,npk*nup1);
% Kanyd2=zeros(1,npk*nup1);
% Kd1e=zeros(1,npk*nup1);
% Kiim=zeros(npk*nup1,npk*nup1);
% Kijm=zeros(npk*nup1,npk*nup1);
% Kdm=zeros(1,npk*nup1);
% Kd1m=zeros(1,npk*nup1);
%
% Kdz=zeros(1,npk*nup1);
% Kiiz=zeros(npk*nup1,npk*nup1);




 if ifalso>0

   bes=real(besselj(mbv1,KK*alon));


 ip=0;
 for imu=pimu
 ipom=(ceil((imu)/pasnu));
 ik1v=Pusas0([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;
 mu=mbv(imu);
  if mu==0
   imu0=1;
  end
 jmu=imu+1;
 ja=jmu+1;
%    for ip1=1:npk
%     ip=ip+1;
    for ip1=ik1v
%     ip=ip1+(ipom-1)*npk;
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
%     for ip2=ip1+1:npk
     for ip2=ip1+1:ik1v(end)
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
   ipom=(ceil(imu/pasnu));
   ik1v=Pusas0([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;
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
%     for ip1=1:npk
%      ip=ip+1;
     for ip1=ik1v
%      ip=ip1+(ipom-1)*npk;
      ip=ip+1;
      Q=KK(ip1)*alon;
      P1=bes(ip1,jmu)^2;
      P2=bes(ip1,jmu-1)*bes(ip1,jmu+1);
      Kdia=.5*(P1-P2);
      Kdz(ip)=Kdia*fbe(3);
%      for ip2=ip1+1:npk
      for ip2=ip1+1:ik1v(end)
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

'dopo ifalso', keyboard


%  disp(' pausa Kij0')
%  pausak


% Kpp=zeros(npk*nup1,npk*nup1);
 Kpp=Kset;
 Kpm=Kpp;
 Ksp=Kpp;
 Ksm=Kpp;
 KAz1=Kpp;
 KBz1=Kpp;
%' qui ret', keyboard

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
%   pmu=pmat+npk*(imu-1)/pasnu;
   for iarm=1:farm
    ipom=(ceil((imu)/pasnu));
    ipom1=(ceil((farm*(imu)+iarm-1)/pasnu));
    ik1v=Pusas0([ldap(ipom1)+1:ldap(ipom1+1)])-(ipom-1)*npk;
    pmu=[ldap(ipom1)+1:ldap(ipom1+1)];

%   ipom=(ceil(imu/pasnu));
%   pmu=[ldap(ipom)+1:ldap(ipom+1)];
%   ik1v=Pusas0([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;

   for jnu=pimu
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
      S=Fi1(imu,inu,:,:)+Fi2(imu,inu,:,:);
      Sr=reshape(S,npk,npk)/fnor;
      Kpp(pmu,pnu)=Sr(ik1v,ik2v);
%      Kpp(pmu,pnu)=reshape(S,npk,npk)/fnor;
      S=Fi1(imu,inu,:,:)-Fi2(imu,inu,:,:);
      Sr=reshape(S,npk,npk)/fnor;
      Kpm(pmu,pnu)=Sr(ik1v,ik2v);
%      Kpm(pmu,pnu)=reshape(S,npk,npk)/fnor;
      S=Fi3(imu,inu,:,:)+Fi4(imu,inu,:,:);
      Sr=reshape(S,npk,npk)/fnor;
      Ksp(pmu,pnu)=Sr(ik1v,ik2v);
%      Ksp(pmu,pnu)=reshape(S,npk,npk)/fnor;
      S=Fi3(imu,inu,:,:)-Fi4(imu,inu,:,:);
      Sr=reshape(S,npk,npk)/fnor;
      Ksm(pmu,pnu)=Sr(ik1v,ik2v);
%      Ksm(pmu,pnu)=reshape(S,npk,npk)/fnor;
      if iztm==1
       S=Fi5(imu,inu,:,:);
       Sr=reshape(S,npk,npk)/fnor;
       KAz1(pmu,pnu)=Sr(ik1v,ik2v);
%       KAz1(pmu,pnu)=reshape(S,npk,npk)/fnor;
       S=Fi6(imu,inu,:,:);
       Sr=reshape(S,npk,npk)/fnor;
       KBz1(pmu,pnu)=Sr(ik1v,ik2v);
%       KBz1(pmu,pnu)=reshape(S,npk,npk)/fnor;
      end

    end
    end
   end
   end
 end %ifalso

'dopo Kz ifalzo -1', keyboard


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

    matpolp0=ones(size(diae));
    matpolp1=ones(size(diae));
    if imu0==1
%     fimam=1:length(KK);
     fimam=1:ldapu(2);
     matpolp0(fimam,:)=0;
     matpolp0(:,fimam)=0;
%     matpolp0(fimam,fimam)=0;
    end

%'prima di pol', keyboard

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
    diam=sqrt((ZMKK));
    prodz=j*diag(KKv./bev.*diam);
    prodzc=conj(prodz);
    prodz=diag(KKv./bev.*diam);
    prodz=prodz(Pusas0);
    prodzc=prodz;
    KMMpz=0.5*prodz*(Ktoiiz+KAz1+pol*KBz1)*prodzc;
    pol=-1;
    KMMmz=0.5*prodz*(Ktoiiz+KAz1+pol*KBz1)*prodzc;
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
    KMMpz=diag(.5*ZMKK.*prodz)*(Ktoiiz+KAz1+pol*KBz1)*diag(prodz_con);
    pol=-1;
    KMMmz=diag(.5*ZMKK.*prodz)*(Ktoiiz+KAz1+pol*KBz1)*diag(prodz_con);
   end
   KAzpdu=[Kzer Kzer; Kzer KMMpz];
   KAzmdu=[Kzer Kzer; Kzer KMMmz];
  else
   KAzpdu=0;
   KAzmdu=0;
  end

  if icirc0==1
   lxivet=ndis;
   ' qui ver icirc0', keyboard
   matv_ci2
   KApdu=KApdu+KApsub/alon^2;
   KAmdu=KAmdu+KAmsub/alon^2;
   KAzpdu=KAzpdu+KAzpsub/alon^2;
   KAzmdu=KAzmdu+KAzmsub/alon^2;
  end


 Kospd (ndis,:,:)=KApdu*alon^2;
 Kosmd (ndis,:,:)=KAmdu*alon^2;
 Koszpd(ndis,:,:)=KAzpdu*alon^2;
 Koszmd(ndis,:,:)=KAzmdu*alon^2;
 disp (' qui uno '), keyboard

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
%    if ~exist('Kospd')  	
%     Kospd=(KAp);
%    end
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

  clear Kospd Kosmd Koszpd Koszmd Fi* A B AB

 if iany>=2
  Kosanp=Kaap/2;
  Kosanm=Kaam/2;
 end

% disp(' qui  xiT  Vett'), keyboard
 %&&&&&&&&&&&&&&&

  if igainshape==1
   isu=[1:lxi];
   if exist('xiT')
    isuT=[1:length(xiT)];
   end
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
% ' Tipo 2 ', keyboard

  if exist('KA')
%   KAp=KA;
%   KAm=KA;

   KApsub=KAp;
   KAmsub=KAm;
  end




else %sha==1

% ' prima civ', keyboard
 lxivet=[1:lxi];
 matv_ci2
% ' dopoi civ', keyboard
% Kosp=Kospsub;
% Kosm=Kosmsub;
% Koszp=Koszpsub;
% Koszm=Koszmsub;
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
si2=2*si;
Znor=[ZEv; ZMv];

%ifp=2
if ifp>=-3, disp(' Kmat_ve '),
 if ifp>1,
  keyboard,
 end,
end


' end kmat_ve2 ', keyboard