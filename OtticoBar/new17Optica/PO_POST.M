

%close all
iploloc=0;
%' po_prou passo', keyboard
nk1max=length(KK);
if iredmat==1
 Pusff=Pusas;
else
 Pusff=Pus;
end
te0g=90;
%if iFFsect==1
% segdir=-1;
%else
% segdir=1;
%end

 segdir=1;

Anf=Are0;
Anh=Arh0;


fkttn=KKt./sqrt(1-(rrfa*KKt).^2);
if ipolar==0
 fkttn=[fkttn; fkttn];
end


lim=1-(rrfa*KKv).^2;
kacce=find(lim<0);
kaccea=[kacce; kacce+length(lim)];

fkttn=fkttn(Pusff);

if iLP==0
SQ=sqrt(1-(KKv).^2);
Znorn=[1./SQ; SQ];
Ynorn=1./Znorn;
%Ynorn(kaccea)=0;
Ynorn=Ynorn(Pusff);

fkttn=KKt./[SQ; SQ];
%fkttn(kaccea)=0;
fkttn=fkttn(Pusff);
 Anh=Anh.*Ynorn;
 Anfz=-segdir*segem*2*(-j)*Anf.*fkttn*fatEz;
 AnfHz=2*(-j)*Anh.*fkttn;
else
Ynorn=1;
fkttn=KKt;
end

fatcut=1;
%fatcut=ones(size(KKt));
%fatcut(kaccea)=0;
%fatcut=fatcut(Pusff);



if iFF==1
%      z0c=input(' distanza ff in cm = ');
      z0=z0c*1e4;
      lambdaff=lambda/rff;
      k0=2*pi/lambdaff;
      teRg=asin(Kff*rrfa)*180/pi;
      fik=find(teRg>te0g);
      if length(fik)>0
       fit=fik;
       Anf(fit)=0;
       figure, plot(abs(Anout)), pausak
      end
     if iFFte==0
% vecchio ff
      rox=z0c*tan(teR);
      R=z0./cos(teR);
%      Fdu=-1./(R.*tan(teR)).*exp(-j*k0*R+j*pi/4);
      Fd0=cos(teR)./(2*pi*lambdaff*R).*exp(-j*k0*R+j*pi/4);
      Fd=diag(Fd0);
     else
% nuovo ff
      Rd=z0;
      Fd0=1/(2*pi*lambdaff*Rd)*exp(-j*k0*Rd+j*pi/4);
      Fd=Fd0*diag(cos(teR));
     end
      Fdz=diag(tan(teR));

      lK=length(KK);
      lKv=diff(ldap);
      ldap1=ldap;
      ldap1=ldap;
%      ldap1(end)=ldap1(end)-1;
      ldap1=ldap1+1;
      Pus1=Pus;
      Pus1=[Pus1 Pus1(end)+1];

    if iLP==0

      Efx=0;
      Efy=0;
      Efz=0;
      Hfx=0;
      Hfy=0;
      Hfz=0;


      Anf=Anf.*fatcut;

%      Anfz=-segdir*segem*2*(-j)*Anh.*fkttn*fatEz;

%      AnfH=Anh.*Ynorn;
      AnfH=Anh;

      AnfHz=2*(-j)*Anf.*Ynorn.*fkttn;


      lkr=1;
      if exist('krv')
        lkr=length(krv);
      end
      lkr1=1;
      Puf=1:nk1max;
      Pudu=repmat(Puf,1,numodi)';
%      for imi=1:length(ibm)*lkr



      for imi=1:(lbv-1)/pasnu*lkr
%      for imi=1:(lbv-1)/pasnu
       imiver=ceil(imi/lkr);
%       [imiver], pausak
%       imi1=imi-1;
%       mm=mbv(ibz(imi));
%       mmp=ibp(imi);
%       mmm=ibm(imi);
%       mmi=ibp(imi)-1;
       imi1=imiver-1;
       mmloc=mbv(ibz(imiver));
       mmp=ibp(imiver);
       mmm=ibm(imiver);
       mmi=ibp(imiver)-1;
       nui=mmloc;
       num=mmloc-1;
       nup=mmloc+1;
%       lK=lKv(imi);
%       pumi=ldap(imi)+1:ldap(imi+1);
       lK=lKv(imiver);
       pumi=ldap(imiver)+1:ldap(imiver+1);
       pumif=Pus1(ldap1(imi)):Pus1(ldap1(imi+1)-1);

       AnEd=Anf(pumi);
       AnMd=Anf(pumi+lKAn/2);

       AnE=zeros(nk1max,1);
       AnE(Pudu(Pusff(pumi)))=AnEd;
       AnM=zeros(nk1max,1);
       AnM(Pudu(Pusff(pumi)))=AnMd;
%       figure, plot(abs(AnM),'r.'); pausak




%       AnE=[zeros(nk1max-length(AnEd),1); AnEd];
%       AnM=[zeros(nk1max-length(AnEd),1); AnMd];

       AnMd=Anfz(pumi+lKAn/2);

       AZ=zeros(nk1max,1);
       AZ(Pudu(Pusff(pumi)))=AnMd;
     if iploloc==1
      figure, plot(abs(AnE),'r.'), hold on, plot(abs(AnM),'g.'),
      plot(abs(AZ),'c.'),
      title(' Campo elettrico E M Z'), pausak
     end

       Imm=Fd*reshape(Im(:,:,mmm),npkf,npkF);
       Imp=Fd*reshape(Im(:,:,mmp),npkf,npkF);
       Imi=Fd*reshape(Im(:,:,mmi),npkf,npkF);

       AD=(AnE-segem*AnM);
       AS=(AnE+segem*AnM);

       Efx=Efx+((Imp*AD.*j^nup)*fmu(mmp,:)+(Imm*AS.*j^num)*fmu(mmm,:));
       Efy=Efy+((-Imp*AD.*j^nup)*gmu(mmp,:)+(Imm*AS.*j^num)*gmu(mmm,:));

       Efz=Efz+(Imi*AZ)*j^nui*fmu(mmi,:);

%        ' po_post ', keyboard


       AnEd=AnfH(pumi);
       AnMd=AnfH(pumi+lKAn/2);
       AnE=zeros(nk1max,1);
       AnE(Pudu(Pusff(pumi)))=AnEd;
       AnM=zeros(nk1max,1);
       AnM(Pudu(Pusff(pumi)))=AnMd;


       AnMd=AnfHz(pumi);
       AZ=zeros(nk1max,1);
       AZ(Pudu(Pusff(pumi)))=AnMd;

     if iploloc==1
      figure, plot(abs(AnE),'r.'), hold on, plot(abs(AnM),'g.'),
      plot(abs(AZ),'c.'),
      title(' Campo magnetico '), pausak
     end

       AD=(AnE-segem*AnM);
       AS=(AnE+segem*AnM);

       Hfy=Hfy+((Imp*AD.*j^nup)*fmu(mmp,:)+(Imm*AS.*j^num)*fmu(mmm,:));
       Hfx=Hfx-((-Imp*AD.*j^nup)*gmu(mmp,:)+(Imm*AS.*j^num)*gmu(mmm,:));

       Hfz=Hfz+(Imi*AZ)*j^nui*gmu(mmi,:);

      end



      Ef=(Efx).^2+(Efy).^2+(Efz).^2;
      mfxy=(max(max(abs(Ef))));
      mfxyc=sqrt(max(max((Ef))));
      Ef=sqrt(abs(Ef/mfxy));

    else   %iLP



      Efx=0;
      Hfy=0;

      AnfH=rrfa*Anh;

      ibm=[1:pasnu:lbv];
      mbvc=mbv(ibm);
      if ipolar==0
       ibm=[ibm ibm];
       mbvc=[mbvc mbvc];
      end


      lPus=0;
%      for imi=1:length(ibm)
      lkr=1;
      if exist('krv')
        lkr=length(krv);
      end
      lkr1=1;
      Puf=1:nk1max;
      Pudu=repmat(Puf,1,numodi)';
      for imi=1:length(ibm)*lkr
%      for imi=1:lkr:length(ibm)*lkr
%      for imi=1:lkr:length(ibm)*lkr
       imiver=ceil(imi/lkr);
       mmloc=mbvc(imiver);

       mmi=ibm(imiver);
       nui=mmloc;
       lK=lKv(imi);
       pumi=ldap(imi)+1:ldap(imi+lkr1);
       pumif=Pus1(ldap1(imi)):Pus1(ldap1(imi+lkr1)-1);

       AnEd=Anf(pumi);
%       figure, plot(abs(Anf)), hold on, plot(pumi,abs(AnEd),'r.'),  pausak
%       AnE=[zeros(nk1max-length(AnEd),1); AnEd];
       AnE=zeros(nk1max,1);
       AnE(Pudu(Pusff(pumi)))=AnEd;
%       figure, plot(abs(AnE),'r.'); pausak

       Imi=Fd*reshape(Im(:,:,mmi),npkf,npkF);

       AD=AnE;

%'Efx', keyboard
       Efx=Efx+(Imi*AD.*j^nui)*fmu(imiver,:);

%       figure, plot(abs(Efx)), pausak


%       AnEd=AnfH(pumi);
%       AnE=[zeros(nk1max-length(AnEd),1); AnEd];

       AnE=zeros(nk1max,1);
       AnE(Pudu(Pusff(pumi)))=AnEd;

       AD=AnE;

       Hfy=Hfy+(Imi*AD.*j^nui)*fmu(imiver,:);
%       ' po_pro', pausak

      end


    end

%'z0'
%keyboard
isal=1;
if isal==0
  fx=find(imag(roxt)==0);
  Xdu=roxt(fx)*fm0*1e-4;
  Ydu=roxt(fx)*gm0*1e-4;
  Ef=Ef(fx,:);

  if iFFte==1
   X=atan(Xdu/z0c)*180/pi;
   Y=atan(Ydu/z0c)*180/pi;
   axli=60;
  else
   X=Xdu;
   Y=Ydu;
   axli=z0c;
  end

   fiazi=atan(sqrt(Xdu.^2+Ydu.^2)/z0c);
   fipla=ones(size(fx))*fian0;
else

  if iFFte==1

   axli=60;
  else
   axli=z0c;
  end

end %isal

% Pvectz=(Ex.*conj(Hy)-Ey.*conj(Hx));
% Pvectx=(Ey.*conj(Hz)-Ez.*conj(Hy));
% Pvecty=-((Ex.*conj(Hz)-Ez.*conj(Hx)));


     if iLP==0
      Pvectx=real(Efy.*conj(Hfz)-Efz.*conj(Hfy));
      Pvecty=-(real(Efx.*conj(Hfz)-Efz.*conj(Hfx)));
      Pvectz=real(Efx.*conj(Hfy)-Efy.*conj(Hfx));
      Pointing=sqrt(Pvectx.^2+Pvecty.^2+Pvectz.^2);

     else
%      Pvectz=(real(Efx.*conj(Hfy)));
      '    Attenzione: cambio in po_pro '
      Pvectz=abs(Efx.*conj(Hfy));
     end
%     Ppo=Pointing/Pointing(1,1);
%     Ppo=Pvectz/Pvectz(1,1);

% figure, plot(X1(:,1),abs(Efx(:,46)/Efx(1,1)).^2,X1(:,1),...
% abs(Ppo(:,46)))
%load sa1
% hold on, plot(X1(:,1),Ec,'.',X1(:,1),Pc,'.')

% Ec=abs(Efx(:,46)/Efx(1,1)).^2;
% Pc=abs(Pvectz(:,46)/Pvectz(1,1));
% save sa Pc Ec

iskf=1;
if iskf==0
  figure
  subplot(1,3,1)
  surf(X,Y,(abs(Efx)/mfxy).^2),
  shading('interp'), view(0,90),
  axis square, axis equal, grid,
  axis([-1 1 -1 1]*axli/2),

  subplot(1,3,2)
  surf(X,Y,(abs(Efy)/mfxy).^2),
  shading('interp'), view(0,90),
  axis square, axis equal, grid,
  axis([-1 1 -1 1]*axli/2),

%  figure
%  surf(Xte,Yte,(abs(Efy)/mfxy).^2),
%  shading('interp'), view(0,90),
%  axis square, axis equal, grid,

  subplot(1,3,3)
  surf(X,Y,(abs(Ef)/mfxy).^2),
  shading('interp'), view(0,90),
  axis square, axis equal, grid,
  axis([-1 1 -1 1]*axli/2),
  xlabel(' x cm '), ylabel(' y  cm'),
  drawnow
end   %iskf
 if ifp>=-3
  pausak
 end
end %iFF

 if iLP==0
 Poim=sqrt(Pvectx.^2+Pvecty.^2+Pvectz.^2);
 else
 Poim=abs(Pvectz);
 end

'end po_post'
%keyboard
