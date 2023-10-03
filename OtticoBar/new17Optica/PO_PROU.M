' po_prou passo', keyboard
nk1max=length(KK);
if iredmat==1
 Pusff=Pusas;
else
 Pusff=Pus;
end
te0g=90;
Aep=Afz(:,iFFsect)+Afzf(:,iFFsect);
Ahp=Afz(:,iFFsect)-Afzf(:,iFFsect);
if iFFsect==1
 segdir=-1;
else
 segdir=1;
end

%Anf=AnFF;
%Anh=AnFFH;

if iant==0
 Anf=Aep.*pes(Pus);
 Anh=Ahp.*pes(Pus);
else
 Anf=Aep;
 Anh=Ahp;
end


fktt=KKt./sqrt(1-KKt.^2);
fkttn=KKt./sqrt(1-(rrfa*KKt).^2);
if ipolar==0
 fktt=[fktt; fktt];
 fkttn=[fkttn; fkttn];
end

fktt=fktt(Pus);
fkttn=fkttn(Pus);

fkn=sqrt(1-(rrfa*KKv).^2);
Znorn=[1./fkn; fkn];
ZnornE=[1./fkn; 1./fkn];
ZnornM=[fkn; fkn];

%Znorn=Znorn(Pus);
%ZnornE=ZnornE(Pus);
%ZnornM=ZnornM(Pus);

%' qui FF', keyboard
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

      Cfir=diag(cos(fian));
      Sfir=diag(sin(fian));

      Efxv=0;
      Efyv=0;
      Efz=0;
      Hfx=0;
      Hfy=0;
      Hfz=0;

      EfxE=0;
      EfyE=0;
      EfxM=0;
      EfyM=0;
      HfxE=0;
      HfyE=0;
      HfxM=0;
      HfyM=0;


      Anfz=-segdir*segem*2*j*Anh.*fktt*fatEz;
      AnfH=rrfa*Anh./Znor(Pus);
      AnfHz=segdir*2*j*rrfa*Anf./Znor(Pus).*fktt;

      AnfHn=Anf./Znorn(Pus);
      AnfHnE=Anf./ZnornE(Pus);
      AnfHnM=Anf./ZnornM(Pus);

%' far field ', keyboard
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
%       figure, plot(abs(AnE),'r.'); pausak
%       figure, plot(abs(AnM),'r.'); pausak




%       AnE=[zeros(nk1max-length(AnEd),1); AnEd];
%       AnM=[zeros(nk1max-length(AnEd),1); AnMd];

       AnMd=Anfz(pumi+lKAn/2);
       AZ=zeros(nk1max,1);
%       AZ=[zeros(nk1max-length(AnEd),1); AnMd];
       AZ(Pudu(Pusff(pumi)))=AnMd;


       Imm=Fd*reshape(Im(:,:,mmm),npkf,npkF);
       Imp=Fd*reshape(Im(:,:,mmp),npkf,npkF);
       Imi=Fd*reshape(Im(:,:,mmi),npkf,npkF);

       AD=(AnE-segem*AnM);
       AS=(AnE+segem*AnM);

       Efxv=Efxv+((Imp*AD.*j^nup)*fmu(mmp,:)+(Imm*AS.*j^num)*fmu(mmm,:));
       Efyv=Efyv+((-Imp*AD.*j^nup)*gmu(mmp,:)+(Imm*AS.*j^num)*gmu(mmm,:));

       EfxE=EfxE+((Imp*AnE.*j^nup)*fmu(mmp,:)+(Imm*AnE.*j^num)*fmu(mmm,:));
       EfyE=EfyE+((-Imp*AnE.*j^nup)*gmu(mmp,:)+(Imm*AnE.*j^num)*gmu(mmm,:));
       EfxM=EfxM+segem*((-Imp*AnM.*j^nup)*fmu(mmp,:)+(Imm*AnM.*j^num)*fmu(mmm,:));
       EfyM=EfyM+segem*((Imp*AnM.*j^nup)*gmu(mmp,:)+(Imm*AnM.*j^num)*gmu(mmm,:));

       Efz=Efz+(Imi*AZ)*j^nui*fmu(mmi,:);

       AnEd=AnfH(pumi);
       AnMd=AnfH(pumi+lKAn/2);

%       AnE=[zeros(nk1max-length(AnEd),1); AnEd];
%       AnM=[zeros(nk1max-length(AnEd),1); AnMd];
       AnE=zeros(nk1max,1);
       AnE(Pudu(Pusff(pumi)))=AnEd;
       AnM=zeros(nk1max,1);
       AnM(Pudu(Pusff(pumi)))=AnMd;

       AnMd=AnfHz(pumi);
%       AZ=[zeros(nk1max-length(AnMd),1); AnMd];
       AZ=zeros(nk1max,1);
       AZ(Pudu(Pusff(pumi)))=AnMd;

       AD=(AnE-segem*AnM);
       AS=(AnE+segem*AnM);

       Hfy=Hfy+((Imp*AD.*j^nup)*fmu(mmp,:)+(Imm*AS.*j^num)*fmu(mmm,:));
       Hfx=Hfx-((-Imp*AD.*j^nup)*gmu(mmp,:)+(Imm*AS.*j^num)*gmu(mmm,:));

       Hfz=Hfz+(Imi*AZ)*j^nui*gmu(mmi,:);

       AnEd=AnfHnE(pumi);
       AnMd=AnfHnE(pumi+lKAn/2);
       AnE=[zeros(nk1max-length(AnEd),1); AnEd];
       AnM=[zeros(nk1max-length(AnEd),1); AnMd];

       AD=AnE;
       AS=AnE;

       HfyE=HfyE+((Imp*AD.*j^nup)*fmu(mmp,:)+(Imm*AS.*j^num)*fmu(mmm,:));
       HfxE=HfxE-((-Imp*AD.*j^nup)*gmu(mmp,:)+(Imm*AS.*j^num)*gmu(mmm,:));

       AnEd=AnfHnM(pumi);
       AnMd=AnfHnM(pumi+lKAn/2);
%       AnE=[zeros(nk1max-length(AnEd),1); AnEd];
%       AnM=[zeros(nk1max-length(AnEd),1); AnMd];
       AnE=zeros(nk1max,1);
       AnE(Pudu(Pusff(pumi)))=AnEd;
       AnM=zeros(nk1max,1);
       AnM(Pudu(Pusff(pumi)))=AnMd;


       AD=-segem*AnM;
       AS=+segem*AnM;

       HfyM=HfyM+((Imp*AD.*j^nup)*fmu(mmp,:)+(Imm*AS.*j^num)*fmu(mmm,:));
       HfxM=HfxM-((-Imp*AD.*j^nup)*gmu(mmp,:)+(Imm*AS.*j^num)*gmu(mmm,:));

      end


      Efzn=-Fdz*((EfxM+EfxE)*Cfir+(EfyM+EfyE)*Sfir);
      Hfzn=Fdz*(-(HfxM)*Sfir+(HfyM)*Cfir);

%      Efx=EfxE+EfxM;
%      Efy=EfyE+EfyM;
%nuove
%      Hfx=HfxE+HfxM;
%      Hfy=HfyE+HfyM;

      Efx=Efxv;
      Efy=Efyv;

      Efxn=Efx;
      Efyn=Efy;


      Hfxn=(-(HfxM-HfxE)*(Cfir.*Sfir)-HfyE*Cfir.^2-HfyM*Sfir.^2);
      Hfyn=((HfyM-HfyE)*(Cfir.*Sfir)+HfxM*Cfir.^2+HfxE*Sfir.^2);


      Ef=(Efx).^2+(Efy).^2+(Efz).^2;
      mfxy=(max(max(abs(Ef))));
      mfxyc=sqrt(max(max((Ef))));
      Ef=sqrt(abs(Ef/mfxy));
%
% figure;
% pograp=[500   50   380   880];
% set(gcf,'Position',pograp)
% subplot(3,1,1)
%            surf(X,Y,abs(Efz/mfxyc)),
%            shading('interp'), view(0,90),
%            colorbar
%            title(' Ef_z ')
%            axis square, axis equal, grid,
%            axis([-1 1 -1 1]*axli/2),
% subplot(3,1,2)
%            surf(X,Y,abs(Efx/mfxyc)),
%            shading('interp'), view(0,90),
%            colorbar
%            title(' Ef_x ')
%            axis square, axis equal, grid,
%            axis([-1 1 -1 1]*axli/2),
% subplot(3,1,3)
%            surf(X,Y,abs(Efy/mfxyc)),
%            shading('interp'), view(0,90),
%            colorbar
%            title(' Ef_y ')
%            axis square, axis equal, grid,
%            axis([-1 1 -1 1]*axli/2),

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



%      Ef=(Efx).^2;
%      mfxy=(max(max(abs(Ef))));
%      mfxyc=sqrt(max(max((Ef))));
%      Ef=sqrt(abs(Ef/mfxy));
%     'po_pro'
%     keyboard



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

%            figure
%            surf(X,Y,fipla0), colorbar
%            shading('interp'), view(0,90),

%      Pvectx=real(Efy.*conj(Hfz)-Efz.*conj(Hfy)).*cos(fipla).*sin(fiazi);
%      Pvecty=real(Efx.*conj(Hfz)-Efz.*conj(Hfx)).*sin(fipla).*sin(fiazi);
%      Pvectz=real(Efx.*conj(Hfy)-Efy.*conj(Hfx)).*cos(fiazi);

%      Pvectx=real(Efy.*conj(Hfz)-Efz.*conj(Hfy))/sqrt(rrfa);
%      Pvecty=-(real(Efx.*conj(Hfz)-Efz.*conj(Hfx)))/sqrt(rrfa);
%      Pvectz=real(Efx.*conj(Hfy)-Efy.*conj(Hfx))/rrfa;

     if iLP==0
      Pvectx=real(Efy.*conj(Hfz)-Efz.*conj(Hfy));
      Pvecty=-(real(Efx.*conj(Hfz)-Efz.*conj(Hfx)));
      Pvectz=real(Efx.*conj(Hfy)-Efy.*conj(Hfx));

      Pvectxn=real(Efyn.*conj(Hfzn)-Efzn.*conj(Hfyn));
      Pvectyn=-(real(Efxn.*conj(Hfzn)-Efzn.*conj(Hfxn)));
      Pvectzn=real(Efxn.*conj(Hfyn)-Efyn.*conj(Hfxn));
     else
%      Pvectz=(real(Efx.*conj(Hfy)));
      '    Attenzione: cambio in po_pro '
      Pvectz=abs(Efx.*conj(Hfy));
     end


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

%'end popro'
%keyboard
