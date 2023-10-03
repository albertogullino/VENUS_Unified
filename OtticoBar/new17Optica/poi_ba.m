iFF=1;
%z0c=10;
te0g=90;
Anf=AnFF;
Anh=AnFFH;

 teRg=repmat(teR,numodi,1);
 Sqff=sqrt((k0/rr)^2-KKv.^2);
 fktff=KKv./Sqff;
 Zte=1./Sqff;
 Ztm=Sqff;
 Znorff=[Zte; Ztm];


if iFF==1
%      z0c=input(' distanza ff in cm = ');
      z0=z0c*1e4;
      k0=2*pi/lambda;
%      teRg=asin(KKtff*rr)*180/pi;
      fik=find(teRg>te0g);
      if length(fik)>0
       fit=fik;
       Anf(fit)=0;
       figure, plot(abs(Anout)), pausak
      end
      Rd=z0;
      Fd=cos(teRd)/(2*pi*lambda*Rd)*exp(-j*k0*Rd+j*pi/4);

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
      Anfz=j*Anh.*fktff*fatEz;
      AnfH=Anh./Znorff(Pus);
      AnfHz=-j*Anf./Znorff(Pus).*fktff;
      for imi=1:(lbv-1)/pasnu
       imi1=imi-1;
       mm=mbv(ibz(imi));
       mmp=ibp(imi);
       mmm=ibm(imi);
       lK=lKv(imi);
       pumi=ldap(imi)+1:ldap(imi+1);
       pumif=Pus1(ldap1(imi)):Pus1(ldap1(imi+1)-1);
       AnEd=Fd(pumif).*Anf(pumi);
       AnMd=Fd(pumif+lKA/2).*Anf(pumi+lKAn/2);
       AnE=[zeros(nk1max-length(AnEd),1); AnEd];
       AnM=[zeros(nk1max-length(AnEd),1); AnMd];

       zer=zeros(nk1max-length(AnEd),1);
       fktiff=[zer; fktff(pumi+lKAn/2)];

       AnMd=Fd(pumif+lKA/2).*Anfz(pumi+lKAn/2);
       AZ=[zeros(nk1max-length(AnEd),1); AnMd];

       Imm=reshape(Im(:,:,mmm),npkf,npk);
       Imp=reshape(Im(:,:,mmp),npkf,npk);
       Imi=reshape(Im(:,:,mmi),npkf,npk);

       AD=(AnE-AnM);
       AS=(AnE+AnM);

       Efx=Efx+((Imp*AD.*j^mmp)*fmu(mmp,:)+(Imm*AS.*j^mmm)*fmu(mmm,:));
       Efy=Efy+((-Imp*AD.*j^mmp)*gmu(mmp,:)+(Imm*AS.*j^mmm)*gmu(mmm,:));
       Efz=Efz+(Imi*(AZ.*fktiff).*j^mmi)*fmu(mmi,:);



%       cos1=exp(i*mmp*pi/2);
%       cos2=exp(i*mmm*pi/2);
%       cos3=exp(i*(mmp-1)*pi/2);
%       Efx=Efx+(AnE-AnM)*fmu(mmp,:)*cos1+(AnE+AnM)*fmu(mmm,:)*cos2;
%       Efy=Efy+(-AnE+AnM)*gmu(mmp,:)*cos1+(AnE+AnM)*gmu(mmm,:)*cos2;
%       Efz=Efz+AnMz*fmu(mmp-1,:)*cos3;

       AnEd=Fd(pumif).*AnfH(pumi);
       AnMd=Fd(pumif+lKA/2).*AnfH(pumi+lKAn/2);
       AnE=[zeros(nk1max-length(AnEd),1); AnEd];
       AnM=[zeros(nk1max-length(AnEd),1); AnMd];

       AnMd=Fd(pumif+lKA/2).*AnfHz(pumi+lKAn/2);
       AZ=[zeros(nk1max-length(AnMd),1); AnMd];

       AD=(AnE-AnM);
       AS=(AnE+AnM);

       Hfy=Hfy+((Imp*AD.*j^mmp)*fmu(mmp,:)+(Imm*AS.*j^mmm)*fmu(mmm,:));
       Hfx=-Hfx-((-Imp*AD.*j^mmp)*gmu(mmp,:)+(Imm*AS.*j^mmm)*gmu(mmm,:));
       Hfz=Hfz+(Imi*(AZ.*fktiff).*j^mmi)*gmu(mmi,:);


%       Hfy=Hfy+(AnE-AnM)*fmu(mmp,:)*cos1+(AnE+AnM)*fmu(mmm,:)*cos2;
%       Hfx=Hfx+(-AnE+AnM)*gmu(mmp,:)*cos1+(AnE+AnM)*gmu(mmm,:)*cos2;
%       Hfz=Hfz+AnMz*gmu(mmp-1,:)*cos3;
      end

%      Efx=[Efx(1,:); Efx];
%      Efy=[Efy(1,:); Efy];
%      Efz=[Efz(1,:); Efz];

%      Hfx=-[Hfx(1,:); Hfx];
%      Hfy=[Hfy(1,:); Hfy];
%      Hfz=[Hfz(1,:); Hfz];


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

    ibm=[1:pasnu:lbv];
    nup1=length(ibm);;
    mbvc=mbv(ibm);

      Ef=0;
      lPus=0;
      for im=1:length(ibm)
       mm=mbvc(im);
       imi1=im-1;
       pumi=ldap(im)+1:ldap(im+1);
       pumif=Pus1(ldap1(im)):Pus1(ldap1(im+1)-1);

% figure, plot([1:length(Anf)],abs(Anf),'.-',pumi,abs(Anf(pumi)),'o')
% pausak
%!%       pumi=ldap(im)+1:ldap(im+1)
%!       pumid=lPus+(1:(ldap(im+1)-ldap(im)));
%!%       pumid=ldap(im)+1:ldap(im+1);
%!       pumi=pumid
%!       lPus=lPus+length(pumid);
%!figure, plot([1:length(Anf)],abs(Anf),'.-',pumi,abs(Anf(pumi)),'o')
%!       pausak
%!%       pumif=Pus(pumi)-Pus1(ldap1(im)+1)+1;
%!%       pumi=(1:lK)+imi1*lK;
       AnEd=Fd(pumif).*Anf(pumi);
       AnE=[zeros(nk1max-length(AnEd),1); AnEd];
       cos1=exp(i*mm*pi/2);
       Ef=Ef+AnE*fmu(im,:)*cos1;
      end
      mfxy=(max(max(abs(Ef))));
      Ef=[Ef(1,:); Ef];
      Ef=abs(Ef)/mfxy;
    end


  fx=find(imag(roxt)==0);
  Xdu=roxt(fx)*fm0*1e-4;
  Ydu=roxt(fx)*gm0*1e-4;
  Ef=Ef(fx,:);
iFFte=1;
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
%   fize=find(Xdu==0);
%   Xdu1=Xdu;
%   Xdu1(fize)=1e-10;
%   fipla=unwrap(atan((Ydu./Xdu1)));
   fipla=ones(size(fx))*fian0;

%            figure
%            surf(X,Y,fipla0), colorbar
%            shading('interp'), view(0,90),

%      Pvectx=real(Efy.*conj(Hfz)-Efz.*conj(Hfy)).*cos(fipla).*sin(fiazi);
%      Pvecty=real(Efx.*conj(Hfz)-Efz.*conj(Hfx)).*sin(fipla).*sin(fiazi);
%      Pvectz=real(Efx.*conj(Hfy)-Efy.*conj(Hfx)).*cos(fiazi);

      Pvectx=real(Efy.*conj(Hfz)-Efz.*conj(Hfy));
      Pvecty=-(real(Efx.*conj(Hfz)-Efz.*conj(Hfx)));
      Pvectz=real(Efx.*conj(Hfy)-Efy.*conj(Hfx));




%            figure
%            surf(X,Y,fiazi), colorbar
%            shading('interp'), view(0,90),
%            figure
%            surf(X,Y,cos(fipla)), colorbar
%            shading('interp'), view(0,90),


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

%plo_poin


