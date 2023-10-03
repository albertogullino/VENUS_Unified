iFF=1;
%z0c=10;
te0g=90;
Anf=AnFF;
if iFF==1
%      z0c=input(' distanza ff in cm = ');
      z0=z0c*1e4;
      k0=2*pi/lambda;
      teRg=asin(KKt*rr)*180/pi;
      fik=find(teRg>te0g);
      if length(fik)>0
       fit=fik;
       Anf(fit)=0;
       figure, plot(abs(Anout)), pausak
      end
      teR=asin(KK(1:length(KK))*rr);
      rox=z0*tan(teR);
      R=z0;
%      F=cos(teR)/R.*exp(-j*k0*R+j*pi/4);
      F=1./(R.*tan(teR)).*exp(-j*k0*R+j*pi/4);

      teRd=asin(KKt*rr);
      Rd=z0;
      Fd=1./(R.*tan(teRd))*exp(-j*k0*Rd+j*pi/4);
      lK=length(KK);
      lKv=diff(ldap);
      ldap1=ldap;
%      ldap1(end)=ldap1(end)-1;
      ldap1=ldap1+1;
      Pus1=Pus;
      Pus1=[Pus1 Pus1(end)+1];
      roxt=[0; rox];
      s2teR=sin(teR).^2;

    if iLP==0

      Efx=0;
      Efy=0;
      Efz=0;
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
       zer=zeros(nk1max-length(AnEd),1);
       fkti=[zer; fkt(pumi+lKAn/2)];
       AnE=[zer; AnEd];
       AnM=[zer; AnMd];

%       cos1=exp(i*mmp*pi/2);
%       cos2=exp(i*mmm*pi/2);
%       cos3=exp(i*(mmp-1)*pi/2);
       xJ=k0*R*s2teR;
       cos1=cos(xJ+pi/2*(mmp-0.5));
       cos2=cos(xJ+pi/2*(mmm-0.5));
       cos3=cos(xJ+pi/2*(mmp-1.5));
       xJ=0;
       cos1=exp(j*(xJ+pi/2*(mmp-0.5)));
       cos2=exp(j*(xJ+pi/2*(mmm-0.5)));
       cos3=exp(j*(xJ+pi/2*(mmp-1.5)));

       Efx=Efx+((AnE-AnM).*cos1)*fmu(mmp,:)+((AnE+AnM).*cos2)*fmu(mmm,:);
       Efy=Efy+(-(AnE-AnM).*cos1)*gmu(mmp,:)+((AnE+AnM).*cos2)*gmu(mmm,:);
%       Efy=Efy+(-AnE+AnM)*gmu(mmp,:)*cos1+(AnE+AnM)*gmu(mmm,:)*cos2;
       Efz=Efz+j*fatEz*(AnM.*fkti.*cos3)*fmu(mmp-1,:);
      end

      Efx=[Efx(1,:); Efx];
      Efy=[Efy(1,:); Efy];
      Efz=[Efz(1,:); Efz];
      Ef=abs(Efx).^2+abs(Efy).^2+abs(Efz).^2;
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

