iFF=1;
z0c=10;
te0g=40;
if iFF==1
%      z0c=input(' distanza ff in cm = ');
      z0=z0c*1e4;
      k0=2*pi/lambda;
      teRg=asin(KKt*rr)*180/pi;
      fik=find(teRg>te0g);
      Anf=Anout;
      if length(fik)>0
       fit=fik;
       Anout(fit)=Anout(fit)*0;
       figure, plot(abs(Anout)), pausak
      end
      teR=asin(KK(2:length(KK))*rr);
      rox=z0*tan(teR);
      R=z0./cos(teR);
      F=-1./(R.*tan(teR)).*exp(-j*k0*R+j*pi/4);
      lK=length(KK);
      Efx=0;
      Efy=0;
      for imi=1:(lbv-1)/pasnu
       imi1=imi-1;
       mm=mbv(ibz(imi));
       mmp=ibp(imi);
       mmm=ibm(imi);
       AnE=F.*Anf((2:lK)+imi1*lK);
       AnM=F.*Anf((2:lK)+(imi1+nup1)*lK);
%%%%%%  cos1=cos((mm+3/2)*pi/2);
%%%%%%  cos2=cos((mm-1/2)*pi/2);
       cos1=exp(i*(mm+3/2)*pi/2);
       cos2=exp(i*(mm-1/2)*pi/2);
       Efx=Efx+(AnE-AnM)*fmu(mmp,:)*cos1+(AnE+AnM)*fmu(mmm,:)*cos2;
       Efy=Efy+(-AnE+AnM)*gmu(mmp,:)*cos1+(AnE+AnM)*gmu(mmm,:)*cos2;
%%%%%%  figure, plot(abs(AnE)), hold on, plot(abs(AnM),'r'), pausak
      end

      Efx=[Efx(1,:); Efx];
      Efy=[Efx(1,:); Efy];
      roxt=[1; rox];

%%%%%% Ef=sqrt(abs(Efx).^2+abs(Efy).^2+abs(Ez).^2);
      Ef=sqrt(abs(Efx).^2+abs(Efy).^2);
      mfxy=(max(max(Ef)));

  fx=find(imag(roxt)==0);
  X=roxt(fx)*fm0*1e-4;
  Y=roxt(fx)*gm0*1e-4;

%  subplot(2,3,4)
  figure
  surf(X,Y,(abs(Efx(fx,:))/mfxy).^2),
  shading('interp'), view(0,90),
  axis square, axis equal, grid,
  axis([-1 1 -1 1]*z0c),
  xlabel(' x cm '), ylabel(' y  cm'),

















  drawnow
%  if ifp>=-3
%   pausak
%  end
end %iFF


