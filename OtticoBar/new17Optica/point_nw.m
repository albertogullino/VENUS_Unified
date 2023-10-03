iFF=1;
%z0c=10;
te0g=90;
%Anf=AnFF;
%Anh=AnFFH;

fktt=KKt./sqrt(1-KKt.^2);

fkttn=KKt./sqrt(1-(rr*KKt).^2);

fkn=sqrt(1-(rr*KKv).^2);
Znorn=[1./fkn; fkn];
fktt=fktt(Pusas);
fkttn=fkttn(Pusas);
Znorn=Znorn(Pusas);

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
      Fd0=1/(2*pi*lambda*Rd)*exp(-j*k0*Rd+j*pi/4);
      Fd=Fd0*diag(cos(teR));

      Fdz=diag(tan(teR));

      lK=length(KK);
      lKv=diff(ldap);
      ldap1=ldap;
      ldap1=ldap;
%      ldap1(end)=ldap1(end)-1;
      ldap1=ldap1+1;
      Pus1=Pus;
      Pus1=[Pus1 Pus1(end)+1];
      Cfir=diag(cos(fian));
      Sfir=diag(sin(fian));

    if iLP==0

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


      Anfz=j*Anh.*fktt*fatEz;
      AnfH=Anh./Znor(Pus);
      AnfHz=-j*Anf./Znor(Pus).*fktt;

 %     Anfz=j*Anh.*fkttn*fatEz;
 %     AnfH=Anh./Znorn(Pus);
 %     AnfHz=-j*Anf./Znorn(Pus).*fkttn;

      AnfHn=Anf./Znorn(Pus);

%      AnfHn=Anf./Znor(Pus);

      for imi=1:(lbv-1)/pasnu
       imi1=imi-1;
       mm=mbv(ibz(imi));
       mmp=ibp(imi);
       mmm=ibm(imi);
       mmi=ibp(imi)-1;
       nui=mm;
       num=mm-1;
       nup=mm+1;
       lK=lKv(imi);
       pumi=ldap(imi)+1:ldap(imi+1);
       pumif=Pus1(ldap1(imi)):Pus1(ldap1(imi+1)-1);
       AnEd=Anf(pumi);
       AnMd=Anf(pumi+lKAn/2);
       AnE=[zeros(nk1max-length(AnEd),1); AnEd];
       AnM=[zeros(nk1max-length(AnEd),1); AnMd];

       AnMd=Anfz(pumi+lKAn/2);
       AZ=[zeros(nk1max-length(AnEd),1); AnMd];


       Imm=Fd*reshape(Im(:,:,mmm),npkf,npk);
       Imp=Fd*reshape(Im(:,:,mmp),npkf,npk);
       Imi=Fd*reshape(Im(:,:,mmi),npkf,npk);

       AD=(AnE-AnM);
       AS=(AnE+AnM);

       Efyv=Efyv+((Imp*AD.*j^nup)*fmu(mmp,:)+(Imm*AS.*j^num)*fmu(mmm,:));
       Efxv=Efxv+((-Imp*AD.*j^nup)*gmu(mmp,:)+(Imm*AS.*j^num)*gmu(mmm,:));

       EfxE=EfxE+((Imp*AnE.*j^nup)*fmu(mmp,:)+(Imm*AnE.*j^num)*fmu(mmm,:));
       EfyE=EfyE+((-Imp*AnE.*j^nup)*gmu(mmp,:)+(Imm*AnE.*j^num)*gmu(mmm,:));
       EfxM=EfxM+((-Imp*AnM.*j^nup)*fmu(mmp,:)+(Imm*AnM.*j^num)*fmu(mmm,:));
       EfyM=EfyM+((Imp*AnM.*j^nup)*gmu(mmp,:)+(Imm*AnM.*j^num)*gmu(mmm,:));

       Efz=Efz+(Imi*AZ)*j^nui*fmu(mmi,:);


       AnEd=AnfH(pumi);
       AnMd=AnfH(pumi+lKAn/2);
       AnE=[zeros(nk1max-length(AnEd),1); AnEd];
       AnM=[zeros(nk1max-length(AnEd),1); AnMd];

       AnMd=AnfHz(pumi);
       AZ=[zeros(nk1max-length(AnMd),1); AnMd];


       AD=(AnE-AnM);
       AS=(AnE+AnM);

       Hfy=Hfy+((Imp*AD.*j^nup)*fmu(mmp,:)+(Imm*AS.*j^num)*fmu(mmm,:));
       Hfx=-Hfx-((-Imp*AD.*j^nup)*gmu(mmp,:)+(Imm*AS.*j^num)*gmu(mmm,:));
       Hfz=Hfz+(Imi*AZ)*j^nui*gmu(mmi,:);

       AnEd=AnfHn(pumi);
       AnMd=AnfHn(pumi+lKAn/2);
       AnE=[zeros(nk1max-length(AnEd),1); AnEd];
       AnM=[zeros(nk1max-length(AnEd),1); AnMd];
       HfxE=HfxE+((Imp*AnE.*j^nup)*fmu(mmp,:)+(Imm*AnE.*j^num)*fmu(mmm,:));
       HfyE=HfyE+((-Imp*AnE.*j^nup)*gmu(mmp,:)+(Imm*AnE.*j^num)*gmu(mmm,:));
       HfxM=HfxM+((-Imp*AnM.*j^nup)*fmu(mmp,:)+(Imm*AnM.*j^num)*fmu(mmm,:));
       HfyM=HfyM+((Imp*AnM.*j^nup)*gmu(mmp,:)+(Imm*AnM.*j^num)*gmu(mmm,:));



      end

       Efxn1=((EfxM-EfxE)*(Cfir.*Sfir)+EfxM*Cfir.^2+EfxE*Sfir.^2);
       Efyn1=((EfyM-EfyE)*(Cfir.*Sfir)+EfyE*Cfir.^2+EfyM*Sfir.^2);

       Efzn=-Fdz*(EfxM*Cfir+EfyM*Sfir);
       Hfxn=(-(HfxM-HfxE)*(Cfir.*Sfir)-HfyE*Cfir.^2-HfyM*Sfir.^2);
       Hfyn=((HfyM-HfyE)*(Cfir.*Sfir)+HfxM*Cfir.^2+HfxE*Sfir.^2);
       Hfzn=Fdz*(-HfxE*Sfir+HfyE*Cfir);


      Efx=EfxE+EfxM;
      Efy=EfyE+EfyM;

      Efxn=Efx;
      Efyn=Efy;
      Efxn=Efxn1;
      Efyn=Efyn1;

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

      Pvectxn=real(Efyn.*conj(Hfzn)-Efzn.*conj(Hfyn));
      Pvectyn=-(real(Efxn.*conj(Hfzn)-Efzn.*conj(Hfxn)));
      Pvectzn=real(Efxn.*conj(Hfyn)-Efyn.*conj(Hfxn));




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



%  zte=
%
