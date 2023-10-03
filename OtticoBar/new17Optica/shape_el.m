maab=max([a b]);
miab=min([a b]);

ng=fix(max([a b])*alim*3*3);
ng1=fix(max([a b])*alim*2*3);
ng2=fix(max([a b])*alim*1*3);

%  [r1,wi]=gauleg(miab,maab,ng);
  [r1,wi]=gauleg(b,a,ng);

  wi=abs(wi);

 rv=r1;
 lrv=length(rv);
 rvp=[r1 fliplr(r1) r1 fliplr(r1) r1(1)]/alon;
 fi=atan(sqrt(((r1/b).^2-1)./(1-(r1/a).^2)));
% fi=atan(sqrt(((r1/maab).^2-1)./(1-(r1/miab).^2)));
 fiv=[fi pi-fliplr(fi) fi+pi 2*pi-fliplr(fi) fi(1)];

 if ifp>=1
  figure, polar(fiv,rvp*alon/kcav0), title(' shape-el'), pausak
  figure, plot(fiv,rvp-1), pausak
 end
 if ifp==-10
  figure, polar(fiv,rvp*alon/kcav0), title(' shape-el'),
  drawnow
 end

 fa1=atan(sqrt(((r1/maab).^2-1)./(1-(r1/miab).^2)));

 if b>a
  sfa=0;
 else
  sfa=pi/2;
 end





 sgimp=ones(size(rv));
% fis=find(rv/a<1);
% sgimp(fis)=-1;

 A=zeros(lrv,nubes+1,nubes+1);
 B=A;

 r=rv/alon;

ivec=0;

if ivec==1

 for imu=pimu
  jmu=imu-meun;
  mu=mbv(imu);
  for inu=pimu
   jnu=inu-meun;
   nu=mbv(inu);
   if (nu+mu)/2-fix((nu+mu)/2)==0
    dmn=mu-nu;
    if dmn==0
     mfat=2*abs(fa1);
%     mfat=2*abs(sfa-fa1);
    else
%     mfat=(sin(dmn*sfa)+sin(dmn*(sfa+pi)))/dmn+(sin(dmn*fa1)+sin(dmn*(fa1+pi)))/dmn;
%     mfat=(sin(dmn*sfa)+sin(dmn*(sfa+pi)))/dmn+(sin(dmn*fa1)+sin(dmn*(fa1+pi)))/dmn;
     mfat1=sin(dmn*(sfa+fa1))/dmn-sin(dmn*sfa)/dmn;
     mfat2=sin(dmn*(fa1+sfa+pi))/dmn-sin(dmn*(sfa+pi))/dmn;
     mfat=mfat1+mfat2;
    end
    mfatd=mfat;
    A(:,jmu,jnu)=mfatd(1:lrv);
    dmn=mu+nu;
    if dmn==0
     mfat=2*abs(fa1);
%     mfat=2*abs(sfa-fa1);
    else
%     mfat=(sin(dmn*sfa)+sin(dmn*(sfa+pi)))/dmn+(sin(dmn*fa1)+sin(dmn*(fa1+pi)))/dmn;
     mfat1=sin(dmn*(sfa+fa1))/dmn-sin(dmn*sfa)/dmn;
     mfat2=sin(dmn*(fa1+sfa+pi))/dmn-sin(dmn*(sfa+pi))/dmn;
     mfat=mfat1+mfat2;
    end
    mfats=mfat;
    B(:,jmu,jnu)=mfats(1:lrv);
    if ifp>10
     figure, plot(r,mfatd,r,mfats)
     title(' A e B in shape-el ')
     disp(' mu, nu='), [mu nu]
     pausak
    end
   end  %if
  end
%  close all
 end

    if ifp>1
     sA=size(A);
     figure, plot(r,reshape(A,sA(1),sA(2)*sA(3))), hold on,
     plot(r,reshape(B,sA(1),sA(2)*sA(3)),'--')
     title(' A e B in shape-el ')
     pausak
    end

% solo per display
if ifp>10
 figure
 for imu=pimu
  jmu=imu-meun;
  mu=mbv(imu);
  for inu=pimu
   jnu=inu-meun;
   nu=mbv(inu);
   if (nu+mu)/2-fix((nu+mu)/2)==0
    mfatd=A(:,jmu,jnu);
    mfats=B(:,jmu,jnu);
     plot(r,mfatd,r,mfats)
     title(' A e B in shape-el ')
     disp(' mu, nu='), [mu nu]
     pausak
   end  %if
  end
 end
% close all
end

end  %ivec


muv=[0:2:2*nubesu];
AB=ones(length(r),length(muv))*NaN;
im=0;
for mu=muv
im=im+1;
   if mu==0
    Ad=2*abs(fa1);
   else
    Ad=2*(sin(mu*(fa1+sfa)))/mu;
   end
   AB(:,im)=Ad';
 fiA=find(abs(AB(:,im))<1e-7);
 if length(fiA)>ng/2
  AB(:,im)=AB(:,im)*0;
 end
end
if ifp>1
 if length(muv)>6
  figure, plot(r,AB(:,1:6)), hold on, plot(r,AB(:,7:length(muv)),'.-'),
 else
  figure, plot(r,AB),
 end
 title(' y0, m2, c4, r6, g8, b10, y12, m14, c16, r18, g20, b22')
 pausak
end

clear Ap Bp
 for imu=pimu
  jmu=imu-meun;
  mu=mbv(imu);
  for inu=pimu
   jnu=inu-meun;
   nu=mbv(inu);
   if (nu+mu)/2-fix((nu+mu)/2)==0
    dmn=abs(mu-nu);
    fim=find(dmn==muv);
    if length(fim)==1
     mfatd=AB(:,fim);
     A(:,jmu,jnu)=AB(:,fim);
    else
     disp('errore A mu in sha_oxi ')
     pausak
    end
    dmn=abs(mu+nu);
    fim=find(dmn==muv);
    if length(fim)==1
     mfats=AB(:,fim);
     B(:,jmu,jnu)=AB(:,fim);
    else
     disp('errore B mu in sha_oxi ')
     pausak
    end
   end  %if
  end
 end
    if ifp>1
     sA=size(A);
     figure, plot(r,reshape(A,sA(1),sA(2)*sA(3))), hold on,
     plot(r,reshape(B,sA(1),sA(2)*sA(3)),'--')
     title(' Ap e Bp in shape-el ')
     pausak
    end
