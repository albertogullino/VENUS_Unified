ng=fix((b)*alim*3*3);
ng1=fix((b)*alim*2*3);
ng2=fix((b)*alim*1*3);

 if igint==1
  [r1,wi]=gauleg(a,b,ng);
%  wi=-wi;
 else
  abi=b-a;
  r11=linspace(b,b-abi/10,ng);
  r12=linspace(b-abi/10,a+abi/10,2*ng);
  r13=linspace(a+abi/10,a+1e-6,ng);
  r1=[r11(1:10) r12 r13(2:11)];
 end
 rv=r1;
 lrv=length(rv);
 rvp=[r1 fliplr(r1) r1 fliplr(r1)]/a;
 fi=atan(sqrt(((r1/b).^2-1)./(1-(r1/a).^2)));
 fiv=[fi pi-fliplr(fi) fi+pi 2*pi-fliplr(fi)];
 if ifp>=1
  figure, polar(fiv,rvp), title(' shape-el'), pausak
  figure, plot(fiv,rvp-1), pausak
 end

% fa1=fliplr(fi);
% rv=fliplr(rv);
 fa1=(fi);

 A=zeros(lrv,nubes+1,nubes+1);
 B=A;

 for imu=2:pasnu:length(mbv)-1
  jmu=imu-1;
  mu=mbv(imu);
  for inu=2:pasnu:length(mbv)-1
   jnu=inu-1;
   nu=mbv(inu);
   if (nu+mu)/2-fix((nu+mu)/2)==0
    dmn=mu-nu;
    if dmn==0
     mfat=2*fa1;
    else
     mfat=(sin(dmn*fa1)+sin(dmn*(fa1+pi)))/dmn;
    end
    mfatd=mfat;
    A(:,jmu,jnu)=mfatd(1:lrv);
    dmn=mu+nu;
    if dmn==0
     mfat=2*fa1;
    else
     mfat=(sin(dmn*fa1)+sin(dmn*(fa1+pi)))/dmn;
    end
    mfats=mfat;
    B(:,jmu,jnu)=mfats(1:lrv);
    R=rv/a;
    if ifp>10
     figure, plot(R,mfatd,R,mfats)
     title(' A e B in shape-el ')
     disp(' mu, nu='), [mu nu]
     pausak
    end
   end  %if
  end
  close all
 end

    if ifp>1
     sA=size(A);
     figure, plot(R,reshape(A,sA(1),sA(2)*sA(3))), hold on,
     plot(R,reshape(B,sA(1),sA(2)*sA(3)),'--')
     title(' A e B in shape-el ')
     pausak
    end

% solo per disply
if ifp>1
 figure
 for imu=2:pasnu:length(mbv)-1
  jmu=imu-1;
  mu=mbv(imu);
  for inu=2:pasnu:length(mbv)-1
   jnu=inu-1;
   nu=mbv(inu);
   if (nu+mu)/2-fix((nu+mu)/2)==0
    mfatd=A(:,jmu,jnu);
    mfats=B(:,jmu,jnu);
     plot(R,mfatd,R,mfats)
     title(' A e B in shape-el ')
     disp(' mu, nu='), [mu nu]
     pausak
   end  %if
  end
 end
 close all
end
