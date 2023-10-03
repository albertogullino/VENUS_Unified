
ng=fix((b)*alim*3*3);
ng1=fix((b)*alim*2*3);
ng2=fix((b)*alim*1*3);
 if b~=a
  if igint==1
   [rv1,wi1]=gauleg(a,b,ng1);
   [rv2,wi2]=gauleg(b,sqrt(a^2+b^2),ng2);
   rv=[rv1 rv2];
   wi=[wi1 wi2];
  else
   nrv1=5*ng;
   nrv2=fix(nrv1*(sqrt(b^2+a^2)-b)/(b-a));
   rv1=linspace(a,b,nrv1);
   rv1=rv1(1:nrv1-1);
   rv2=linspace(b,sqrt(a^2+b^2),nrv2);
   rv=[rv1 rv2];
  end

  r1=rv2;
  fa1=atan(sqrt((r1/b).^2-1));
  r21=rv1;
  fa21=pi/2-atan(sqrt((r21/a).^2-1));
  r22=rv2;
  fa22=pi/2-atan(sqrt((r22/a).^2-1));
  fa2=[fa21 fa22];
  r2=[r21 r22];
  lrv=length(rv);
  fat=[fa1 fliplr(fa2)];
  rut=[r1 fliplr(r2)];
 else

  if igint==1
   [rv2,wi]=gauleg(b,sqrt(a^2+b^2),ng);
   rv=[rv2];
  else
   nrv1=50;
   rv2=linspace(b,sqrt(a^2+b^2),nrv1);
   rv=[rv2];
  end
  r1=rv2;
  fa1=atan(sqrt((r1/b).^2-1));
  fa2=pi/2-atan(sqrt((r1/b).^2-1));
  lrv=length(rv);
  fat=[fa1 fliplr(fa2)];
  rut=[r1 fliplr(r1)];
 end

 if ifp>=1
  figure, polar(fat,rut/a), title(' shape-sq'), pausak
 end


 if ifp>=1 & (ifr==1)
  figure, plot(rut/a,fat)
  pausak
  figure, polar(fat,rut/a)
  pausak
 end
 A=zeros(lrv,nubes+1,nubes+1);
 B=A;
% C=A;
% D=A;

 for imu=pimu
  jmu=imu-meun;
  mu=mbv(imu);
  for inu=pimu
   jnu=inu-meun;
   nu=mbv(inu);
   if (nu+mu)/2-fix((nu+mu)/2)==0
    dmn=mu-nu;
    if dmn==0
     if b~=a
      mfa1=2*fa21;
      mfa2=2*(fa22-fa1);
      mfat1=[mfa1 mfa2];
      fif=find(abs(mfat1)<1e-10);
      mfat1(fif)=0;
     else
      mfat1=[2*(fa2-fa1)];
     end
%      mfat2=zeros(1,length(mfat1));
    else
     if b~=a
      mfa1=(sin(dmn*fa21)+sin(dmn*(fa21+pi)))/dmn;
      mfa2=(sin(dmn*fa22)-sin(dmn*fa1)+sin(dmn*(fa22+pi))-sin(dmn*(fa1+pi)))/dmn;
      mfat1=[mfa1 mfa2];
%      mfa1=(cos(dmn*fa21)+cos(dmn*(fa21+pi)))/dmn;
%      mfa2=(cos(dmn*fa22)-cos(dmn*fa1)+cos(dmn*(fa22+pi))-cos(dmn*(fa1+pi)))/dmn;
%      mfat2=-[mfa1 mfa2];
     else
      mfat1=(sin(dmn*fa2)-sin(dmn*fa1)+sin(dmn*(fa2+pi))-sin(dmn*(fa1+pi)))/dmn;
%      mfat2=-(cos(dmn*fa2)-cos(dmn*fa1)+cos(dmn*(fa2+pi))-cos(dmn*(fa1+pi)))/dmn;
     end
    end
    mfd1=mfat1;
%    mfd2=mfat2;
    A(:,jmu,jnu)=mfat1(1:lrv);
%    C(:,jmu,jnu)=mfat2(1:lrv);
    dmn=mu+nu;
    if dmn==0
     if b~=a
      mfa1=2*fa21;
      mfa2=2*(fa22-fa1);
      mfat=[mfa1 mfa2];
      fif=find(abs(mfat)<1e-10);
      mfat(fif)=0;
     else
      mfat=[2*(fa2-fa1)];
     end
     mfat1=mfat;
%     mfat2=zeros(1,length(mfat1));
    else
     if b~=a
      mfa1=(sin(dmn*fa21)+sin(dmn*(fa21+pi)))/dmn;
      mfa2=(sin(dmn*fa22)-sin(dmn*fa1)+sin(dmn*(fa22+pi))-sin(dmn*(fa1+pi)))/dmn;
      mfat1=[mfa1 mfa2];
%      mfa11=(cos(dmn*fa21)+cos(dmn*(fa21+pi)))/dmn;
%      mfa21=(cos(dmn*fa22)-cos(dmn*fa1)+cos(dmn*(fa22+pi))-cos(dmn*(fa1+pi)))/dmn;
%      mfa10=(cos(dmn*fa21*0)+cos(dmn*(fa21*0+pi)))/dmn;
%      mfa20=(cos(dmn*fa22*0)-cos(dmn*fa1*0)+cos(dmn*(fa22*0+pi))-cos(dmn*(fa1*0+pi)))/dmn;
%      mfa1=mfa11-mfa10;
%      mfa2=mfa21-mfa20;
%      mfat2=-[mfa1 mfa2];
     else
      mfat1=(sin(dmn*fa2)-sin(dmn*fa1)+sin(dmn*(fa2+pi))-sin(dmn*(fa1+pi)))/dmn;
%      mfat21=(cos(dmn*fa2)-cos(dmn*fa1)+cos(dmn*(fa2+pi))-cos(dmn*(fa1+pi)))/dmn;
%      mfat20=(cos(dmn*fa2*0)-cos(dmn*fa1*0)+cos(dmn*(fa2*0+pi))-cos(dmn*(fa1*0+pi)))/dmn;
%      mfat2=-(mfat21-mfat10);
     end
    end
    B(:,jmu,jnu)=mfat1(1:lrv);
%    D(:,jmu,jnu)=mfat2(1:lrv);
    mfs1=mfat1;
%    mfs2=mfat2;
    R=rv/a;
    if ifp>1
%     figure, plot(R,mfd1,R,mfs1,R,mfd2,R,mfs2)
     figure, plot(R,mfd1,R,mfs1)
     title(' A e B ')
     disp(' mu, nu='), [mu nu]
     pausak
    end
   end  %if
  end
%  close
 end
%keyboard
