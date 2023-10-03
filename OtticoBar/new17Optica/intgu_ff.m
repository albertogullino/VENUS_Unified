% Kff=linspace(0,max(KK),nrff)';
%
% Kff=linspace(0,k_ret,nrff)';
% Kff=linspace(0,0.1,nrff)';

 nup1=1;
 npkF=length(KK);
 Kiie=zeros(nrff,npkF*nup1);
 Kde=zeros(1,npkF*nup1);

 aFFvero=max(xvero)*2/3;
% aFFvero=max(xvero);

 %aFFvero=max([axtot aytot])*1.25;
 cce1=cce;
 if cce==1000
  cce1=0;
 end
aFFvero=(max(max(aytot))+abs(cce1))*1.5;
aFFvero=(max(max(aytot)));
% aFFvero=8;
aFF=aFFvero*kcav0;
%ifp=-10
% ' aFF', keyboard
 bes=real(besselJ(mbv1,KK*aFF));
 besf=real(besselJ(mbv1,Kff*aFF));

 if iLP==0
  imuv=[pimu(1)-1:pimu(end)+1];
 else
  imuv=pimu;
 end


 for imu=imuv
 jmu=imu+1;
    for ip1=1:nrff
     Q=Kff(ip1)*aFF;

     for ip2=1:npkF
      Q1=KK(ip2)*aFF;
      if abs(Q-Q1)<=1e-5
        P1=besf(ip1,jmu)^2;
        P2=besf(ip1,jmu+1)*besf(ip1,jmu-1);
        Kadia=0.5*(P1-P2);
      else
        Q1=KK(ip2)*aFF;
        P1=besf(ip1,jmu)*bes(ip2,jmu-1);
        P2=bes(ip2,jmu)*besf(ip1,jmu-1);
        Kadia=1/(Q^2-Q1^2)*(Q1*P1-Q*P2);
      end
      Kiie(ip1,ip2)=Kadia;
     end
    end

    Im(:,:,imu)=Kiie*aFFvero^2;

 end %mu


