 nup1=fix(nubes/pasnu)+1;
 Kiie=zeros(npk*nup1,npk*nup1);
 Kzer=Kiie;
 Kije=zeros(npk*nup1,npk*nup1);
 Kde=zeros(1,npk*nup1);
 Kiim=zeros(npk*nup1,npk*nup1);
 Kijm=zeros(npk*nup1,npk*nup1);
 Kdm=zeros(1,npk*nup1);
 Kd1m=zeros(1,npk*nup1);

 aFFvero=max(xvero);
 aFF=aFFvero*kcav0;
   bes=real(besselj(mbv1,KK*aFF));


 for imu=[pimu(1)-1 pimu pimu(end)+1]
 jmu=imu+1;
 ip=0;
    for ip1=1:npk
     ip=ip+1;
     Q=KK(ip1)*aFF;
     P1=bes(ip1,jmu)^2;
     P2=bes(ip1,jmu+1)*bes(ip1,jmu-1);
     Kdia=0.5*(P1-P2);
     Kde(ip)=Kdia;

     for ip2=ip1+1:npk
      Q1=KK(ip2)*aFF;
      P1=bes(ip1,jmu)*bes(ip2,jmu-1);
      P2=bes(ip2,jmu)*bes(ip1,jmu-1);
      Kadia=1/(Q^2-Q1^2)*(Q1*P1-Q*P2);
      Kiie(ip,ip+ip2-ip1)=Kadia;
     end
    end

    Im(:,:,imu)=(Kiie+Kiie'+diag(Kde))*aFFvero^2;

 end %mu

%    Ktoiie=(Kiie+Kiie'+diag(Kde))/4;
%    Ktoije=-(Kije+Kije'+diag(Kd1e))/4;
%    Ktoiim=(Kiim+Kiim'+diag(Kdm))/4;
%    Ktoijm=-(Kijm+Kijm'+diag(Kd1m))/4;

