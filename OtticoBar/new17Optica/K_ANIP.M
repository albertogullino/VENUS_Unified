clear Ktoian Kiie Kde

 nup1=fix(nubes/pasnu)+1;
% ' nup1 ', keyboard

if ~exist('Pus0')
 Pus0=Pus;
end
 lPU=length(Pus0);
 lPUi=length(Pus0)/nup1;
 Kset=zeros(lPU,lPU);
 Kseti=zeros(lPUi,lPUi);
 Kiie=Kseti;

 alon=P.Ry*kcav;

   bes=real(besselj(mbv1,KK*alon));


 ip=0;
 icm=0;
 for imu=pimu
 icm=icm+1;
 ipom=(ceil((imu)/pasnu));
 ik1v=Pusas0([ldap(ipom)+1:ldap(ipom+1)])-(ipom-1)*npk;
 mu=mbv1(imu);
 jmu=imu;
    for ip1=ik1v
     ip=ip1;
     Q=KK(ip1)*alon;
     P1=bes(ip1,jmu)^2;
     P2=bes(ip1,jmu+1)*bes(ip1,jmu-1);
     Kdia=.5*(P1-P2);
     Kde(ip)=Kdia;
     for ip2=ip1+1:ik1v(end)
      Q1=KK(ip2)*alon;
      P1=bes(ip2,jmu)*bes(ip1,jmu-1);
      P2=bes(ip1,jmu)*bes(ip2,jmu-1);
      Kadia=1/(Q1^2-Q^2)*(Q*P1-Q1*P2);
      Kiie(ip,ip2)=Kadia;
     end
    end
    Ktoian{icm}=(Kiie+Kiie'+diag(Kde))/4;
%     'vedi'
%     pausak
 end %mu

  if isi==0
   diae=diag(ZEKK(ik1v));
   diam=diag(ZMKK(ik1v));
   diae1=1/2;
   diam1=1/2;
  else
   diae=diag(sqrt(ZEKK(ik1v)));
   diam=diag(sqrt(ZMKK(ik1v)));
   diae(1,1)=diae(1,1)*fatde1;
   diam(1,1)=diam(1,1)*fatde1;
   diae1=diae/2;
   diam1=diam/2;
  end

 icm=0;
 KaEE=Kset;
 KaEM=Kset;
 KaME=Kset;
 KaMM=Kset;

nupd=is_even(nube);
if nupd==0

 pumo=1:nup1-1;
 if length(pumo)==0
  icm=1;
  Kana=Ktoian{icm};
  KEEp=diae*Kana*diae1;
  KEMp=diae*Kana*diam1;
  KMEp=diam*Kana*diae1;
  KMMp=diam*Kana*diam1;
   KaEE0=KEEp;
   KaEM0=KEMp;
   KaME0=KMEp;
   KaMM0=KMMp;
 end

 for icm=pumo
  Kana=Ktoian{icm};
  KEEp=diae*Kana*diae1;
  KEMp=diae*Kana*diam1;
  KMEp=diam*Kana*diae1;
  KMMp=diam*Kana*diam1;
  KaEE((icm-1)*lPUi+ik1v,icm*lPUi+ik1v)=KEEp;
  KaEM((icm-1)*lPUi+ik1v,icm*lPUi+ik1v)=KEMp;
  KaME((icm-1)*lPUi+ik1v,icm*lPUi+ik1v)=-KMEp;
  KaMM((icm-1)*lPUi+ik1v,icm*lPUi+ik1v)=-KMMp;
  if icm==1
   KaEE0=KEEp;
   KaEM0=KEMp;
   KaME0=KMEp;
   KaMM0=KMMp;
  end
   KaEE((icm)*lPUi+ik1v,(icm-1)*lPUi+ik1v)=KEEp;
   KaEM((icm)*lPUi+ik1v,(icm-1)*lPUi+ik1v)=-KEMp;
   KaME((icm)*lPUi+ik1v,(icm-1)*lPUi+ik1v)=KMEp;
   KaMM((icm)*lPUi+ik1v,(icm-1)*lPUi+ik1v)=-KMMp;
 end

 KEEd=KaEE;
 KEEd(ik1v,ik1v)=KaEE0;
 KEMd=KaEM;
 KEMd(ik1v,ik1v)=KaEM0;
 KMEd=KaME;
 KMEd(ik1v,ik1v)=KaME0;
 KMMd=KaMM;
 KMMd(ik1v,ik1v)=KaMM0;
 Kanpp=[KEEd segem*KEMd; segem*KMEd KMMd]*alon^2;

 KEEd=KaEE;
 KEEd(ik1v,ik1v)=-KaEE0;
 KEMd=KaEM;
 KEMd(ik1v,ik1v)=-KaEM0;
 KMEd=KaME;
 KMEd(ik1v,ik1v)=-KaME0;
 KMMd=KaMM;
 KMMd(ik1v,ik1v)=-KaMM0;
 Kanpm=[KEEd segem*KEMd; segem*KMEd KMMd]*alon^2;

else

  icm=1;
  Kana=Ktoian{icm};
  KEEp=diae*Kana*diae1;
  KEMp=diae*Kana*diam1;
  KMEp=diam*Kana*diae1;
  KMMp=diam*Kana*diam1;
  KaEE0=KEEp;
  KaEM0=KEMp;
  KaME0=KMEp;
  KaMM0=KMMp;

 for icm=2:nup1-1
  Kana=Ktoian{icm};
  KEEp=diae*Kana*diae1;
  KEMp=diae*Kana*diam1;
  KMEp=diam*Kana*diae1;
  KMMp=diam*Kana*diam1;
  KaEE((icm-1)*lPUi+ik1v,icm*lPUi+ik1v)=KEEp;
  KaEM((icm-1)*lPUi+ik1v,icm*lPUi+ik1v)=KEMp;
  KaME((icm-1)*lPUi+ik1v,icm*lPUi+ik1v)=-KMEp;
  KaMM((icm-1)*lPUi+ik1v,icm*lPUi+ik1v)=-KMMp;
  KaEE((icm)*lPUi+ik1v,(icm-1)*lPUi+ik1v)=KEEp;
  KaEM((icm)*lPUi+ik1v,(icm-1)*lPUi+ik1v)=-KEMp;
  KaME((icm)*lPUi+ik1v,(icm-1)*lPUi+ik1v)=KMEp;
  KaMM((icm)*lPUi+ik1v,(icm-1)*lPUi+ik1v)=-KMMp;
 end

 KEEd=KaEE;
 KEMd=KaEM;
 KEMd(ik1v+lPUi,ik1v)=-2*KaEM0;
 KMEd=KaME;
 KMEd(ik1v,ik1v+lPUi)=-KaME0;
 KMMd=KaMM;
 KMMd(ik1v,ik1v+lPUi)=-KaMM0;
 KMMd(ik1v+lPUi,ik1v)=-2*KaMM0;
 Kanpp=[KEEd segem*KEMd; segem*KMEd KMMd]*alon^2;

 KEEd=KaEE;
 KEEd(ik1v,ik1v+lPUi)=KaEE0;
 KEEd(ik1v+lPUi,ik1v)=2*KaEE0;
 KEMd=KaEM;
 KEMd(ik1v,ik1v+lPUi)=KaEM0;
 KMEd=KaME;
 KMEd(ik1v+lPUi,ik1v)=2*KaME0;
 KMMd=KaMM;

 Kanpm=[KEEd segem*KEMd; segem*KMEd KMMd]*alon^2;



end


if ifp==-10
' fine K_anis', keyboard
end