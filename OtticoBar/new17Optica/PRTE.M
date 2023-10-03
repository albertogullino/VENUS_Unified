load
 T=T';
 T0=T(1,:);
 rodisT=rdis;
 for iz=1:length(zdis)-1
  fi=find(zeta>=zdis(iz) & zeta<zdis(iz+1));
  Tz(iz)=mean(T0(fi));
  mea=mean(T(:,fi)');
  for ir=1:length(rdis)-1
   fi=find(ro_in>=rdis(ir) & ro_in<rdis(ir+1));
   Tr(ir)=mean(mea(fi));
  end
  Tr(ir+1)=0;
  Tdis(:,iz)=Tr';
 end


 zp=[zdis(1:ndz) zdis(2:ndz+1)-1e-4];
 tp=[Tz Tz];
 [zpu,iso]=sort(zp);
 tpu=tp(iso);

 figure, plot(zpu,tpu,zeta,T0)
 figure, plot(rdis,Tdis,'.-')
