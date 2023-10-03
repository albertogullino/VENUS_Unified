clear
 load
 T0=T(1,:);

 rodisT=rdis;
 for iz=1:length(zdis)-1
% for iz=1:length(zdis)
  fi=find(zeta>=zdis(iz) & zeta<zdis(iz+1));
  Tz(iz)=mean(T0(fi));
  mea=mean(T(:,fi)');
  for ir=1:length(rdis)-1
   fi=find(ro_in>=rdis(ir) & ro_in<rdis(ir+1));
   Tr(ir)=mean(mea(fi));
  end
%  Tr(ir+1)=0;
  Tdis(:,iz)=Tr';
 end
 fi=find(isnan(Tdis)==1);
 Tdis(fi)=0;
 rdis=rdis(1:end-1);

 figure, plot(zpu,tpu,zeta,T0)

 dT=max(T0)-min(T0);
 sdT=dT/(length(zdis)-1);

icd=1;
clear Tr
t0=T0(1);
ki=1;
 for k=2:length(T0)
  if abs(T0(k)-t0)>=sdT
   Td(icd)=mean(T0(ki:k));
   zd(icd)=zeta(k);
   t0=T0(k);
  mea=mean(T(:,ki:k)');
  for ir=1:length(rdis)-1
   fi=find(ro_in>=rdis(ir) & ro_in<rdis(ir+1));
   Tr(ir)=mean(mea(fi));
  end
  Tr(ir+1)=mea(end);
  Tdis1(:,icd)=Tr';
  figure, plot(rdis,Tdis1), pausak
   ki=k;
   icd=icd+1;
  end
 end
 Td(end+1)=mean(T0(ki:end));

 zp=[[0 zd(1:end)] [zd(1:end)-1e-4 zeta(end) ]];
 tp=[Td Td];
 [zpu1,iso]=sort(zp);
 tpu1=tp(iso);

 figure, plot(zpu,tpu,zeta,T0,zpu1,tpu1)

