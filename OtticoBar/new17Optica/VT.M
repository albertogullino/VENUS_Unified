load t
igau=0;
if length(T)~=1

if ifp==-10
' discretizzo Temp '
keyboard
end
sT=size(T);
 if length(sT)==2

 ' qui dentro', keyboard
 T0=T(1,:);
 rodisT=rdis;
 for iz=1:length(zdis)-1
  fi=find(zeta>=zdis(iz) & zeta<zdis(iz+1));
  Tz(iz)=mean(T0(fi));
  mea=mean(T(:,fi)');
%  iz
%  pausak
  for ir=1:length(rdis)-1
   fi=find(ro_inT>=rdis(ir) & ro_inT<rdis(ir+1));
   Tr(ir)=mean(mea(fi));
  end
  Tr(ir+1)=0;
  Tdis(:,iz)=Tr';
 end
 fi=find(isnan(Tdis)==1);
 Tdis(fi)=0;


 zp=[zdis(1:ndz) zdis(2:ndz+1)-1e-4];
 tp=[Tz Tz];
 [zpu,iso]=sort(zp);
 tpu=tp(iso);
 if ifp==-10
 figure, plot(zpu,tpu,zeta,T0)
 figure, plot(rdis,Tdis,'.-')
 keyboard
 end
 elseif length(sT)==3
  ' OK' , keyboard
 end

 igau=4;
else
 rodisT=0;
 zdis=0;
 Tdis=0;
end

