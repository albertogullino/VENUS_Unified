
dz=diff(hz(1:2));
dFp=diff(Fi);
dFs=diff(dFp);
zd=dFp(1:end-1).*dFp(2:end);
fmi=find(zd<0 & dFs>0);
fma=find(zd<0 & dFs<0);
' parte commentata per perdite '
%fma(1)=L_i(1);
f0=find(diff([fma; 1])~=0);
fma=fma(f0);
Fip=-Fi.*imag(perm)';
Fic=Fi;
 for k=1:length(fma)-1
  fiz=find(hz>hz(fma(k)) & hz<=hz(fma(k+1)) );
  me=sum(Fip(fiz))/diff(hz(fiz([1 end])))*dz;
  fmepu(k)=me;
  zmeu(k)=mean(hz(fiz));
  mec=sum(Fic(fiz))/diff(hz(fiz([1 end])))*dz;
  fmepuc(k)=mec;
%  figure, plot(hz(fiz),Fip(fiz),mean(hz(fiz)),me,'ro'), pausak
 end

zmu=[fma(1) zmeu];
fmepu=[0 fmepu];
ffi=interp1(zmu,fmepu,hz,'linear');

fmepuc=[0 fmepuc];
ffic=interp1(zmu,fmepuc,hz,'linear');

%if ifp>0 | ifp==-10
if ifp>0
 figure, plot(hz,-imag(perm),'w',hz,Fip,zmu,fmepu,'g',hz,ffi,'r'), pausak
 figure, plot(hz,Fic,zmu,fmepuc,'g',hz,ffic,'r'), pausak
end
Iper=sum(ffi(find(isnan(ffi)==0)))/1000;
zper=([fma(1) zmeu]-fma(1))/1000;
Prof_p=fmepu/Iper;

Iperc=sum(ffic(find(isnan(ffi)==0)))/1000;
Prof_c=fmepuc/Iperc;
dglo.zp=zper;
dglo.fp=Prof_p;
dglo.fc=Prof_c;

%'subperd', keyboard
