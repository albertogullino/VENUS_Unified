uuFunc=ztot'*0;
uuF0=ztot'*0;
ztom=floor(ztot'*1000*1e4)*1e-4;
dJ=floor(dJ*1e4)*1e-4;
lJ=length(dJ);
for il=1:lJ-1
 if fJ(il+1)==-1
  fiz=find(ztom>dJ(il) & ztom<=dJ(il+1));
  uuFunc(fiz)=1;
%  keyboard
 end
 if fJ(il+1)==-1 & il==leqw
  uuF0(fiz)=1;
 end
end

Fi2=real((Fidu.*indz).^2);
%figure, plot(ztot,uuFunc,ztot,uuF0,'.',ztot,Fi2/max(Fi2)), pausak
%figure, plot(ztot,uuFunc.*Fi2,'.',ztot,uuF0.*Fi2,'ro'), pausak
% figure, plot(hz/1000,Fi/max(Fi),ztot,(Fidu/max(Fidu)).^2)

dz=[diff(ztot); 0];
fN=Fi2*dz;
udz=uuFunc'.*dz;
[du,ima]=max(udz);
fi=find(udz>0);
udz(fi)=udz(fi(1));


uLong0=(Fi2.*uuF0/fN)*dz;
uLong=(Fi2/fN)*udz;

%du=Fi2.*uuF0;
%figure, plot(ztot,du/max(du),'.')
%sua=sum(du);
%du=Fi2.*uuFunc;
%hold on, plot(ztot,du/max(du),'ro')
%sut=sum(du);
%R2=sut/sua


%Fi3=real(relPerm).*Fi';
%du=Fi3.*uF0;
%figure, plot(hz/1000,du/max(du),'.')
%suan=sum(du);
%du=Fi3.*uFunc;
%hold on, plot(hz/1000,du/max(du),'ro')
%sutn=sum(du);
%R3=sutn/suan

uLong=(Fi2/fN)*udz;

uL(1)=uLong0;
uL(2)=uLong;


