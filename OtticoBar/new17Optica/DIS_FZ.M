
 Ares=(Acoz(pu1,:)+Acoz(pu2,:));
 dz=.5;
 Zma=200;
% dz=5;
% Zma=1000;
 if iFFsect==1
  if ~exist('dz')
   dz=diff(zi(1:2));
  end
  sectu=[0:dz:Zma];
  sectu=-fliplr(sectu(2:end));
  ztot=[sectu zi];
 else
  if ~exist('dz')
   dz=diff(zi(end-1:end));
  end
  sectu=[0:dz:Zma];
  sectu=sectu(2:end);
  ztot=[zi zi(end)+sectu];
 end
 lambdam=lambda*(1+ze);
 if iFFsect==1
  rff=nv(1,1);
 else
  rff=nv(end,1);
 end
 beta0=2*pi/lambdam*rff;
 sqr=sqrt(1-(KKv(Pusas)*rr/rff).^2);
 Au=Ares(:,1);
 for ks=1:length(sectu)
  Lu=sectu(ks);
  Az=(Au.*exp(-j*beta0*abs(Lu)*sqr));
  Azz(:,ks)=Az;
 end
 Azu=[Azz Ares];
%   AnFs=Andu.*Fs;
%   Mvefm=diag(AnFs)*Mvefm0;
%   Exdu=besm'*Mvefm;
 if exist('if90')
  bes_an=diag(Mvefm0(:,if90))*besm;
 else
  bes_an=besm;
 end
 Esr=(bes_an'*Azu);
 Es=abs(Esr/max(max(Esr))).^2;
 filen=find(shavet==7);
 if length(filen)>1
  filen=1:filen(end);
  sz=sum(Litn(filen))+.5;
 else
  sz=0.5;
 end
% fiz=find(ztot<sz);
 fiz=1:length(ztot);
 [du,pu1]=min(abs(xvero-aiat/2));
 [du,pu2]=min(abs(xvero-aiat));

 fiz0=find(ztot==0);
 fiz0=fiz0(1);
% figure
 if exist('soloplot')
  if soloplot==1
   return
  end
 end

 figure, plot(ztot(fiz),Es([1 pu1 pu2],fiz)),
 xlabel(' zeta toward device (micron)')
 ylabel(' Es ')
 title(' in dis-fz: 3 radial coord: 0, Ox/2, Ox ')
 pausak

 fiz0=find(ztot==0);
 fiz0=fiz0(1);

 fiex=find(ztot<0);
 [du,fima]=max(Es(1,fiex));
 Enor=diag(1./Es(1,fiz0));
 E0=Es(:,fiz0)*Enor;
 puma=fiex(fima);
 Em=Es(:,puma)*Enor;
 figure, h=plot(xvero,[E0 Em]);
 set(h,'linewidth',5);
 title(' thick: cuts at lens (yellow) and at focus (magenta) ')
 st_z=10;
 zf=st_z;
 Ecut=[];
 while zf<Zma
  [du,fizf]=min(abs(ztot+zf));
%  fizf
  zf=zf+st_z;
  Em=Es(:,fizf)*Enor;
  Ecut=[Ecut Em];
 end
 hold on
 plot(xvero,Ecut),
 a=axis;
 a(2)=1.5*aiat;
 axis(a)
 pausak

isk=1;
if isk==0
 map(Es,ztot,xvero)
 axis([ztot(1),ztot(end),0,10])
 pausak
 [du,pu1]=min(abs(xvero-aiat/2));
 [du,pu2]=min(abs(xvero-aiat));
 map(Es(:,fiz),ztot(fiz),xvero)
 axis([ztot(fiz(1)),ztot(fiz(end)),0,10])
 pausak
end

 %map(Es(:,fiz)*diag(1./Es(1,fiz)),ztot(fiz),xvero)

 map(Es(:,fiz)*diag(1./max(Es(:,fiz))),ztot(fiz),xvero)

 axis([ztot((1)),ztot((end)),0,10])
 hold on,
 xi1=xvero;
 xi0=xvero;
 zi1=ones(size(xvero))*0;
 zi0=ones(size(xvero));
 for k=1:3
  xi1=[xi1; xi0];
  zi1=[zi1; zi0*sum(Litn(1:k))];
 end
 hi=ones(size(xi1));
 plot3(zi1',xi1',hi','w');

 pausak

% fiz0=find(ztot==0);
% figure
% plot(xvero,Es(:,fiz0)*diag(1./max(Es(:,fiz0))))
% pausak

 return

 Esre=real(Esr/max(max(Esr)));
 map(Esre,zi,xvero)
 axis([0,zi(end),0,10])
 pausak
 sog=5e-2;
 fi=find(abs(Esre)>sog);
 Ess=abs(Esre);
 Ess(fi)=sog;
 map(Ess,zi,xvero)
 axis([0,zi(end),0,10])
 pausak
