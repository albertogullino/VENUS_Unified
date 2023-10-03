%save wo Anout Gas iFFsect lambda ze nv rr KKv besm Mvefm0 xvero aiat Pusas


 Ares=Anout.*(1+Gas);

% Ares=(Acoz(pu1,:)+Acoz(pu2,:));
 dz=.5;
 Zma=200;
 zi=0;
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
 [du,pu1]=min(abs(xvero-aiat/2));
 [du,pu2]=min(abs(xvero-aiat));

 fiz0=find(ztot==0);
% figure
 if exist('soloplot')
  if soloplot==1
   return
  end
 end

 figure, plot(ztot,Es([1 pu1 pu2],:)),
 xlabel(' zeta toward device (micron)')
 ylabel(' Es ')
 title(' in dis-fz: 3 radial coord: 0, Ox/2, Ox ')
 pausak



 fiz0=find(ztot==0);

 fiex=find(ztot<0);
 [du,fima]=max(max(Es(:,fiex)));
 [du,fima0]=max(Es(:,fiz0));
 Enor=1./Es(fima0,fiz0);
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

 map(Es,ztot,xvero)
 axis([ztot(1),ztot(end),0,10])
 pausak
 [du,pu1]=min(abs(xvero-aiat/2));
 [du,pu2]=min(abs(xvero-aiat));

 dx=diff(xvero);
 dx=[dx dx(end)];
 dx(end)=dx(end)/2;
 xdx=xvero.*dx;
 xdx=dx;
 figure,
 for kz=length(ztot):-1:1
  Ep=Es(:,kz)'.*xdx;
  Ed=Es(:,kz);
%  fis=find(Es(:,kz)>.1);
  ema=max(Ep);
  fis=find(Ep/ema>.1);
%%  plot(xvero,Es(:,kz),xvero(fis),Es(fis,kz),'.'), pausak
  xme=mean(xvero(fis));
  su=sum(dx(fis));
  xve(kz)=su;
  Eve(kz)=sum(Ep(fis))/(sum(dx(fis)));
%  plot(xvero,Ed,xvero(fis),Ed(fis),'.',su,Eve(kz),'ro'), pausak
 end
 figure, plot(ztot,Eve)
 figure, plot(ztot,xve)
 figure, plot(ztot,Eve.*xve)
