function T=tr_i(k,lambda,rr,itetm,L,r)
inc=0;
Lg=L;
ifp=0;
verca=0;
%' conat', keyboard
if length(r)>1

 ND=100;
 n0=r(1);
 ns=r(2);
 ln=log(ns/n0);
 xp=linspace(0,L,ND);
 np0=n0*exp(xp/L*ln);
 nl=linspace(n0,ns,ND);
% figure, plot(xp,np0,'g.-',xp,nl), pausak
 inc=1;
end
be=(sqrt(1-k.^2));

if inc==0
 be0=-j*2*pi/lambda*rr*L;
 if itetm==1
  Z=1./be;
  Kz=0;
 elseif itetm==2
  Z=be;
  depez=1-(rr/r).^2;
  Kz=[1 -1; 1 -1]*Z*depez/2*(k^2./(1-k^2));
 end
 Idel=Z/2;
  if inc==0
   depep=(r/rr).^2-1;
  else
   depepi=((rv/rr).^2-1);
   depep=sum(depepi)*dx;
%   figure, plot(xx,depepi);
%   title('depep')
%   pausak
  end
 Kt=[1 1; -1 -1]*Idel*depep;
 P=be0*[be 0; 0 -be]+ be0*(Kt + Kz);
 T=expm(P);
else

 ln=log(ns/n0);
 gam=ln/L;
 k0=2*pi/lambda;
 x0=k0*n0/gam;
 xs=k0*ns/gam;
% n1=n0;
% n2=ns;
 n1=rr;
 n2=rr;
 kr0=k0*rr;
 k00=k0*n0;
 k0s=k0*ns;
 k01=k0*n1;
 k02=k0*n2;
 kt=k*kr0;
 be1=sqrt(k01^2-kt^2);
 be2=sqrt(k02^2-kt^2);
 bk0=be1/k00;
 bks=be2/k0s;
 bk1=be1/k01;
 bk2=be2/k02;
 z1=be1/(k0*n1^2);
 z2=be2/(k0*n2^2);


 if itetm==1
  nu0=sqrt((kt/gam)^2);
  nu=nu0+[-1 0 1];
  j0=(besselj(nu,x0));
  js=(besselj(nu,xs));
  y0=(bessely(nu,x0));
  ys=(bessely(nu,xs));

  j0p=(j0(1)-j0(3))/2;
  y0p=(y0(1)-y0(3))/2;
  jsp=(js(1)-js(3))/2;
  ysp=(ys(1)-ys(3))/2;


  J0=bk0*j0(2)-j*j0p;
  Js=js(2)+j*jsp/bks;
  Y0=bk0*y0(2)-j*y0p;
  Ys=ys(2)+j*ysp/bks;
  J0a=bk0*j0(2)+j*j0p;
  Jsa=js(2)-j*jsp/bks;
  Y0a=bk0*y0(2)+j*y0p;
  Ysa=ys(2)-j*ysp/bks;

  Den=(j0(2)*y0p-j0p*y0(2));
  T(1,1)=+J0*Ys-Y0*Js;
  T(1,2)=-J0a*Ys+Y0a*Js;
  T(2,1)=+J0*Ysa-Y0*Jsa;
  T(2,2)=+Y0a*Jsa-J0a*Ysa;
  Tdu=-T*j/(2*Den);
  if ifp==1
   ' verigica TE<', keyboard
  end
  if verca==1
   T=Tdu;
   z=linspace(0,L,50);
   zn=linspace(-L,0,50);
   C=1;
   D=-T(2,1)/T(2,2)*C;
   F=(T(1,1)-T(1,2)*T(2,1)/T(2,2))*C;
   A=-j/Den*(-C*Y0+D*Y0a);
   B=-j/Den*(C*J0-D*J0a);
   Fv=A*js(2)+B*ys(2);
   f1=C*exp(-j*be1*zn)+D*exp(j*be1*zn);
   x=k0*n0*exp(z*gam)/gam;
   Nu0=nu(2);
   fg=A*(besselj(Nu0,x))+B*(bessely(Nu0,x));
   f2=F*exp(-j*be2*z);
   i1=-j*be1*(C*exp(-j*be1*zn)-D*exp(j*be1*zn));
   x=k0*n0*exp(z*gam)/gam;
   Nu0=nu([1 3]);
%   ' ver', keyboard
   jb=(besselj(Nu0',x)).';
   jp=0.5*(jb(1,:)-jb(2,:));
   yb=(bessely(Nu0',x)).';
   yp=0.5*(yb(1,:)-yb(2,:));
   ig=gam*x.*(A*jp+B*yp);
   i2=-j*be2*F*exp(-j*be2*z);
   zt=[zn z z+L];
   vt=[f1 fg f2];
   it=[i1 ig i2];
   ntot=[ones(size(z))*n1 x*gam/k0 ones(size(z))*n2];
   figure, plot(zt,abs(vt),zt,real(vt),zt,imag(vt)),
   hold on, plot(zt,ntot,'w.'),
   title(' tensione TE')
   pausak
   figure, plot(zt,abs(it),zt,real(it),zt,imag(it)),
   hold on, plot(zt,ntot,'w.'),
   title(' corrente TE ')
   pausak
   keyboard
  end

 else

  nu0=sqrt((kt/gam)^2+1);
  nu=nu0+[-1 0 1];
  j0=(besselj(nu,x0));
  js=(besselj(nu,xs));
  y0=(bessely(nu,x0));
  ys=(bessely(nu,xs));

  j0p=(j0(1)-j0(3))/2;
  y0p=(y0(1)-y0(3))/2;
  jsp=(js(1)-js(3))/2;
  ysp=(ys(1)-ys(3))/2;
  K0=j0p+j0(2)/x0;
  Ks=jsp+js(2)/xs;
  H0=y0p+y0(2)/x0;
  Hs=ysp+ys(2)/xs;
  J0=j0(2);
  Js=js(2);
  Y0=y0(2);
  Ys=ys(2);
%  Csi=bk2*n2/n0*Ks-j*ns/n0*Js;
%  Zet=bk2*n2/n0*Hs-j*ns/n0*Ys;
%  Tet=H0+j*bk1*n0/n1*Y0;
%  Tau=K0+j*bk1*n0/n1*J0;

  Csi=Ks/z2-j*ns*Js;
  Zet=Hs/z2-j*ns*Ys;
  Tet=z1*n0*Y0-j*H0;
  Tau=z1*n0*J0-j*K0;

  Csia=Ks/z2+j*ns*Js;
  Zeta=Hs/z2+j*ns*Ys;
  Teta=z1*n0*Y0+j*H0;
  Taua=z1*n0*J0+j*K0;

  Den=(J0*H0-Y0*K0);
  T(1,1)=-Csi*Tet+Zet*Tau;
  T(1,2)=-Csi*Teta+Zet*Taua;
  T(2,1)=-Csia*Tet+Zeta*Tau;
  T(2,2)=(-Csia*Teta+Zeta*Taua);
  Tdu=T/(2*Den*n0);

%  T(1,1)=Csi*Tet-Zet*Tau;
%  T(1,2)=-(Csi*Tet'-Zet*Tau');
%  T(2,1)=Csi'*Tet-Zet'*Tau;
%  T(2,2)=-(Csi*Tet-Zet*Tau)';
%  Tdu=T*j/(2*Den);

  if ifp==1
    ' verigica TM<', keyboard
  end
  if verca==1
   T=Tdu;
   z=linspace(0,L,50);
   zn=linspace(-L,0,50);
   C=1;
   D=-T(2,1)/T(2,2)*C;
   F=(T(1,1)-T(1,2)*T(2,1)/T(2,2))*C;
%   A=j/Den*(-C*Tau+D*Tau');
%   B=j/Den*(C*Tet-D*Tet');
%   Fv=-j*(ns/n0)*(A*Js+B*Ys);
   B=1/(Den*n0)*(C*Tau+D*Taua);
   A=-1/(Den*n0)*(C*Tet+D*Teta);
%   Fv=(A*Ks+B*Hs)/z2;
%   Fv1=-j*ns*(A*Js+B*Ys);
%   rappf=Fv/F
%'   controllo rapporti'
%   pausak
   f1=j*(C*exp(-j*be1*zn)-D*exp(j*be1*zn));
   Nu0=nu(2);
   x=k0*n0*exp(z*gam)/gam;
%   nbar=exp(z*gam);
   nbar=n0*exp(z*gam);
   fg=nbar.*(A*(besselj(Nu0,x))+B*(bessely(Nu0,x)));
   f2=j*F*exp(-j*be2*z);

   Nu0=nu;
   jb=(besselj(Nu0',x)).';
   jp=0.5*(jb(1,:)-jb(3,:));
   yb=(bessely(Nu0',x)).';
   yp=0.5*(yb(1,:)-yb(3,:));
   K=jp+jb(2,:)./x;
   H=yp+yb(2,:)./x;

   v1=j*z1*(C*exp(-j*be1*zn)+D*exp(j*be1*zn));
   vg=j*(A*K+B*H);
   v2=j*z2*F*exp(-j*be2*z);
   zt=[zn z z+L];
   it=[f1 fg f2];
   vt=[v1 vg v2];
   ntot=[ones(size(z))*n1 x*gam/k0 ones(size(z))*n2];
   figure, plot(zt,abs(vt),zt,real(vt),zt,imag(vt)),
   hold on, plot(zt,ntot,'w.'),
   title(' tensione TM ')
   pausak
   figure, plot(zt,abs(it),zt,real(it),zt,imag(it)),
   hold on, plot(zt,ntot,'w.'),
   title(' corrente TM ')
   pausak
   keyboard
  end
 end


%  nout=n2;
%  ning=n1;
%
%  ru=rr/nout;
%  ri=rr/ning;
%  K=k;
%   bei=(sqrt(1-(ri*K.').^2));
%   beu=(sqrt(1-(ru*K.').^2));
%   ber=(sqrt(1-(K.').^2));
%  if itetm==1
%   Zi=1./(ning*bei);
%   Zu=1./(nout*beu);
%   Zr=1./(rr*ber);
%  else
%   Zi=bei./ning;
%   Zu=beu./nout;
%   Zr=ber./rr;
%  end

%  Tdu=Tdu*sqrt(Zi/Zu);

%  Gau=(Zu-Zr)./(Zu+Zr);
%  Gai=(Zi-Zr)./(Zi+Zr);
%  pz=(Zu.*Zr);
%  if imag(pz)==0
%   pz=abs(pz);
%  end
%  Tu=2*sqrt(pz)./(Zu+Zr);
%  fi=find(imag(Tu)~=0);
%  Tu(fi)=0;
%  pz=(Zi.*Zr);
%  if imag(pz)==0
%   pz=abs(pz);
%  end
%  Ti=2*sqrt(pz)./(Zi+Zr);
%  fi=find(imag(Ti)~=0);
%  Ti(fi)=0;
%
%  Tti=1/Ti*[1 -Gai; -Gai 1];
%  Ttu=1/Tu*[1 Gau; Gau 1];
%  T=Ttu*Tdu*Tti;

  T=Tdu;

%  ' controllo T ', keyboard
end
if ifp==1
'controllo T'
pausak
end

