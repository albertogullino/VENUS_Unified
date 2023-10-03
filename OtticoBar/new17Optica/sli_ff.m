axli=30;
 roxt=1e4*tan(X(:,1)*pi/180); z0c=1;
 Ef2=Ef.^2;
 fi=find(Ef2<=1e-6);
 Ef2(fi)=1e-6;
 [du,im]=max(max(Ef2));
 fim=find(Ef2==du);
 fim0=fim(1);
 x00=abs(X(fim0));
 y00=abs(Y(fim0));
 x0=-max(max(X))*.8; x1=-x0; y0=y00*.9; y1=y00*1.1;
 fi_x=find(( X>x0&X<x1) & (Y>y0&Y<y1) );
 y0=-max(max(X))*.8; y1=-y0; x0=x00*.9; x1=x00*1.1;
 fi_y=find(( X>x0&X<x1) & (Y>y0&Y<y1) );
 roM=max(roxt)*1e-4;
 Xi=linspace(-1,1,2000)*roM*.8;
 Yix=z0c*tan(x00*pi/180);
 Yiy=z0c*tan(y00*pi/180);
 Fix=atan(Xi/Yix);
 Fiy=atan(Xi/Yiy)+pi/2;
 fim=find(Fix<0);
 Fix(fim)=Fix(fim)+2*pi;
 fim=find(Fiy<0);
 Fiy(fim)=Fiy(fim)+2*pi;
 Rix=sqrt(Yix.^2+Xi.^2);
 Riy=sqrt(Yiy.^2+Xi.^2);
 rX=roxt*1e-4;
 figure, polar(Fix,Rix), hold on, polar(Fiy,Riy,'r'),

 Zx=interpn(rX,fian,log(Ef2),Rix,Fix,'spline');
 FFx=exp(Zx);
 FFx=FFx/max(FFx);
 Zy=interpn(rX,fian,log(Ef2),Riy,Fiy,'spline');
 FFy=exp(Zy);
 FFy=FFy/max(FFy);
 Xid=atan(Xi/z0c)*180/pi;
 Yid=(atan(Xi/z0c))*180/pi;

 figure,
 plot(Xid,FFx,'r',Yid,FFy,'b'), hold on,
 plot(X(fi_x),Ef2(fi_x),'.',Y(fi_y),Ef2(fi_y),'g.')
 figure, plot(Xid,FFx,'r',Yid,FFy,'b'),
 axis([[-1 1]*axli/2 0 1])
 pausak, close


 fila=find(abs(FFx-.5)<1e-2);
 du=find(diff([0 fila])>1);
 fila=fila(du);
 lalo=diff(Xid(fila(1:2)));
 dilo=diff(Xid(fila(2:3)))+lalo;
 tilx=[' lobe= ',num2str(lalo),'  apert= ',num2str(dilo)];
 figure,
 plot(Xid,FFx,'r',Xid(fila),FFx(fila),'y.'), %grid
 title(tilx)
 axis([[-1 1]*axli/2 0 1])
 pausak, close

 FFm=sqrt(FFx.*FFy);
 FFm=FFm./max(FFm);

 fila=find(abs(FFm-.5)<1e-2);
 du=find(diff([0 fila])>1);
 fila=fila(du);
 lalo=diff(Xid(fila(1:2)));
 dilo=diff(Xid(fila(2:3)))+lalo;
 lalom=lalo;
 tilm=[' lobe= ',num2str(lalo),'  apert= ',num2str(dilo)];
 figure,
 plot(Xid,FFm,'r',Xid(fila),FFm(fila),'y.'), %grid
 title(tilm)
 axis([[-1 1]*axli/2 0 1])
 pausak, close

 fila=find(abs(FFy-.5)<1e-2);
 du=find(diff([0 fila])>1);
 fila=fila(du);
 lalo=diff(Xid(fila(1:2)));
 dilo=diff(Xid(fila(2:3)))+lalo;
 tily=[' lobe= ',num2str(lalo),'  apert= ',num2str(dilo)];


 figure,
 subplot(211)
 plot(Xid,FFx,'r'), %grid
 tilxm=[tilx,'  apertM= ',num2str(lalom)];
 title(tilxm)
 axis([[-1 1]*axli/2 0 1])
 subplot(212)
 plot(Yid,FFy), %grid
 title(tily)
 axis([[-1 1]*axli/2 0 1])
 pausak



% n=[4 6 8 10];
% n1=[4 6 8 10];
% ap=[2.39 1.62 1.2 1.0];
% ap_ea=[2.58 1.69 1.24];
% div_ea=[5.44 5.5 5.504 5.504 ];
% apm=[2.1 1.44 1.07 0.89];
% du=ap.*n;
% co=polyfit(n,du,0);
% ap0=co./n1;
% figure, plot(n1,ap0,n,ap,'o',n1,apm,'d')

