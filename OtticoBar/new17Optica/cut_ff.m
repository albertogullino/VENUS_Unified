fi0=find(abs((fian0-pi))<1e-4);
xc1=flipud(X(:,fi0));
ff1=flipud(Ef(:,fi0));
ff1(end)=ff1(end)*1.12;
maff1=ff1(end);
fi0=find(fian0==0);
xc1=[xc1; X(2:end,fi0)];
ff1=[ff1; Ef(2:end,fi0)];
xc1f=linspace(-30,30,200);
ff1f=spline(xc1,ff1,xc1f);

fi0=find(abs((fian0-pi/2*3))<1e-4);
xc2=flipud(Y(:,fi0));
ff2=flipud(Ef(:,fi0));
ff2(end)=ff2(end)*1.12;
maff2=ff2(end);
fi0=find(abs((fian0-pi/2))<1e-4);
xc2=[xc2; Y(2:end,fi0)];
ff2=[ff2; Ef(2:end,fi0)];
xc2f=linspace(-30,30,200);
ff2f=spline(xc2,ff2,xc2f);
figure, plot(xc1,ff1/maff1,'.',xc2,ff2/maff2,'.'), hold on
plot(xc1f,ff1f/maff1,xc2f,ff2f/maff2)
