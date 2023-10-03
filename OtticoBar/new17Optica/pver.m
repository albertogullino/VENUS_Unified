r0ai=0
  cte(1)=(8*R*Delta)^2;
  cte(2)=0;
  cte(3)=-128*(R*Delta)^2;
  cte(4)=16*R*r0ai*Delta;
  cte(5)=64*(Delta*R)^2+16*(1+Delta)*Delta*R^2;
  cte(6)=-16*R*r0ai*Delta;
  cte(7)=-16*R^2*Delta*(1+Delta);
  cte(8)=2*r0ai*R*(1+Delta);
  cte(9)=(R^2*(1+Delta)^2+r0ai^2);
 rop=r0ai+R*(1+Delta*cos(4*fi)).*exp(j*fi);
 fx=fi;
 x=cos(fx);
 y=polyval(cte,x);
 Rte=R*(1+Delta*cos(4*fi));
 M=sqrt(r0ai^2+Rte.^2+2*r0ai*Rte.*cos(fi));
 figure, plot(fx,sqrt(y),fi,abs(rop),fi,M), pausak
 figure, plot(fi,M-abs(rop)), pausak
 figure, plot(fi,sqrt(y)-abs(rop)), pausak

fi0=.4;
 rop0=abs(r0ai+R*(1+Delta*cos(4*fi0)).*exp(j*fi0));
 x0=cos(fi0);
 y0=sqrt(polyval(cte,x0));
 Rte0=R*(1+Delta*cos(4*fi0));
 M0=sqrt(r0ai^2+Rte0.^2+2*r0ai*Rte0*cos(fi0));
 [rop0 y0 M0]
