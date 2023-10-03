return
Poim=sqrt(Pvectx.^2+Pvecty.^2+Pvectz.^2);
 Xpr=sin([teR])*cos(fian);
 Ypr=sin([teR])*sin(fian);
 Zpr=cos([teR])*ones(size(fian));
 Px_n=Pvectx.*Xpr;
 Py_n=Pvectx.*Ypr;
 Pz_n=Pvectz.*Zpr;
 Poim_n=sqrt(Px_n.^2+Py_n.^2+Pz_n.^2);

figure, plot(Poim(:,46)), hold on, plot(Poim_n(:,46),'r'), pausak


figure, plot(Poim(:,46)), hold on, plot(Pvectx(:,46),'r')
plot(Pvecty(:,46),'g'), plot(Pvectz(:,46),'b'), pausak
plot(Poim_n(:,46),'.'),  plot(Px_n(:,46),'r.')
plot(Py_n(:,46),'g.'), plot(Pz_n(:,46),'b.'), pausak
map(Poim,X1,Y1), pausak
map(Pvectx,X1,Y1), pausak
map(Pvecty,X1,Y1), pausak
map(Pvectz,X1,Y1), pausak
