Ex1=Ex;
Ey2=Ey;

  s1=abs(abs(Ex1.^2)+abs(Ey2.^2));
  s2=-(conj(Ex1).*Ey2+conj(Ey2).*Ex1);
  s3=j*(conj(Ex1).*Ey2-conj(Ey2).*Ex1);

  figure
  pograp=[624    15   506   961];
  set(gcf,'Position',pograp);

XP1=xvero_tot;
YP1=ztot;
        map(s1,XP1,YP1), pausak
        map(s2,XP1,YP1), pausak
        map(s3,XP1,YP1)
