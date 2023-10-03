Ex1=E2xo;
Ey2=E2yo;

  s1=abs(abs(Ex1.^2)+abs(Ey2.^2));
  s2=-(conj(Ex1).*Ey2+conj(Ey2).*Ex1);
  s3=j*(conj(Ex1).*Ey2-conj(Ey2).*Ex1);

  figure
  set(gcf,'Position',pograp);
ibar=1;
titl='S1';
  subplot(3,1,1)
        map_fnew(XP,YP,s1,aax,Cug.x,Cug.y,Cug.z,titl,ibar)
titl='S2';        
  subplot(3,1,2)
        map_fnew(XP,YP,s2,aax,Cug.x,Cug.y,Cug.z,titl,ibar)
titl='S3';        
  subplot(3,1,3)
        map_fnew(XP,YP,s3,aax,Cug.x,Cug.y,Cug.z,titl,ibar)
