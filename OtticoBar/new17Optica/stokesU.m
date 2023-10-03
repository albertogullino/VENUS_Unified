Ex1=Ex;
Ey2=Ey;

  s1=abs(Ex1.^2)-abs(Ey2.^2);
%  s2=real(Ex1.*Ey2);
%  s3=imag(Ex1.*Ey2);
  s2=real(Ex1.*conj(Ey2));
  s3=imag(Ex1.*conj(Ey2));

  figure
  set(gcf,'Position',pograp);

  subplot(3,1,1)
        map_fnew(XP,YP,s1,aax,Cug.x,Cug.y,Cug.z,'S1',ibar)
  subplot(3,1,2)
        map_fnew(XP,YP,s2,aax,Cug.x,Cug.y,Cug.z,'S2',ibar)
  subplot(3,1,3)
        map_fnew(XP,YP,s3,aax,Cug.x,Cug.y,Cug.z,'S3',ibar)
