pograp0=[138     3   269   801];
peso=-[0:.1:1]*(1-j)/sqrt(2);

aax=8;
for k=1:length(peso)
pes=peso(k)

Ex1=Exm{1};
Ey2=Eym{1};



  s1=abs(Ex1.^2)-abs(Ey2.^2);
  s2=real(Ex1.*conj(Ey2));
  s3=imag(Ex1.*conj(Ey2));

Ex1=pes*Exm{2};
Ey2=pes*Eym{2};

  s1=s1+abs(Ex1.^2)-abs(Ey2.^2);
  s2=s2+real(Ex1.*conj(Ey2));
  s3=s3+imag(Ex1.*conj(Ey2));

  figure
  set(gcf,'Position',pograp0);

  subplot(5,1,1)
        map_fnew(XP,YP,s1,aax,Cug.x,Cug.y,Cug.z,'S1',ibar)
  subplot(5,1,2)
        map_fnew(XP,YP,s2,aax,Cug.x,Cug.y,Cug.z,'S2',ibar)
  subplot(5,1,3)
        map_fnew(XP,YP,s3,aax,Cug.x,Cug.y,Cug.z,'S3',ibar)
  subplot(5,1,4)
        map_fnew(XP,YP,abs(Ex1),aax,Cug.x,Cug.y,Cug.z,'Ex',ibar)      
  subplot(5,1,5)
        map_fnew(XP,YP,abs(Ey2),aax,Cug.x,Cug.y,Cug.z,'Ey',ibar)              
 pausak        
end        
