%pograp0=[138     3   269   801];
peso=-[0:.1:1]*(1-j)/sqrt(2);
%peso=[0:.1:1];

aax=8;
for k=1:length(peso)
pes=peso(k)

Ex1=Exm{1}+pes*Exm{2};
Ey2=Eym{1}+pes*Eym{2};

Ex1=Exm{1};
Ey2=pes*Eym{2};

Ex1=Exm{2};
Ey2=pes*Eym{1};


%  s1=abs(Ex1.^2)-abs(Ey2.^2);
%  s2=real(Ex1.*conj(Ey2));
%  s3=imag(Ex1.*conj(Ey2));
s1=0;
s2=0;
s3=0;

%Ex1=pes*Exm{2};
%Ey2=pes*Eym{2};

  s1=s1+abs(Ex1.^2)-abs(Ey2.^2);
  s2=s2+real(Ex1.*conj(Ey2));
  s3=s3+imag(Ex1.*conj(Ey2));

  figure
  set(gcf,'Position',pograp);

  subplot(3,1,1)
        map_fnew(XP,YP,s1,aax,Cug.x,Cug.y,Cug.z,'S1',ibar)
  subplot(3,1,2)
        map_fnew(XP,YP,s2,aax,Cug.x,Cug.y,Cug.z,'S2',ibar)
  subplot(3,1,3)
        map_fnew(XP,YP,s3,aax,Cug.x,Cug.y,Cug.z,'S3',ibar)
 pausak  
 
   figure
   set(gcf,'Position',pograp);
 
   subplot(4,1,1)
         map_fnew(XP,YP,real(Exm{1}),aax,Cug.x,Cug.y,Cug.z,'real Ex1',ibar)
   subplot(4,1,2)
         map_fnew(XP,YP,real(Eym{1}),aax,Cug.x,Cug.y,Cug.z,'real Ey1',ibar)
   subplot(4,1,3)
         map_fnew(XP,YP,imag(Exm{1}),aax,Cug.x,Cug.y,Cug.z,'imag Ex1',ibar)
   subplot(4,1,4)
         map_fnew(XP,YP,imag(Eym{1}),aax,Cug.x,Cug.y,Cug.z,'imag Ey1',ibar)   

   figure
   set(gcf,'Position',pograp);
 
   subplot(4,1,1)
         map_fnew(XP,YP,real(Exm{2}),aax,Cug.x,Cug.y,Cug.z,'real Ex2',ibar)
   subplot(4,1,2)
         map_fnew(XP,YP,real(Eym{2}),aax,Cug.x,Cug.y,Cug.z,'real Ey2',ibar)
   subplot(4,1,3)
         map_fnew(XP,YP,imag(Exm{2}),aax,Cug.x,Cug.y,Cug.z,'imag Ex2',ibar)
   subplot(4,1,4)
         map_fnew(XP,YP,imag(Eym{2}),aax,Cug.x,Cug.y,Cug.z,'imag Ey2',ibar)   
 
end        
