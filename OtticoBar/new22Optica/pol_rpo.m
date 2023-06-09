
 alma=4;
 poy=200;
 pox=100;
 Exdo=Pol.Ex;
 Eydo=Pol.Ey;

 MY=max(max(Eydo));
 MX=max(max(Exdo));
 if abs(MX)>abs(MY)
  N=MX;
  cht=' Ex ';
  Estro=Exdo/N;
  Esmal=Eydo/N;
 else
  N=MY;
  cht=' Ey ';
  Estro=Eydo/N;
  Esmal=Exdo/N;
 end

Phim=zeros(size(Estro))*NaN;
Xim=zeros(size(Estro))*NaN;

%figure
for kce=1:length(cce)
 xinf=real(cce(kce))-2;
 xsup=real(cce(kce))+2;
 yinf=imag(cce(kce))-2;
 ysup=imag(cce(kce))+2;
% [xinf xsup yinf ysup]
 fis=find(XP>=xinf & XP<=xsup & YP>=yinf & YP<=ysup);
%figure, plot3(XP(fis),YP(fis),ones(size(fis)),'r.'), view(2), axis equal
% fisx=find(XP(>=xinf & XP<=xsup );
% fisy=find(YP>=yinf & YP<=ysup);
 E=Estro(fis); 
 e=Esmal(fis); 
 Pt=abs(E).^2+abs(e).^2;
 P=abs(E).^2;
% plot3(XP(fis),YP(fis),Pt,'.'), view(2), axis equal
% pausak
 mE=mean(E);
 me=mean(e);
 Sto=e.*conj(E)./P;
 Phim(fis)=mean(real(Sto));
 Xim(fis)=mean(imag(Sto));
end

 figure
 pograp=[100   200   550   700];
 set(gcf,'Position',pograp)
 subplot(2,1,1)
 fang=180/pi;

 map_fnew(XP,YP,fang*Phim,aax,Cug.x,Cug.y,Cug.z,'Phi',ibar,iaoff)

 subplot(2,1,2)
 map_fnew(XP,YP,fang*Xim,aax,Cug.x,Cug.y,Cug.z,'Xi',ibar,iaoff)
 pausak

 MY=max(max(Eydo));
 MX=max(max(Exdo));
 if abs(MX)>abs(MY)
  N=MX;
  Estro=Exdo/N;
  Esmal=Eydo/N;
 else
  N=MY;
  Estro=Eydo/N;
  Esmal=Exdo/N;
 end

 figure
 Ptot1=abs(Estro).^2;
 Ptot2=abs(Esmal).^2;
 Ptot=Ptot1+Ptot2;
  titl=['Total power: dominant field ',cht];
 map_fnew(XP,YP,Ptot,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)


 figure;
 pograp=[pox   poy   680   680];
 set(gcf,'Position',pograp)
 subplot(2,2,1)
 titl='real Estro';
 map_fnew(XP,YP,real(Estro),aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
 titl='imag Estro';
 subplot(2,2,3)
 map_fnew(XP,YP,imag(Estro),aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)

 subplot(2,2,2)
 titl='real Esmal';
 map_fnew(XP,YP,real(Esmal),aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)

 titl='imag Esmal';
 subplot(2,2,4)
 map_fnew(XP,YP,imag(Esmal),aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)







alfav=linspace(0,alma,alma+1);
%alfav=4;
supo0=30;

for ia=1:length(alfav)
 alfa=alfav(ia);
 supo=supo0*ia;
 salfa=sin(alfa*pi/180);
 Pup=abs(Esmal).^2+abs(salfa*Estro).^2-2*salfa*real(Esmal.*conj(Estro));

 alfa=-alfav(ia);
 salfa=sin(alfa*pi/180);
 Pum=abs(Esmal).^2+abs(salfa*Estro).^2-2*salfa*real(Esmal.*conj(Estro));

 figure
 pograp=[pox+supo   poy   550   700];
 set(gcf,'Position',pograp)
 subplot(2,1,1)
 titlm=num2str(alfa);
 titlp=num2str(-alfa);
 map_fnew(XP,YP,Pum,aax,Cug.x,Cug.y,Cug.z,titlm,ibar,iaoff)

 subplot(2,1,2)
 map_fnew(XP,YP,Pup,aax,Cug.x,Cug.y,Cug.z,titlp,ibar,iaoff)
 pausak

end


return
 figure
 Ptot=abs(Estro).^2+abs(Estro).^2;alfav=linspace(0,2,11);
 titl='Total power';supo0=30;
 map_fnew(XP,YP,Ptot,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)

for ia=1:length(alfav)
 alfa=alfav(ia);
 supo=supo0*ia;
 salfa=sin(alfa*pi/180);
 calfa=cos(alfa*pi/180);
 Pup=abs(Esmal*salfa).^2+abs(calfa*Estro).^2+2*salfa*calfa*real(Esmal.*conj(Estro));

 alfa=-alfav(ia);
 salfa=sin(alfa*pi/180);
 calfa=cos(alfa*pi/180);
 Pum=abs(Esmal*salfa).^2+abs(calfa*Estro).^2+2*salfa*calfa*real(Esmal.*conj(Estro));

 figure
 pograp=[pox+supo   poy   550   700];
 set(gcf,'Position',pograp)
 subplot(2,1,1)
 titlm=num2str(alfa);
 titlp=num2str(-alfa);
 map_fnew(XP,YP,Pum,aax,Cug.x,Cug.y,Cug.z,titlm,ibar,iaoff)

 subplot(2,1,2)
 map_fnew(XP,YP,Pup,aax,Cug.x,Cug.y,Cug.z,titlp,ibar,iaoff)
 pausak

end


% s1=abs(Ex1.^2)-abs(Ey2.^2);
% s2=-(conj(Ex1).*Ey2+conj(Ey2).*Ex1);
% s3=j*(conj(Ex1).*Ey2-conj(Ey2).*Ex1);


firat=Esmal.*conj(Estro)/abs(N^2);

Xi=imag(firat);
Phi=real(firat);
%fiI=find(isinf(abs(Xi)));
%Xi(fiI)=0;
%fiI=find(isinf(abs(Phi)));
%Phi(fiI)=0;
%
%fiI=find(abs(Xi)>.05);
%Xi(fiI)=0;
%fiI=find(abs(Phi)>.05);
%Phi(fiI)=0;

 subplot(2,3,3)
 titl=' Phi';
 map_fnew(XP,YP,Phi,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)

 titl='Xi';
 subplot(2,3,6)
 map_fnew(XP,YP,Xi,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)




 titl=' Alpha';
 figure
 map_fnew(XP,YP,Phi*180/pi,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)

 figure
 pograp=[300   200   550   700];
 set(gcf,'Position',pograp)
 subplot(2,1,1)
 titl=' + ';
 map_fnew(XP,YP,Xi.^2+(Phi-max(max(Phi))).^2,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)

 subplot(2,1,2)
 titl=' - ';
 map_fnew(XP,YP,Xi.^2+(Phi-min(min(Phi))).^2,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
