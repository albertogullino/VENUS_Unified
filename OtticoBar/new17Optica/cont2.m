for f10=fsovi
 kp=f10;

 ip=1;
 ys=[];
 fs=[];
 while kp<=min([length(Fi1du) length(ya1)])-1
% while kp<=length(ya1)-1
  de=ya1(kp+1)-ya1(kp)
  if isnan(de)==0
  if de>0
   ys(ip)=ya1(kp);
   ys(ip+1)=ya1(kp+1);
   fs(ip)=Fi1du(kp);
   fs(ip+1)=Fi1du(kp+1);
  else
   break
  end
  end
  ip=ip+1;
  kp=kp+1;
 end
%keyboard
 kp=f10;
 fd=[];
 yd=[];
 ip=1;
 'KP', keyboard
% kp=min([length(ya1) length(Fi1du)]);
% fina=find(ya1>0 & isnan(ya1)==0);
% kp=fina(end);
 while kp>1
  de=ya1(kp)-ya1(kp-1)
  pausak
  if isnan(de)==0 & de>0
   yd(ip+1)=ya1(kp-1);
   yd(ip)=ya1(kp);
   fd(ip+1)=Fi1du(kp-1);
   fd(ip)=Fi1du(kp);
  else
   if isnan(de)==0 
    break
   end 
  end
  kp=kp-1;
  ip=ip+1;
 end
 if length(fd)>0 & length(fs)>0
  Fi=[fliplr(fd) fs(2:length(fs))];
  ya=[fliplr(yd) ys(2:length(fs))];
 end
 if length(fd)>0 & length(fs)==0
  Fi=[fliplr(fd)];
  ya=[fliplr(yd)];
 end
 if length(fd)==0 & length(fs)>0
  Fi=fs;
  ya=ys;
 end

 if exist('ya')==0
  ya=[];
 end

 if length(ya)>=1
  dya=[1 diff(ya)./diff(Fi)];
  f2mv=(find(dya>0 & ya<0));
  f2pv=(find(dya>0 & ya>0));
  f2m1=length(f2mv);
  f2p1=length(f2pv);
  ie=0;
%  f1=find(dya>0 & abs(ya)<300);
  f1=find(dya>0 & abs(ya)<alMAX);
%  f1=find(dya>0 & abs(ya)<2e7);
  st=diff(Fint1(1:2));
  if f2m1*f2p1==0
   f2=find(dya>0 & abs(ya)<5);
%   f2=find(dya>0 & abs(ya)<100);
   ie=1;
  else
%   f2=find(dya>0 & abs(ya)<30);
   f2=find(dya>0);
  end
%  if icontr>=1
%   if exist('he'), close(he), clear he, end
%   he=figure; plot(Fi,ya,'*',Fi1,ya1), grid, pausak
%  end
  ya2=ya(f1);
  fa2=Fi(f1);
  dya3=[1 diff(ya2)./diff(fa2)];
  f2mv=(find(dya3>0 & ya2<0));
  f2pv=(find(dya3>0 & ya2>0));
  f2m=length(f2mv);
  f2p=length(f2pv);
  clear ya
  nosol=0;
 else
  nosol=1;
  f1=0;
  f2=0;
 end

if length(f1)>=2 & length(f2)>0 & nosol==0
%  [du,izp]=min(abs(ya2));
%  zes=fa2(izp);
  ly=length(ya2);
  fiy=find(ya2(1:ly-1).*ya2(2:ly)<0);
  puy=[fiy fiy+1];
  cozs=polyfit(fa2(puy),ya2(puy),1);
  zes=roots(cozs);
 if f2m>=2 & f2p>=2
  [du,iso]=sort(abs(ya2(f2mv)));
  yas1=ya2(f2mv(iso(1:2)));
  xas1=fa2(f2mv(iso(1:2)));
  [du,iso]=sort(abs(ya2(f2pv)));
  yas2=ya2(f2pv(iso(1:2)));
  xas2=fa2(f2pv(iso(1:2)));
  yf=[yas1 yas2];
  xf=[xas1 xas2];
  [du,isce]=sort(xf);
  xfu=xf(isce(1:3));
  yfu=yf(isce(1:3));
  coz=polyfit(xfu,yfu,length(xfu)-1);
  zep=roots(coz);
  ifz=find(imag(zep)==0);
  zep1=zep(ifz);
  if length(zep1)>1
   [du,izp1]=min(abs(zep1-zes));
   ze=zep1(izp1);
   if ze>fa2(puy(1)) & ze<fa2(puy(2))
   else
    [du,izp1]=max(abs(zep1-zes));
    ze=zep1(izp1);
   end
  else
   ze=zep1;
  end
 elseif f2m>=2 & f2p<2
  [du,iso]=sort(abs(ya2(f2mv)));
  yas1=ya2(f2mv(iso(1:2)));
  xas1=fa2(f2mv(iso(1:2)));
  yf=[yas1 ya2(f2pv)];
  xf=[xas1 fa2(f2pv)];
  coz=polyfit(xf,yf,length(xf)-1);
  zep=roots(coz);
  ifz=find(imag(zep)==0);
  zep1=zep(ifz);
  if length(zep1)>1
   [du,izp1]=min(abs(zep1-zes));
   ze=zep1(izp1);
   if ze>fa2(puy(1)) & ze<fa2(puy(2))
   else
    [du,izp1]=max(abs(zep1-zes));
    ze=zep1(izp1);
   end
  else
   ze=zep1;
  end
 elseif f2p>=2 & f2m<2
  [du,iso]=sort(abs(ya2(f2pv)));
  yas1=ya2(f2pv(iso(1:2)));
  xas1=fa2(f2pv(iso(1:2)));
  yf=[yas1 ya2(f2mv)];
  xf=[xas1 fa2(f2mv)];
  coz=polyfit(xf,yf,length(xf)-1);
  zep=roots(coz);
  ifz=find(imag(zep)==0);
  zep1=zep(ifz);
  if length(zep1)>1
   [du,izp1]=min(abs(zep1-zes));
   ze=zep1(izp1);
   if ze>fa2(puy(1)) & ze<fa2(puy(2))
   else
    [du,izp1]=max(abs(zep1-zes));
    ze=zep1(izp1);
   end
  else
   ze=zep1;
  end
 else
  [du,iso]=sort(abs(ya2));
  yas=ya2(iso);
  xas=fa2(iso);
  yf=yas(1:2);
  xf=xas(1:2);
  coz=polyfit(xf,yf,1);
  ze=roots(coz);
 end
  if icontr>=1
%   he=figure; subplot(211), plot(xf,yf,'*',Fi1,ya1,ze,0,'wo'), grid,
%   ax=axis; ax(3:4)=[-5 2]; axis(ax);
%   pausak
   if iins==4
    ichg=input(' cambio punti ? [0/1] ');
    if isempty(ichg)==1
     ichg=0;
    end
    if ichg==1
     pu=input(' punti = ');
     yf=Fi1(pu);
     xf=ya1(pu);
     coz=polyfit(xf,yf,1);
     ze=roots(coz);
     plot(xf,yf,'*',Fi1,ya1,ze,0,'wo'),
    end
   end
%   close(he), clear he
end
end
end
