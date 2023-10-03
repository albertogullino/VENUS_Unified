ipa=1
if ipa==1
nsc=0;
Fint1=Fint(1:pak1:length(Fint));
%for nso=1:nma
%' nmasce ', keyboard
save sacont

for nso=nmascel
clear Fi1du
%disp(' acktung !!!!!!!!!!!! dissus '), keyboard
% for nso=2:nma

   ps=1:pak1:sMem(2);
   ya1=aou(nso,ps)-alpha_th;

%&&&&&&&&&&&&&
   fi=find(isnan(ya1)==0);
%   ya1=aou(nso,fi);
   ya1=ya1(fi);
%   ta=tetmu(nso,fi);
%   ga=gou(nso,fi);
%   Fi1=Fint1(fi);
   ta=tetmu(nso,ps);
   ta=ta(fi);
   ga=gou(nso,ps);
   ga=ga(fi);
   Fi1=Fint1(fi);
%   'ga', keyboard
   fisa=fi;
  if icontr==2
   h1=figure;
   subplot(121), plot(Fi1,ya1), grid,
   subplot(122), semilogy(Fi1,ga), grid,
  end
 dya1=[1 diff(ya1)./diff(Fi1)];
 dya1p=[diff(ya1)./diff(Fi1) 1];
 ly=[1:length(dya1)-1];
 dya2=[1 ya1(ly).*ya1(ly+1)];
% f10v=find(dya1>0 & dya1p>0 & dya2<0)
 f10v=find(dya1>0 & dya2<0);
% f10v=find(dya2<0);
  if length(f10v)>0
     fpu=f10v(1);
     fpA=fpu+[-2:2];
     fiA=find(fpA>=1 & fpA<=length(Fint));
     fpA=fpA(fiA);
     Adispla=abs(Amsort(:,fpA,nso));   
%     'A1', keyboard
  end

  if icontr==2
   figure(h1),
   subplot(121), hold on, plot(Fi1(f10v),ya1(f10v),'wo'), ,
   subplot(122), hold on, semilogy(Fi1(f10v),ga(f10v),'wo'), ,
   title(' controllo 5')
   disp(' controllo ')
   pausak
   close(h1)
  end


  ifnz=0;
  fsovi=[];
  if length(f10v)>1

   lumax=15;
%  lumax=270;
   if length(ya1)>lumax

    Fi1s=Fi1;
    ya1s=ya1;
    gasa=ga;

    zev=Fi1(f10v);

    gra=ceil(length(ya1)-1);
    pfo=linspace(Fi1(1),Fi1(end),100);

    co=polyfit(Fi1,log10(ga),gra);
    gf=10.^(polyval(co,pfo));
    g0v=10.^(polyval(co,zev));
    [du,igmi]=min(g0v);


    if icontr==2
     h1=figure;
     subplot(121), plot(Fi1s,ya1s,'.-',Fi1s(f10v),ya1s(f10v),'wo',zev,0,'r+'),
     subplot(122), semilogy(Fi1s,gasa,'.-',Fi1s(f10v),gasa(f10v),'wo'),
     disp(' controllo 10')
     pausak
    end
    ksi=0;
    if icontr==2
    hdu=figure;
    end
    clear zeu g0u
    for ks=1:length(zev)
     zel=zev(ks);
     [du,fi]=sort(abs(Fi1-zel));
     pu=[1:7];
     grl=length(pu)-1;
     fi=sort(fi(pu));
     fz=Fi1(fi);
     az=ya1(fi);

     co=polyfit(fz,az,grl);
     zed=roots(co);
     fia=find(abs(imag(zed))<=abs(real(zed)) & ...
         (real(zed)<fz(end) & real(zed)>fz(1)));
     if length(fia)>0
      ksi=ksi+1;
      zedu=zed(fia);
      [du,izeus]=min(abs(zedu-zel));
      zeus=zedu(izeus);
      zeu(ksi)=zeus;
      co=polyfit(fz,log10(ga(fi)),grl);
      g0k=10.^(polyval(co,zeus));
      g0u(ksi)=g0k;
       if icontr==2
        figure(hdu)
        subplot(121), plot(Fi1s,ya1s,'.-',Fi1s(f10v),ya1s(f10v),'wo',fz,az,'c.',zeus,0,'r*'),
        subplot(122), semilogy(Fi1s,gasa,'.-',Fi1s(f10v),gasa(f10v),'wo',fz,ga(fi),'c.',zeus,g0k,'r*'),
        disp(' controllo ')
        pausak
       end

     end
    end
     
    [g0s,ifs]=min(g0u);
    zes=zeu(ifs);

    ze=zes;
    gg0=g0s;
    [du,ivero]=min(abs(Fi1-ze));
    fsovi=ivero;
     
    if icontr==2
     close(hdu)
     h1=figure;
     subplot(121), plot(Fi1s,ya1s,'.-',Fi1s(f10v),ya1s(f10v),'wo',zeu,0,'r+',zes,0,'go'),
     subplot(122), semilogy(Fi1s,gasa,'.-',Fi1s(f10v),gasa(f10v),'wo',zeu,g0u,'r+',zes,g0s,'go'),
     disp(' controllo ULTIMO')
     pausak

     [du,igmi]=min(gasa(f10v));
     fpu=f10v(igmi);
     fpA=fpu+[-2:2];
     fiA=find(fpA>=1 & fpA<=length(Fint));
     fpA=fpA(fiA);     
     Adispla=abs(Amsort(:,fpA,nso));
%          'A5', keyboard

%     figure, semilogy(Adispla), pausak

     close(h1)
    end

    fr_pun=[];
    for kze=1:length(zeu)
     zelo=zeu(kze);
     [du,ize]=min(abs(Fint-zelo));
     fr_pun=[fr_pun ize];
    end
%    'passo 1', keyboard
     clear As
     for iaf=1:length(Fint)
      Pun=pou(nso,iaf);
      if Pun>0
       Ada=reshape(Anu(:,Pun,iaf),sA(1),1);
      else
       Ada=zeros(sA(1),1);
      end
      [du,ima]=max(abs(Ada));
      As(:,iaf)=Ada*sign(real(Ada(ima)));
     end
    lep=size(pou);     
 for iaf=1:lep(2)
  Pun=pou(nso,iaf);
  if Pun>0
   Ada=reshape(Anu(:,Pun,iaf),sA(1),1);
  else
   Ada=zeros(sA(1),1);
  end
%  [du,ima]=max(abs(Ada));
%  As(:,iaf)=Ada*sign(real(Ada(ima)))/median(abs(Ada));
  As(:,iaf)=Ada/median(abs(Ada));
 end
 
 Apre=As(:,1);
 corp=Apre'*Apre;
for ksor=2:lep(2)
  Apo=As(:,ksor);
  corpi=Apre'*Apo;
  Apre=Apo;
  corp=corpi;
  if real(corpi)<0
   As(:,ksor)=-As(:,ksor);
   Apre=-Apo;
   corp=-corpi; 
  end
 end
 
%    field_fr

   else

    Fi1s=Fi1;
    ya1s=ya1;
    gasa=ga;

    if f10v(1)>2
     pumi=f10v(1)-2;
    else
     pumi=f10v(1)-1;
    end
    if f10v(end)<length(Fi1)
     puma=f10v(end)+1;
    else
     puma=f10v(end);
    end
    puf=pumi:puma;
    Fi1=Fi1(puf);
    ya1=ya1(puf);
    ga=ga(puf);
%    'ga2', keyboard


    gra=ceil(length(ya1)-1);
    co=polyfit(Fi1,ya1,gra);
    pfo=linspace(Fi1(1),Fi1(end),100);
    af=polyval(co,pfo);
    zed=roots(co);
    fia=find(imag(zed)==0 & (real(zed)<Fi1(end) & real(zed)>Fi1(1)));

    if icontr==2
     figure;
     plot(Fi1s,ya1s,Fi1,ya1,'w.',pfo,af,zed(fia),0,'ro'),
%     plot(Fi1s,ya1s,Fi1,ya1,'w.',pfo,af,Fi1s(f10v),0,'ro'),
    end

    df=diff(Fi1(1:2));
    if length(fia)==0
     fia=find(imag(zed)==0 & (real(zed)<Fi1(end)+df & real(zed)>Fi1(1)-df ));
    end
    zev=zed(fia);

%     figure, plot(Fi1,ya1,'g',pfo,af,'r',zev,zev*0,'w*'), pausak

    co=polyfit(Fi1,log10(ga),gra);
    gf=10.^(polyval(co,pfo));
    g0v=10.^(polyval(co,zev));
    [du,igmi]=min(g0v);
     fpu=zev(igmi);
    [du,fpu]=min(abs(fpu-Fint));
     fpA=fpu+[-2:2];
     fiA=find(fpA>=1 & fpA<=length(Fint));
     fpA=fpA(fiA);     
     Adispla=abs(Amsort(:,fpA,nso));
     if ifp==-10
     'A6', keyboard
     end
%     figure, plot(Fi1,ga,'c',ze,gg0,'wo',pfo,gf,'m',zev,g0v,'w*'), pausak
    if icontr==2
     h1=figure;
     subplot(121), plot(Fi1s,ya1s,Fi1s(f10v),ya1s(f10v),'wo',zev,0,'r+'),
     subplot(122), semilogy(Fi1s,gasa,Fi1s(f10v),gasa(f10v),'wo'),
     disp(' controllo questo ')
     pausak
    end
    ksi=0;
    for ks=1:length(zev)
     zel=zev(ks);
     [du,fi]=sort(abs(Fi1-zel));
     pu=[1:3];
     grl=length(pu)-1;
     fi=sort(fi(pu));
     fz=Fi1(fi);
     az=ya1(fi);
     co=polyfit(fz,az,grl);
     zed=roots(co);
%     fia=find(imag(zed)==0 & (real(zed)<fz(end) & real(zed)>fz(1)));
     fia=find(abs(imag(zed))<=abs(real(zed)) & ...
         (real(zed)<fz(end) & real(zed)>fz(1)));
     if length(fia)>0
      ksi=ksi+1;
      zedu=zed(fia);
      [du,izeus]=min(abs(zedu-zel));
      zeus=zedu(izeus);
      zeu(ksi)=zeus;
      co=polyfit(fz,log10(ga(fi)),grl);
      g0k=10.^(polyval(co,zeus));
      g0u(ksi)=g0k;
     end
    end
    if icontr==2
     figure(h1)
     subplot(122), hold on, plot(zeu,g0u,'r+'),
     disp(' controllo ')
     pausak
     close(h1)
    end

    [g0s,ifs]=min(g0u);
    zes=zeu(ifs);

    fr_pun=[];
    for kze=1:length(zeu)
     zelo=zeu(kze);
     [du,ize]=min(abs(Fint-zelo));
     fr_pun=[fr_pun ize];
    end
%    'passo 2', keyboard

%    field_fr


    ze=zes;
    gg0=g0s;
    [du,ivero]=min(abs(Fi1-ze));
    fsovi=ivero;
%   else
%    [du,imi]=min(ga(f10v));
%    fsovi=f10v(imi);
   end

  elseif length(f10v)==1
   fsovi=f10v;
  elseif length(f10v)==0
   fsovi0=find(abs(ya1)<alMAX);
%   fsovi0=find(abs(ya1)<2e7);
   npma=2;
   if length(fsovi0)>=npma
    fsovi=[];
     [du,idu]=min(abs(ya1));
     fsovi00=ya1(idu);
     if fsovi00>0
      ifnz=-1;
     else
      ifnz=1;
     end
   end
  end

sau=size(aou);
%'sau', keyboard
if sau(1)==1
 icas=length(find(aou>0))==length(Fint) | length(find(aou<0))==length(Fint);
 ya1=aou;
 Fi1du=Fint(1:length(ya1));
% 'Fi1 3', keyboard
 dya1=[1 diff(ya1)./diff(Fi1du)];
 dya1p=[diff(ya1)./diff(Fi1du) 1];
 ly=[1:length(dya1)-1];
 dya2=[1 ya1(ly).*ya1(ly+1)];
 f10v=find(dya1>0 & dya2<0);

 if icas==1 | length(f10v)==0
   ya=aou;
   yg=gou;
  if icas==1
   dya=diff(ya)./diff(Fint);
   fim=find(dya>0);
   if length(fim)==0
    ipu0=1;
    puas=[1 2];
   else
    if fim(end)+1==length(Fint)
      puas=[1:length(Fint)];
      if length(find(ya>0))==0
       ipu0=length(Fint);
      else
       ipu0=1;
      end
    else
     if fim(1)==1
      ipu0=1;
      puas=[fim fim(end)+1];
     else
      ipu0=length(Fint);
      puas=[fim(1)-1 fim];
     end
    end
   end
  else
   ipu0=1;
   puas=[1 2];
  end

   fiAz=ipu0;
   yas=ya(puas);
   ygs=yg(puas);
   xgn=Fint(puas);
   coz=polyfit(xgn,yas,1);
   ze=roots(coz);
   [mi,imi]=min(abs(ze-Fint));
   if mi>abs(diff(Fint(1:2))*2)
    ze=Fint(imi);
   end

   cog=polyfit(xgn,log10(ygs),length(xgn)-1);
   gg0=10^(polyval(cog,ze));
   ta=tetmu(nso,puas);
   cot=polyfit(xgn,ta,length(xgn)-1);
   tt0=(polyval(cot,ze));
   if gg0/vg<GMA
    nsc=nsc+1;
    ns=[ ns nso];
    ipu=[ipu ipu0];
    fso=[fso ze];
    gso=[gso gg0];
    'gso 1', keyboard
    aso=[aso 0];
    tso=[tso tt0];
   end
   sA=length(puA);
   An0=reshape(Anu(:,pou(nso,fiAz),fiAz),sA(1),1);
   Anso=[Anso An0];
   
' capeoi', keyboard   

   if icampi>=1

%    fieval
%   disp('camdu in diss_nst 1'), keyboard
%ifp=-10
    lep=size(pou);
    Fint=Fint(1:lep(2));
    iLP=iLP1;
    fie_new
    iLP=iLPr;

%    'ferma', keyboard
    if iLP==1
     rtetm=0.5;
     nuazi=0;
     polca=0;
     polratio=0;
     mrad=0;
    end

    rtetmv=[rtetmv rtetm];
    maziv=[maziv nuazi];
    mradv=[mradv mrad];
    polcav=[polcav polca];
    polrat=[polrat polratio];
        ' cont diss', keyboard
 %  disp('polcav'), pausak
   end





%   Dla_new=Dlam_mod;
%   dF=diff(Dlam_mod(1:2))/(Dlam_mod(3)-1);
%   dFv=[-dF dF];
%   Dla_new(1:2)=ze*1000*lambda+dFv;


 end

end  %sau

%'fsovi', keyboard
if ~exist('Fi1du')==1 
 Fi1du=Fi1;
end
for f10=fsovi
 kp=f10;

 ip=1;
 ys=[];
 fs=[];
 while kp<=min([length(Fi1du) length(ya1)])-1
% while kp<=length(ya1)-1
  de=ya1(kp+1)-ya1(kp);
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
% 'KP', keyboard
% kp=min([length(ya1) length(Fi1du)]);
% fina=find(ya1>0 & isnan(ya1)==0);
% kp=fina(end);
 while kp>1
  de=ya1(kp)-ya1(kp-1);
%  pausak
  if de>0
   yd(ip+1)=ya1(kp-1);
   yd(ip)=ya1(kp);
   fd(ip+1)=Fi1du(kp-1);
   fd(ip)=Fi1du(kp);
  else
   break
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

if (ie==1 & min(abs(xf-ze))<st) | ie==0
%  if icontr==1
%   he=figure; plot(Fi,ya,'*',Fi1,ya1), grid, pausak
%   close(he), clear he
%  end

 [du,fd]=sort(abs(Fi-ze));
 fd0=fd(1);
 [du,fd]=sort(abs(Fint1-ze));
 ffA=(fd(1:2));
 FiA=Fint1(ffA);

 [du,fd]=sort(abs(xf-ze));
 clear ff
 for kf=1:length(xf)
  ff(kf)=find(xf(kf)==Fi1);
 end
 ff=sort(ff);
 xg=Fi1(ff);
 yg=ga(ff);
 ygm=max(yg);
 yt=ta(ff);

%disp('key')
%keyboard
% if ygm>20*min(yg) & length(yg)>2
%  f=find(yg~=ygm);
%  yg=yg(f);
%  xg=xg(f);
%  yt=yt(f);
% end

% cog=polyfit(xg,yg,length(xg)-1);
% gg0=polyval(cog,ze);
  if length(xg)>2
   xgma=max(abs(xg));
   xgn=xg/xgma;
   cog=polyfit(xgn,log10(yg),length(xg)-1);
   gg0=10^(polyval(cog,ze/xgma));
  else
   cog=polyfit(xg,yg,length(xg)-1);
   gg0=polyval(cog,ze);
  end
 cot=polyfit(xg,yt,1);
 tt0=polyval(cot,ze);
  if icontr>=1
   falam=lambda*1000;
   he=figure; subplot(211), plot(xf*falam,yf,'*',Fi1du*falam,ya1,ze*falam,0,'wo'), grid,
   figure(he); subplot(212)
   semilogy(xg*falam,yg,ze*falam,gg0,'wo'), grid,
   title(' control 35')
   pausak
%   close(he), clear he
  end
% disp(' verifica '), pausak

% pausak


 if gg0/vg<GMA
  nsc=nsc+1;
  ns=[ ns nso*fd0./fd0];
%  ipu=[ipu fi(fd0)];
  ipu=[ipu fd0];
  fso=[fso ze];
  gso=[gso gg0];
      'gso 2', keyboard

%  disp('gso'), pausak
  aso=[aso 0];
  tso=[tso tt0];
%  fip=find(pou(nso,:)~=0);
%  for ifip=1:length(fip)
%   V(:,ifip)=(Anu(:,pou(nso,fip(ifip)),fip(ifip)));
%  end
%  sA=size(Amu);
  if icontr>=1
%   he=figure; subplot(211), plot(xf,yf,'*',Fi1,ya1,ze,0,'wo'), grid,
%   ax=axis; ax(3:4)=[-5 2]; axis(ax);
%   pausak
%   close(he), clear he
  end
  lF=length(Fint1);

  fiAz=find(FiA(1)==Fint);
  if exist('Kvet')
   KKdu=Kve(:,fiAz);
   if KKdu(1)<1e-10
    KKdu(1)=1e-10;
   end
   fiK=find(KKdu~=0);
   KK=KKdu(fiK);
   npk=length(KK);
  end
  pbA=[1:dnum*length(KK)];
  sA=length(puA);

%  An0=reshape(Anu(:,pou(nso,fiAz),fiAz),sA(1),1);
clear As
lep=size(pou);

% for iaf=1:length(Fint)
 for iaf=1:lep(2)
  Pun=pou(nso,iaf);
  if Pun>0
   Ada=reshape(Anu(:,Pun,iaf),sA(1),1);
  else
   Ada=zeros(sA(1),1);
  end
%  [du,ima]=max(abs(Ada));
%  As(:,iaf)=Ada*sign(real(Ada(ima)))/median(abs(Ada));
  As(:,iaf)=Ada/median(abs(Ada));
 end
 
 Apre=As(:,1);
 corp=Apre'*Apre;
for ksor=2:lep(2)
  Apo=As(:,ksor);
  corpi=Apre'*Apo;
  Apre=Apo;
  corp=corpi;
  if real(corpi)<0
   As(:,ksor)=-As(:,ksor);
   Apre=-Apo;
   corp=-corpi; 
  end
 end
 
  duded=((Fint-ze));
  fimaz=find(duded>0);
  isoz=fimaz(1)+[-1 0];
  DeFi=(diff(Fint(isoz)));
  dude=abs(Fint(isoz)-ze);
%  An0=(dude(1)/DeFi)*As(:,isoz(1))+(1-dude(2)/DeFi)*As(:,isoz(2));
  An0=(dude(2)/DeFi)*As(:,isoz(1))+(dude(1)/DeFi)*As(:,isoz(2));
 if ifp==-1
  KKd=1:length(An0);
  figure, semilogy(KKd,abs(As),KKd,abs(An0),'w.-'), pausak
 end
  if icam_fr==1
   KKd=1:length(An0);
   haut=figure;
%   plot(KKd,abs(As),KKd,abs(An0),'w.-'), pausak
%   ' An0 disply '
   fr_pun =1:length(Fint);
   field_fr
   keyboard
  end

  if icontr>=1
   figure, semilogy(abs(An0)/median(abs(An0)),'ro'), hold on,
   semilogy(Adispla),
   semilogy(abs(As(:,isoz(1:2))),'.')
%   figure, semilogy(abs(An0),'ro'), hold on,
%   semilogy(Adispla),   
   expla=' Coeff to be interpolated,  Red dots: solution ';
   title(expla)
   
   
   pausak
  end
%  'prima Anso'
 % keyboard
fiAv=length(find(isnan(An0))); 

if fiAv==0
  Anso=[Anso An0];
else 
 gso=gso(1:end-1);
 fso=fso(1:end-1);
end  
%  'dopo Anso'
%  keyboard

  if ifp~=-4
  disp('campi in diss_new')
  end
%' capeoi questa', keyboard   

  if icampi>=1 & fiAv==0
%%%%%%%%% routine plot campi
%disp(' prima di fieval in dissus '), keyboard

%   fieval
%   disp('camdu in diss_nst 2'), keyboard
   fa1d=fou(nso,ps);
   fa1d=fa1d(fi);
%   figure, plot(fa1d,ya1), pausak, plot(fa1d,ga), pausak
%   disp('fie_new in diss_nst 1'), keyboard
%ifp=-10
    Fint=Fint(1:lep(2));
    iLP=iLP1;
    fie_new
    iLP=iLPr;   

    'ferma fie?new', keyboard
   if iLP==1
    rtetm=0.5;
    nuazi=0;
    polca=0;
    polratio=0;
    mrad=0;
   end
   M2v=[M2v M2];
   rtetmv=[rtetmv rtetm];
   maziv=[maziv nuazi];
   mradv=[mradv mrad];
   polcav=[polcav polca];
   polrat=[polrat polratio];
%  disp('polcav'), pausak
  end


 end %if gg0
end %ie
end %length(f1)
end %for f10

if iins>1 & nsc==0
 fst_d=5;
else
 fst_d=1;
end

if abs(ifnz)==1 & length(fsovi)==0 & iins>=1

  ifi=length(find(ya1<0));
  if length(ya1)==ifi
   iso=[ifi-1 ifi];
   isop=iso;
  elseif ifi==0
   iso=[1 2];
   isop=iso;
  else
   [du,fid]=sort(abs(ya1));
   isop=sort(fid(1:2));
%   isop=fid;
   if isop(1)<length(isop)
    iso=[isop(1) isop(1)+1];
   else
    iso=[isop(1)-1 isop(1)];
   end
  end
  isoA=isop;
  yt=ta(iso);
  yg=ga(iso);
  yf=ya1(iso);
  xf=Fi1(iso);
  xg=xf;
  coz=polyfit(xf,yf,npma-1);
  zevt=sort(roots(coz));
%  pausak
  iz=find(imag(zevt)==0);
  if length(iz)==0
   isacc=0
  elseif length(iz)==1
   ze=zevt;
   if min(abs(ze-xf))<2*diff(sort(xf))*fst_d
    isacc=1;
   else
    isacc=0;
   end
  elseif length(iz)==2
   zev=sort(zevt(iz));
   xfs=(sort(xf));
   stv=diff(sort(xf));
   st2=stv(1);
   st1=stv(1);
   if zev(2)<Fint(1)
    zev=zev(2);
    st1=stv(1)*5;
   elseif zev(2)>Fint(length(Fint))
    zev=zev(1);
    if iins>1
     st1=stv(1)*fst_d;
    end
   end
   izea=find(zev>xfs(1)-st1 & zev<xfs(length(xfs))+st2);
   zevs=zev(izea);
   isacc=1;
   if length(zevs)>1
    if ifnz==-1
     ze=zevs(1);
    else
     ze=zevs(2);
    end
   elseif length(zevs)==1
     ze=zevs;
   elseif length(zevs)==0
    isacc=0;
   end
  end
 if nsc==0
  isacc=1;
 end
 if iins>2
  if icontr>=1
   he=figure; plot(xf,yf,'*',Fi1,ya1,ze,0,'wo'), grid,
%   ax=axis; ax(3:4)=[-5 2]; axis(ax);
   pausak
   isacc=input(' accetti soluzione ? [0/1]');
   if isempty(isacc)==1
    isacc=1;
   end
   close(he), clear he
  end
 end

 if isacc==1

  if iins<=2
   if icontr>=1
    he=figure; plot(xf,yf,'*',Fi1,ya1,ze,0,'wo'), grid,
    ax=axis; ax(3:4)=[-5 2]; axis(ax);
    pausak
    close(he), clear he
   end
  end
  ygm=max(yg);
  if ygm>10*min(yg) & length(yg>2)
   f=find(yg~=ygm);
   yg=yg(f);
   xg=xg(f);
   yt=yt(f);
  end
  if length(xg)>2
   cog=polyfit(xg,log10(yg),length(xg)-1);
   gg0=10^(polyval(cog,ze));
  else
   cog=polyfit(xg,yg,length(xg)-1);
   gg0=polyval(cog,ze);
  end
  cot=polyfit(xg,yt,1);
  tt0=polyval(cot,ze);
  if icontr>=1
   he=figure; semilogy(xg,yg,ze,gg0,'wo'), pausak
   close(he), clear he
  end
   if gg0/vg<GMA
    nsc=nsc+1;
    ns=[ ns nso];
    [du,fiz]=min(abs(Fint-ze));
    ipu=[ipu fiz];
    fso=[fso ze];
    gso=[gso gg0];
        'gso 3', keyboard

 %  disp('gso'), pausak
    aso=[aso 0];
    tso=[tso tt0];
 %   sA=size(Amu);
    sA=length(puA);

 %   An1=reshape(Amu(:,nso,fisa(isoA)),sA(1),length(isoA));
 %   An0=An0d(:,1);
 %   if icontr==1
 %    figure, plot(abs(An0d)), hold on, plot(abs(An0),'w*'), pausak
 %   end
   fib=find(pou(nso,:)~=0);
   Fz=Fint./Fint*1e10;
   Fz(fib)=Fint(fib);
   [du,fiAz]=min(abs(ze-Fz));

'alsp metodo', keyboard
    An0=reshape(Anu(:,pou(nso,fiAz),fiAz),sA(1),1);
    Anso=[Anso An0];


   disp('campi in diss_new altro')

 %  imod=1;
 %  disp('camdu in dissu3')
 %  keyboard
   if icampi>=1

%    fieval
%   disp('camdu in diss_nst 3'), keyboard
%ifp=-10
    Fint=Fint(1:lep(2));
    iLP=iLP1;
    fie_new
    iLP=iLPr;    
    'ferma', keyboard

    if iLP==1
     rtetm=0.5;
     nuazi=0;
     polca=0;
     polratio=0;
     mrad=0;
    end

    rtetmv=[rtetmv rtetm];
    maziv=[maziv nuazi];
    mradv=[mradv mrad];
    polcav=[polcav polca];
    polrat=[polrat polratio];
    
    ' cont diss', keyboard
 %  disp('polcav'), pausak
   end

%%%%%%%%%%%%%%%%%%%%%%%%%

   if exist('Kvet')
    KKdu=Kve(:,fiAz);
    KKdu(1)=1e-10;
    fiK=find(KKdu~=0);
    KK=KKdu(fiK);
    npk=length(KK);
   end
%   pbA=[1:dnum*length(KK)];
%   An0c=An0(pbA);

   if icontr>=1
    figure, plot(abs(An0),'w.'),
    pausak
   end
  end  %gg0
 end % isacc
end  % if ifnz

end % nso

%'qui cont', keyboard
avv=aou;
Gvv=gou;
F=fou;
if exist('isofg')==0
 isofg=1;
end
gsoM=gso.*M2v;
if isofg==1
[du,iso]=sort(fso);
else
[du,iso]=sort(gsoM);
end
%iso=iso(1:nmasce);

%'qui nmaze', keyboard
gso=gso(iso);
if length(gso)>nmasce
 gso(nmasce+1:end)=-gso(nmasce+1:end);
end 

if exist('Gsov')
%'Gsov', keyboard
%if length(Gsov)>length(gso)
for ksov=1:length(Gsov)
 if length(find(Gsov(ksov)/2==-gso))==1 
  Gsov(ksov)=-1;
 end
end
end
%end
fso=fso(iso);

%fso=fso(1:end-1);
%gso=gso(1:end-1);

%fou=fou(iso,:);
%aou=aou(iso,:);
%gou=gou(iso,:);

ns=ns(iso);
ipu=ipu(iso);
aso=aso(iso);


if icampi==0 | length(icampi)==0
 polcav=zeros(size(gso));
 polrat=polcav;
 rtetmv=polcav;
 maziv=polcav;
 mradv=polcav;
end
polcav=polcav(iso);
polrat=polrat(iso);
rtetmv=rtetmv(iso);
maziv=maziv(iso);
mradv=mradv(iso);
%gamso=gamso(iso);
%gamoso=gamoso(iso);
tso=tso(iso);
Am=Anso(:,iso);

nsv(1:length(fso),its)=ns';
ipuv(1:length(fso),its)=ipu';
fsov(1:length(fso),its)=fso';
gsov(1:length(fso),its)=gso';
rtetmsov(1:length(fso),its)=rtetmv';
mazisov(1:length(fso),its)=maziv';
mradsov(1:length(fso),its)=mradv';
polasov(1:length(fso),its)=polcav';
polratv(1:length(fso),its)=polrat';
%gamsov(1:length(fso),its)=gamso';
%gamosov(1:length(fso),its)=gamoso';
asov(1:length(fso),its)=aso';
tsov(1:length(fso),its)=tso';
sm=size(Am);
Amv(its,1:sm(1),1:sm(2))=Am;

 %pv0=avv.*Gvv/vg;
 %z0=aso.*gso/vg;
 pv0=avv;
 fie1=find(tso>=.5);
 z0e=aso(fie1);
 fim1=find(tso<.5);
 z0m=aso(fim1);
 sG=size(Gvv);
end

if isadiss==1
 save diss
end 
%' fine diss_nst', keyboard
