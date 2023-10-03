% Ex=real(E2x/(max(max(E2x))));
% Ey=real(E2y/(max(max(E2y))));

 Ex=real(Ex);
 Ey=real(Ey);
 Exm=max(max(Ex));
 Eym=max(max(Ey));
 if Exm/Eym>1
  tyl=1;
 else
  tyl=-1;
 end
 sE=size(Ex);
 pasnu=1;
 nubesu=2;
 lfil=sE(2);
 icapol=1;


 aitot=D3(1,1);
 if isnan(aitot)==1
  aitot=4;
 end
 xvero=xro;
 x=xro;
 ive=1;
 fian0=fian;
fixx=find(xvero<1.2*max(aitot) &  xvero>0.5*max(aitot));
fixxdu=find(xvero<1.2*max(aitot));
fixx0=fixxdu(length(fixxdu));
mxy=1;



if icapol==1
 if i2D==3
   mX=max(max(Ex));
   mY=max(max(Ey));
   if iLP==0
    polratio=(mX/mY)^2;
   else
    polratio=0;
   end
   [du,ifia90]=min(abs(fian-pi/2));
   if mX>mY
    polca=1;
    mX0=sum(Ex(:,1));
    [mX90,i90]=max(sum(Ex(:,2:lfil)));
    ifia90=i90;
     ya=real(Ex(fixx0,:)/mxy);
     if mX0<mX90
      imacchie=-1;
      yy=Ex(:,ifia90)/mxy;
     else
      imacchie=1;
      yy=Ex(:,1)/mxy;
     end
   else
    polca=2;
    mX0=sum(Ey(:,1));
    [mX90,i90]=max(sum(Ey(:,2:lfil)));

    ya=real(Ey(fixx0,:)/mxy);
    ifia90=i90;
     if mX0<mX90
      imacchie=-1;
      yy=Ey(:,ifia90)/mxy;
     else
      imacchie=1;
      yy=Ey(:,1)/mxy;
     end
   end
 else
   yy=Ex;
 end

   yys=real(yy);
   yy=yys(fixx);
   xx=x(fixx)';

   yp=yy;
   lyp=length(yp);
   yp0=yp(2:lyp).*yp(1:lyp-1);
   fim=find(yp0<0);


   mrad=(length(fim));

 if i2D==3
   xx=fian0';
   xx=xx(1:2);
   yas=ya;
   yami=min(ya);
   yama=max(ya);
   ray=yami/yama;

   if ray<0.1
    ive=0;
    if ive==1
     ya=[ya(2) ya];
     yp=diff(ya')./diff(xx);
     fi0=find(abs(yp)<=1e-8);
     if length(fi0)>0
      yp(fi0)=-1e-3;
     end
     lyp=length(yp);
     yp0=[ 1e-3; yp(2:lyp).*yp(1:lyp-1); 1e-3];

     ys=[diff([0; yp])./diff(xx); 0];
     fim=find(yp0<=0 & ys<0);
     fim1=find(yp0<=0 & ys>0);
     nuazi=length(fim)/2;
    else
     indaz=[0:pasnu:2*(nubesu)];
%     MUV=[0:pasnu:2*(nubesu)];
%     if is_even(mm_ver)==0
%      indaz=MUV;
%     else
%      indaz=MUV+1;
%     end
     clear miazs miazc
     for in=1:length(indaz)
      inar=indaz(in);
      miazs(in)=mean(sin(fian*inar).*yas);
      miazc(in)=mean(cos(fian*inar).*yas);
     end
     miaz=abs([miazc miazs]);
     [du,ima]=max(miaz);
     if ima>length(indaz)
      nuazi=-indaz(ima-length(indaz));
     else
      nuazi=indaz(ima);
     end
    end
   else
    nuazi=0;
   end

   if ive==1
    if nuazi~=0
     if nuazi==1 & ya(fim1(1))>0
      if ya(fim(1))/ya(fim1(1))<20
       nuazi=0;
      end
     end
     if fim(1)~=2
      polca=-polca;
     end
    end
   end
  else  %i2D
   nuazi=mm;
  end

  if ifp>-4
   disp('[nuazi mrad]'), disp([nuazi mrad])
  end



else
  polca=0;
  nuazi=0;
  mrad=0;
  rtetm=-1;
  polratio=0;
end   %icapol
ord(2,ip)=nuazi;
ord(1,ip)=mrad;
ty(ip)=tyl;

%' qui mod_chap ', keyboard
