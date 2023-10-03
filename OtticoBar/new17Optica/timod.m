function [mrad,nuazi]=timod(Ex,aax,fian,xvero)

nubesu=2;
pasnu=1;
 sE=size(Ex);
 lfil=sE(2);
x=xvero;
mxy=1;
fian0=fian;

fixx=find(xvero<1.1*aax);
fixx0=fixx(length(fixx));


   mX=max(max(Ex));
   Ex=Ex/mX;
   [du,ifia90]=min(abs(fian-pi/2));
    polca=1;
    su=sum(Ex,1);
    [ma,ima]=max(su);
    [mi,imi]=min(su);
    if ima>imi
     ifia90=ima;
     mX90=ma;
     mX0=mi;
     sgn=-1;
    else
     ifia90=imi;
     mX90=mi;
     mX0=ma;
     sgn=1;
    end
    
     ya=real(Ex(fixx0,:)/mxy);
     if mX0<mX90
      imacchie=-1;
      yy=Ex(:,ifia90)/mxy;
     else
      imacchie=1;
      yy=Ex(:,1)/mxy;
     end

   yys=real(yy);
   yy=yys(fixx);
   xx=x(fixx)';
   yp=yy;
   lyp=length(yp);
   yp0=yp(2:lyp).*yp(1:lyp-1);
   fim=find(yp0<0);


   mrad=(length(fim));


   xx=fian0';
   xx=xx(1:2);
   if abs(mX0/mX90-1)<.1
    yas=ya;
   else
    yas=ya-mean(ya);
   end
   yami=min(ya);
   yama=max(ya);
   ray=yami/yama;
ray=0;
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
     MUV=[0:pasnu:(nubesu)];
     indaz=MUV;
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

nuazi=sgn*nuazi/2;

%if ifp==-100
%'[mrad nuazi] '
%[mrad nuazi]
%O' qui mod_cha ', keyboard
%end
