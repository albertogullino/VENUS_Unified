% Ex=real(E2x/(max(max(E2x))));
% Ey=real(E2y/(max(max(E2y))));

 Ex=real(E2x);
 Ey=real(E2y);
 An=Andu;
 sE=size(Ex);
 lfil=sE(2);


%if ifalso>=0
 aax=1.2*max([avero; bvero]);
 aarad=aax/1.2;
%else
% aax=1.1*rax;
% aarad=aax/1.1;
%end
fixx=find(xvero<1.1*max(max(aitot)));
fixx0=fixx(length(fixx));
mxy=1;


if i2D==3
   mX=max(max(Ex));
   mY=max(max(Ey));
   if iLP==0
    polratio=(mX/mY)^2;
   else
    polratio=0;
   end
%' prima ya'
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
%   yp=diff(yy)./diff(xx);
%   lyp=length(yp);
%   yp0=[ 0; yp(2:lyp).*yp(1:lyp-1); 0];
%   ys=[diff([0; yp])./diff(xx); 0];
%   fim0=find(yp0<0 & yy>.01 & ys<0);
%   fim=find(yp0<=0 & yy>.01 & ys<0);
%   fim=find(yp0<0 & abs(yy)>.01);
%   if yy(1)>.02
%    fim=[1; fim];
%   end
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
     MUV=[0:pasnu:2*(nubesu)];
     if is_even(mm_ver+iLP)==0
      indaz=MUV;
     else
      indaz=MUV+1;
     end
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

if icapol==1

   if iLP==0
    l=fix(length(An)/2);
    le=1:l;
    lm=l+1:l*2;
    A=abs(An);
    Ate=A(le);
    Atm=A(lm);
    se=sum(Ate);
    sm=sum(Atm);
    rtetm=se/(se+sm);
   else
    rtetm=-1;
   end

else
  polca=0;
%  nuazi=0;
%  mrad=0;
  rtetm=-1;
  polratio=0;
end   %icapol

if ifp==-100
'[mrad nuazi] '
[mrad nuazi]
' qui mod_cha ', keyboard
end
clear ya
%'dopo mod_cha', keyboard