function nlap=flapu(fu,xul);

global pif2 nri Mdelf Nco Nri Mdxr Mdx2 xa xau
global imeto pufi pufu Npdif idiffu

  if imeto==0
    Mpd=fu(pif1,:);
    MDer1d=diff(diff(Mpd'))'.*Mdx2;
    MDer1=[MDer1d MDer1d(:,Nco)];
    MDer2=[diff(Mpd')'.*Mdxr];
    MDerd=MDer2+MDer1;
    MDer=[MDerd(:,1) MDerd];

  elseif imeto==1

    for jf=1:Nri
     F=fu(jf,:);
     cf=polyfit(xul(pufi),F(pufi),2);
     Fav(jf,1)=polyval(cf,xa);
     if idiffu==1
      cf=polyfit(xul(pufu),log(F(pufu)),2);
      Fa=exp(polyval(cf,xau));
     else
      Fa=F(Nri);
     end
     Fauv(jf,1)=Fa;
    end

    Fd=[fu(:,1) Fav fu(:,2:Npdif)];
    dp=diff(Fd')';
    cDer1=dp.*Mdxr;

    pd=[fu Fauv];
    dp1=diff(pd')';
    dp1=[Fauv*0 dp1];

    cDer2=diff(dp1')'.*Mdx2;
    MDer=cDer2+cDer1;
%    disp(' flapfau ')
%    keyboard

  end



 Mpd=fu(pif2,nri);
 MDer2=diff(diff(Mpd));
 MDer0=MDer2.*Mdelf;
 MDef=[MDer0(:,1) MDer0];

 nlap=MDer+MDef;
% keyboard
