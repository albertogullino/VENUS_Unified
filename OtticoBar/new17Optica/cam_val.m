iskim=0;
  if i2D==3



    sumfield
%' dopo sufue', keyboard


%    E2xd=abs(Exdu.^2);
%    E2yd=abs(Eydu.^2);
%    E2zd=abs(Ezdu.^2);
    ENx=max(max(Exdu));
    ENy=max(max(Eydu));

    if abs(ENx)+abs(ENy)==0
     iskim=1;
     'campo nullo (pol. errata ed iLP=1)!', keyboard
     Exdu=zeros(length(xvero),length(fian));
     Eydu=Exdu;
     Ezdu=Exdu;
     Hzdu=Exdu;
     Hxdu=Exdu;
     Hydu=Exdu;
     Etot=abs(Exdu).^2+abs(Eydu).^2;
     E=sqrt(Etot);
%       'iskim', keyboard
     return
    end

    if abs(ENx)>abs(ENy)
     ensx=1/ENx;
     if ENy~=0
      ensy=abs(ENy/ENx)/ENy;
     end
     if iztm==1
      ENz=max(max(Ezdu));
      ensz=abs(ENz/ENx)/ENz;
     end
    else
     ensy=1/ENy;
     ensx=abs(ENx/ENy)/ENx;
     if iztm==1
      ENz=max(max(Ezdu));
      ensz=abs(ENz/ENy)/ENz;
     end
    end
%    E2xd=real(Exdu*ensx);
    E2xd=(Exdu*ensx);
    if iLP==0
%     E2yd=real(Eydu*ensy);
     E2yd=(Eydu*ensy);
     if iztm==1
%      E2zd=real(Ezdu*ensz);
      E2zd=(Ezdu*ensz);
     end
    else
     E2yd=0;
     E2zd=0;
    end
    Etot=abs(Exdu).^2+abs(Eydu).^2+abs(Ezdu).^2;
    Etot=abs(Exdu).^2+abs(Eydu).^2;
    E=sqrt(Etot);
  else
   AnFs=Andu.*Fs;
   Edu=AnFs(puAc).'*besm(puAc,:);
   Etot=abs(Edu.^2)';
   E=sqrt(Etot);
   EN=max(max(Edu));
%   E2xd=real(Edu/EN);
   E2xd=(Edu/EN);
   E2yd=zeros(size(Etot));
  end
  mxy=max(max(Etot));
  %'Etot', keyboard
  if mxy==0
   iskim=1;
  end

   if kmo==20
    Andus=Andu;
    lAn=length(Andu);
    Andu(lAn/2+1:end)=0;
    sumfield
    pograd=[ 40 400 800 300];
    zete=0;
    ducam

    Andu=Andus;
    Andu(1:lAn/2)=0;
    sumfield
    pograd=[ 40 70 800 300];
    zete=1;
    ducam
    Andu=Andus;

   end


%'fine camval', keyboard