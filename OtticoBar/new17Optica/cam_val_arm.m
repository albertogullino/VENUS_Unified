iskim=0;
  if i2D==3



    sumfield
%' dopo sufue', keyboard


%    E2xd=abs(Exdu.^2);
%    E2yd=abs(Eydu.^2);
%    E2zd=abs(Ezdu.^2);
    ENx=max(max(Exdu));
    ENy=max(max(Eydu));
    %'iskim', keyboard
    if abs(ENx)+abs(ENy)==0
     iskim=1;
     return
    end
if imazver==1
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
    if kmo==1 
     Ensx=ensx;
     Ensy=ensy;
     Ensz=ensz;
    end
else   
     ensx=Ensx;
     ensy=Ensy;
     ensz=Ensz;
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


%'fine camval arm', keyboard