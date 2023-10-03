  if i2D==3

    sumfield
    E2xd=real(Exdu*ensx);
    if iLP==0
     E2yd=real(Eydu*ensy);
     if iztm==1
      E2zd=real(Ezdu*ensz);
     end
    else
     E2yd=0;
     E2zd=0;
    end
    Etot=abs(Exdu).^2+abs(Eydu).^2+abs(Ezdu).^2;
    E=sqrt(Etot);
  end
