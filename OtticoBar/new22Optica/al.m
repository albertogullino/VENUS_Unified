

 h=figure;
 assi=15;
 for k=1:length(zcut)
  [du,ip]=min(abs(zcut(k)-ztot));


 Aue=Azue(:,ip);
 Auez=Azuez(:,ip);
 Auh=Azuh(:,ip);
 Auhz=Azuhz(:,ip);

  if iLP==1

   if exist('if90')
    bes_an=diag(Mvefm(:,if90))*besm;
   else
    bes_an=besm;
   end
   Esr=(bes_an'*Azu);
   Es=abs(Esr/max(max(Esr))).^2;

  else  %vettoriale

   Mfme=diag(Aue)*[ Mvefm;  segem*Mvefm];
   Mgme=diag(Aue)*[ Mvegm;  segem*Mvegm];
   Mfpe=diag(Aue)*[ Mvefp; -segem*Mvefp];
   Mgpe=diag(Aue)*[-Mvegp;  segem*Mvegp];

   Mfmh=diag(Auh)*[ Mvefm;  segem*Mvefm];
   Mgmh=diag(Auh)*[ Mvegm;  segem*Mvegm];
   Mfph=diag(Auh)*[ Mvefp; -segem*Mvefp];
   Mgph=diag(Auh)*[-Mvegp;  segem*Mvegp];

   Mez =diag(Auez)*[ 0*Mvez;         Mvez];
   Mhz =diag(Auhz)*[  Mvhz;        0*Mvhz];

   Ex=besp'*Mfpe+besm'*Mfme;
   Ey=besp'*Mgpe+besm'*Mgme;

   Hy=besp'*Mfph+besm'*Mfmh;
   Hx=-(besp'*Mgph+besm'*Mgmh);


   if iztm==1
    Ez=besz'*Mez;
    Hz=besz'*Mhz;
   else
    Ez=0;
    Hz=0;
   end
   Pvectz=(Ex.*conj(Hy)-Ey.*conj(Hx));
   Pvectx=(Ey.*conj(Hz)-Ez.*conj(Hy));
   Pvecty=-((Ex.*conj(Hz)-Ez.*conj(Hx)));
   Pointingr=sqrt(real(Pvectx).^2+real(Pvecty).^2+real(Pvectz).^2);
   Pointingi=sqrt(imag(Pvectx).^2+imag(Pvecty).^2+imag(Pvectz).^2);

   if ipoi>=1
    if ipoi==2
     Edu=Pointingi;
    else
     Edu=Pointingr;
    end

   else
    if iex==1
     if ieh==1
      Edu=Ex;
     else
      Edu=Hx;
     end
    elseif iex==2
     if ieh==1
      Edu=Ey;
     else
      Edu=Hy;
     end
    elseif iex==3
     if ieh==1
      Edu=Ez;
     else
      Edu=Hz;
     end
    end

    Edu=(Edu/max(max(Edu)));
   end
%
  end

  Esez{k}=Edu;
  map(abs(Edu),YP,XP,h,-1)
%  axis([-1 1 -1 1]*assi(k))
  title(' real part' )
  axis([-1 1 -1 1]*assi)
  axis equal
 if ipoi==0 & length(zcut)<4
  pausak
  map(imag(Edu),YP,XP)
  title(' imag part' )
%  axis([-1 1 -1 1]*assi(k))
  axis([-1 1 -1 1]*assi)
  axis equal
  pausak
 end
  if k==1
  title([' Section = ',num2str(-zcut(k)),' micron'])
  else
  title([num2str(-zcut(k))])
  end
  pausak
 end

