
 for k=1:2:length(zcut)
%  subplot(3,3,k)
  [du,ip]=min(abs(zcut(k)-ztot));

 fatEz=fatEzv(ip);
 rfip=Nztot(ip);
 beta0=2*pi/lambdam*rfip;
 sqr=sqrt(1-(KKt(Pusas)*rr/rfip).^2);
 lim=1-(KKv*rr/rfip).^2;
 kacce=find(lim<0);
 kaccea=[kacce; kacce+length(lim)];

 if iLP==0
  SQ=sqrt(lim);
  Znorn=[1./SQ; SQ]./Nztot(ip);
  Ynorn=1./Znorn;
  Ynorn(kaccea)=0;
  Ynorn=Ynorn(Pusff);

  fkttn=KKt./[SQ; SQ];
  fkttn(kaccea)=0;
  fkttn=fkttn(Pusff);
 else
  fkttn=KKt(Pusff);
  Ynorn=Nztot(ip);
 end

% fatcut=ones(size(KKt));
% fatcut(kaccea)=0;
% fatcut=fatcut(Pusff);


 Aue=Azue(:,ip);
 Auez=Azuez(:,ip);
 Auh=Azuh(:,ip);
 Auhz=Azuhz(:,ip);

  if iLP==1

   Mfme=diag(Aue)*Mvefm0;
   Esr=besm'*Mfme;
   Edu=abs(Esr/max(max(Esr))).^2;

  else  %vettoriale

   Mfme=diag(Aue)*[ Mvefm0;  segem*Mvefm0];
   Mgme=diag(Aue)*[ Mvegm0;  segem*Mvegm0];
   Mfpe=diag(Aue)*[ Mvefp0; -segem*Mvefp0];
   Mgpe=diag(Aue)*[-Mvegp0;  segem*Mvegp0];

   Mfmh=diag(Auh)*[ Mvefm0;  segem*Mvefm0];
   Mgmh=diag(Auh)*[ Mvegm0;  segem*Mvegm0];
   Mfph=diag(Auh)*[ Mvefp0; -segem*Mvefp0];
   Mgph=diag(Auh)*[-Mvegp0;  segem*Mvegp0];

   Mez =diag(Auez)*[ 0*Mvez0;         Mvez0];
   Mhz =diag(Auhz)*[  Mvhz0;        0*Mvhz0];

%   Mvegp=diag(Az)*Mgp;
%   Ex=besp'*Mvefp+besm'*Mvefm;
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