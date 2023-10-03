 while istr<=length(Li)

   if ifp>=1
    disp(' while sopra: istr'), [istr]
   end
  icop=fmlsi(istr,1);
  ncop=fmlsi(istr,2);
  istr0=istr;
  if itutmir==1
   if ncop~=0
    Ncop=ncop;
   else
    Ncop=1;
   end
  else
   Ncop=1;
  end
  iziv=1;
  ncopF=1;
  istrctv=0;
  ncopco=0;

  ifplaco=1;

  for istrct=1:Ncop

   istr=istr0;
   if ncopco==0
    Tc=Ide;
   end
   ncopco=ncopco+1;
   for istrc=1:icop
    if ifp>=1
     disp(' for: istr, istrc'), [istr istrc], pausak
    end
    if igau==4
     isomv=isomv+Li(istr);
     [du,izi]=min(abs(isomv-zedis(2:length(zedis))));
      KTe=reshape(KTemp(:,:,izi),si);
      KTez=reshape(KTempz(:,:,izi),si);
    end
    if icaco==1

     eval(emme)

     if istrc==1
      Tcm=Tc;
     end
    end
    istr=istr+1;
   end  %istrc
   if igau==4 & itutmir==1
    if Ncop>1
     if izi~=iziv | istrct==Ncop
      icaco=1;
     else
      icaco=0;
     end
     if icaco==1
      ncopF=ncopco;
      ncopco=0;
     end
    else
     icaco=1;
     ncopco=0;
    end
    iziv=izi;
    istrctv=istrct;
   else
    icaco=1;
    ncopF=ncop;
   end

   if icaco==1

    if ncop<=1
     if imem==0
      T=Tc*T;
     else
      T=prodmat(Tc,T,ifplatot);
     end
    else

     if ncopF-fix(ncopF)==0
      if imem==0
       Tco=Tc^ncopF;
      else
       Tco=powermat(Tc,ncopF,ifplaco);
      end
%       ncopF
     else
      if imem==0
       Tco=Tcm*Tc^fix(ncopF);
      else
       Pow=powermat(Tc,fix(ncopF),ifplaco);
       Tco=prodmat(Tcm,Pow,ifplaco);
      end
     end
     if ncop>1

      if imem==0
       Tmir=Tco*Tmir;
      else
       Tmir=prodmat(Tco,Tmir,ifplatot);
      end
     end
     if imem==0
      T=Tco*T;
     else
      T=prodmat(Tco,T,ifplatot);
     end
    end
   end  %icaco
%   disp(' istrct '), istrct
%   pausak
  end  %istrct

  if length(find(istfie==istr))==1
   icmem=icmem+1;
   Tmefb(:,:,icmem)=T;
%   disp(' memorizzo  Tmefb')
%   keyboard
  end
 end   %istr
