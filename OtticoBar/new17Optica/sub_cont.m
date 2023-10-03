
 Andu=Andu0;
 cam_val
 E2x=E2xd;
 E2y=E2yd;

 figure,
 if iLP==0
  pograp=[200   50   600   900];
  set(gcf,'Position',pograp);
  subplot(2,1,1)
 end
 ibar=1;
 iaoff=0;
   titl=' E_x  Output';
  map_fnew(XP,YP,E2x,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
  if iLP==0
   titl=' E_y  Output';
   subplot(2,1,2)
   map_fnew(XP,YP,E2y,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
  end
  pausak
  nmodi=length(ibm);

    pa=2;
    ngy=3;
    ngx=ceil(nmodi/ngy/pa);
    vsel0=zeros(size(Andu0));
    vk=1:pa*length(KK);

    figx=figure,
    pograp=[50   50   300*ngy   300*ngx];
    set(gcf,'Position',pograp);
    if iLP==0
     figy=figure,
     pograp=[150   50   300*ngy   300*ngx];
     set(gcf,'Position',pograp);
    end

   for imc=1:pa:nmodi
    pusel0=(imc-1)*length(KK)+vk;

    if iLP==0
     pusel=[pusel0 pusel0+nmodi*length(KK)];
    else
     pusel=pusel0;
    end
    vsel=vsel0;
    vsel(pusel)=1;
    Andu=Andu0.*vsel;
%    figure, plot(abs(Andu)), pausak
    cam_valc
    E2x=E2xd;
    E2y=E2yd;
    ibar=1;
    iaoff=0;
      titl=[' E_x  Output for mode  ',num2str(ibm(imc))];
     figure(figx)
     subplot(ngx,ngy,ceil(imc/pa))
     map_fnew(XP,YP,E2x,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
     if iLP==0
      titl=[' E_y  Output for mode  ',num2str(ibm(imc))];
      figure(figy)
      subplot(ngx,ngy,ceil(imc/pa))
      map_fnew(XP,YP,E2y,aax,Cug.x,Cug.y,Cug.z,titl,ibar,iaoff)
     end
   end
